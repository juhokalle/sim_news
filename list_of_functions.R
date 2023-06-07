# SIMULATIONS
simu_y = function(model, n.obs, rand.gen = stats::rnorm, n.burnin = 0, ...)
{
  
  d = dim(model$sys)
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  
  if ( (m*(p+1)) == 0 ) stop('illegal ARMA system (m=0 or p<0)')
  
  # check 'n.obs'
  n.obs = as.integer(n.obs[1])
  if (n.obs <= 0) stop('"n.obs" must be a positive integer!')
  
  # check 'n.burnin'
  n.burnin = as.integer(n.burnin[1])
  if (n.burnin < 0) stop('"n.burnin" must be a non negative integer!')
  
  # get ARMA parameters
  a = unclass(model$sys$a)
  b = unclass(model$sys$b)
  
  # convert ARMA system
  #    a[0] y[t] + a[1] y[t-1] + ... = b[0] u[t] + b[1] u[t-1] + ...
  # to system of the form
  #    y[t] = a[1] y[t-1] + ... + b[0] u[t] ...
  # and create parameter matrices a = (a[p], ..., a[1]) and b = (b[0], ..., b[q])
  a0 = matrix(a[,,1], nrow = m, ncol = m)
  
  dim(b) = c(m, n*(q+1))
  if ((n*(q+1)) > 0) {
    b = solve(a0, b)
  }
  
  # note for solve_ARMA_cpp we have to reshuffle the AR parameter as:  a = (a[p],...,a[1])
  if (p > 0) {
    a = a[,,(p+1):2, drop = FALSE]
    dim(a) = c(m, m*p)
    a = -solve(a0, a)
  } else {
    a = matrix(0, nrow = m, ncol = 0)
  }
  
  # generate disturbances
  u = matrix(rand.gen((n.obs+n.burnin)*n), nrow = n, ncol = n.obs+n.burnin)
  u = model$sigma_L %*% t(solve(chol(var(t(u[, (n.burnin+1):(n.burnin+n.obs)]))))) %*% u
  # outputs
  y = matrix(0, nrow = m, ncol = n.obs+n.burnin)
  svarmawhf::solve_ARMA_cpp(a, b, u, y, 1)
  
  return(list(y = t(y[, (n.burnin+1):(n.burnin+n.obs), drop = FALSE]),
              u = t(u[, (n.burnin+1):(n.burnin+n.obs), drop = FALSE])))
}

sim_news <- function(beta, rho, 
                     no_sim = FALSE, nobs = 250, nu = 3, 
                     rg_fun = function(x) stats::rt(x, nu))
{
  theta <- 1/(1-beta*rho)
  ar_pol <- array(c(diag(2),                   # lag 0
                    -rho, -rho*theta, 0, 0),   # lag 1
                  dim = c(2,2,2))
  ma_pol <- array(c(diag(2), # lag 0
                    0, -theta/beta, 0, 1/beta,       # lag 1
                    -1/beta^2, -theta/beta^2, 1/(theta*beta^2), 1/beta^2), # lag 2
                  dim = c(2,2,3))
  bmat <- matrix(c(1, theta, 0, theta*beta^2), 2, 2)
  re_mod <- armamod(sys = lmfd(ar_pol, ma_pol), sigma_L = bmat)
  if(no_sim) return(re_mod)
  data_out <- simu_y(model = re_mod,
                     rand.gen =  rg_fun,
                     n.burnin = 1000,
                     n.obs = nobs)
  return(list(y = data_out, 
              mod = re_mod, 
              prms = c("beta" = beta, "rho" = rho, "nobs" = nobs, "nu" = nu))
         )
}

mixed_marg_dists <- function(dim_out, nu)
{
  function(x) c(replicate(x/dim_out, c(rnorm(2), stats::rt(dim_out-2, nu))))
}

get_simu_model <- function(mod_tmpl, imp_mat)
{
  root_flag <- TRUE
  scl_prm <- 1
  n_unst <- k_kappa2nunst(q = mod_tmpl$input_orders$MAorder,
                          dim_out = mod_tmpl$input_orders$dim_out,
                          k = mod_tmpl$input_orders$k,
                          kappa = mod_tmpl$input_orders$kappa)
  
  while(root_flag && scl_prm > 0){
    mod_prms <- rnorm(mod_tmpl$n_par, mean = 0, sd = scl_prm)
    arma_whf <- fill_tmpl_whf_rev(mod_prms, mod_tmpl)
    arma_std <- armamod_whf(mod_prms, mod_tmpl)
    Bmat <- diag(diag(imp_mat))
    ma_polm <- arma_std$polm_ma%r%(solve(unclass(arma_std$polm_ma)[,,1])%*%imp_mat%*%solve(Bmat))
    # Check root constraints
    ar_flag <- any(abs(eigen(companion_matrix(arma_whf$polm_ar))$val) > 1)
    ma_bwd_flag <- any(abs(eigen(companion_matrix(arma_whf$polm_ma_bwd))$val) > 1)
    ma_fwd_flag <- any(abs(eigen(companion_matrix(arma_whf$polm_ma_fwd))$val) > 1)
    ma_flag <- !(sum(abs(eigen(companion_matrix(ma_polm))$val)>1)==n_unst)
    root_flag <- ar_flag || ma_bwd_flag || ma_fwd_flag || ma_flag
    scl_prm <- scl_prm-.05
  }
  armamod(lmfd(arma_std$polm_ar, ma_polm), Bmat)
}

# BOOTSTRAP
mb_boot <- function(y, prms, tmpl, b.length=10, nboot=500)
{
  # Most of this is a direct copy of svars::mb.boot
  u <- get_residuals_once(params_deep = prms, tmpl = tmpl, data_long = y)
  armamod <- armamod_whf(theta = prms, tmpl = tmpl)
  obs <- nrow(y)
  k <- ncol(y)
  errors <- list()
  N <- obs/b.length
  blocks <- array(NA, c(b.length, k, obs - b.length + 1))
  
  for (i in 0:(obs - b.length)) {
    blocks[, , (i + 1)] <- u[(i + 1):(i + b.length), ]
  }
  for (i in 1:nboot) {
    epsilon.star <- matrix(0, b.length * (ceiling(N) + 1), 
                           ncol(u))
    epsilon.star <- list()
    for (kk in 1:(ceiling(N) + 1)) {
      epsilon.star[[kk]] <- blocks[, , floor(runif(1, 1, obs - b.length + 2))]
    }
    epsilon.star <- do.call("rbind", epsilon.star)
    for (s in 1:b.length) {
      b.mean <- colSums(epsilon.star[1:(s + (obs - b.length)), ])/(obs - b.length + 1)
      for (j in 0:floor(N)) {
        epsilon.star[j * b.length + s, ] <- epsilon.star[j * b.length + s, ] - b.mean
      }
    }
    epsilon.star <- epsilon.star[1:obs, ]
    errors[[i]] <- epsilon.star
  }
  ystar <- lapply(errors, function(x) solve_de(lmfd(armamod$polm_ar, armamod$polm_ma), x)$y)
  return(ystar)
}

# SHOCK ID & LABELLING
get_llf <- function(p, q, kappa, k, dtype, sd)
{
  total_data <- readRDS("local_data/total_data.rds")
  data_i <- total_data %>% 
    filter(mp_type==dtype) %>% 
    pull(data_list) %>% .[[1]]
  dim_out <- dim(data_i)[2]
  tmpl_i <- tmpl_whf_rev(dim_out = dim_out,
                         ARorder = p,
                         MAorder = q,
                         kappa = kappa,
                         k = k,
                         shock_distr = "sgt")
  ll_whf_factory(data_wide = t(data_i), tmpl = tmpl_i, shock_distr = sd)
  
}

id_news_shox <- function(irf_arr, policy_var)
{
  dim_in <- dim(irf_arr)[2]
  news_ix <- which.min(abs(unclass(irf_arr)[policy_var,,1]))
  rot_mat <- diag(dim_in)[, c((1:dim_in)[-news_ix], news_ix)]
  irf_tmp <- irf_arr%r%rot_mat
  # Impose positive responses to the shocks: 1) conventional shock(s)
  diag0 <- unclass(irf_tmp)[1:(dim_in-1), 1:(dim_in-1), 1]
  if(dim_in > 2) diag0 <- diag(diag0)
  rot_mat <- rot_mat%*%diag(c(sign(diag0), 1))
  
  # 2) news shock
  peak_ix <- which.max(abs(unclass(irf_tmp)[policy_var, dim_in, ]))
  if(unclass(irf_tmp)[policy_var, dim_in, peak_ix] < 0){
    rot_mat <- rot_mat%*%diag(c(rep(1, dim_in-1), -1))
  }
  rot_mat
}

rotmat <- function(x, n)
{
  rot_ix <- combn(n,2)
  final_mat <- diag(n)
  
  for(jj in 1:ncol(rot_ix)){
    rot_mat <- matrix(c(cos(x[jj]), -sin(x[jj]), sin(x[jj]), cos(x[jj])), 2, 2)
    temp_mat <- diag(n)
    temp_mat[rot_ix[,jj], rot_ix[,jj]] <- rot_mat
    final_mat <- final_mat %*% temp_mat
  }
  final_mat
}

ff <- function(x, zero_ix, input_mat)
{
  row_ix <- if(is.null(dim(zero_ix))) zero_ix[1] else zero_ix[,1]
  col_ix <- if(is.null(dim(zero_ix))) zero_ix[2] else zero_ix[,2]
  sum(sqrt(diag(input_mat[row_ix,] %*% rotmat(x, ncol(input_mat))[, col_ix])^2))
}

id_policy_shox <- function(irf_arr, policy_var, n_shox)
{
  
  dim_out <- dim(irf_arr)[1]
  n_ahead <- dim(irf_arr)[3]
  fevd_obj <- get_fevd(irf_arr, int_var = policy_var)[[1]]
  combs <- combn(dim_out, n_shox)
  res_mt <- matrix(NA, n_ahead, ncol(combs))
  
  for(i in 1:n_ahead){
    for(j in 1:ncol(combs)){
      res_mt[i, j] <- sum(fevd_obj[i, combs[,j]])
    }
  }
  combs[,which.max(colMeans(res_mt))]
}

optim_zr <- function(input_mat, zr_ix, opt_it = TRUE)
{
  n_var <- dim(input_mat)[1]
  
  if(opt_it){
    # chooses the zero restriction based on minimum distance (measured by F-rnom) 
    # from unconstrained estimator
    opt_ix <- rep(0, n_var)
    rotmat_ix <- array(NA, dim = c(n_var, n_var, n_var))
    n_rest <- length(zr_ix)
    for(jj in 1:n_var){
      # UNCOSTRAINED OPTIM: BFGS
      zero_ix = if(n_rest==1) c(zr_ix, jj) else matrix(c(zr_ix, rep(jj, n_rest)), n_rest, n_rest)
      optim(par = rep(0, (n_var^2-n_var)/2),
            fn = ff, 
            zero_ix = zero_ix,
            input_mat = input_mat,
            method = "BFGS") -> opt_obj
      
      rotmat_ix[,,jj] <- rotmat(opt_obj$par, n_var)
      opt_ix[jj] <- norm(input_mat - input_mat%*%rotmat_ix[,,jj], "F")
      return(rotmat_ix[,,which.min(opt_ix)])
    }
  } else {
    # ex ante fixed the zero restriction
    optim(par = rep(0, (n_var^2-n_var)/2),
          fn = ff, 
          zero_ix = zr_ix,
          input_mat = input_mat,
          method = "BFGS") -> opt_obj
    
    return(rotmat(opt_obj$par, n_var))
  }
}

get_fevd <- function (irf_arr, int_var = NULL, by_arg = NULL)
{
  irf_arr <- unclass(irf_arr)
  nvar <- dim(irf_arr)[1]
  n.ahead <- dim(irf_arr)[3]
  fe <- list()
  if(is.null(int_var)){
    for (i in 1:nvar) {
      fe[[i]] <- as.data.frame(t(irf_arr[i, , ]))
    }
  } else{
    fe[[1]] <- as.data.frame(t(irf_arr[int_var, , ]))
  }
  fe2 <- fe
  seq_grid <- if(is.null(by_arg)) 1:n.ahead else seq(1,n.ahead, by = by_arg)
  for (i in 1:length(fe)) {
    for (j in seq_grid) {
      fe2[[i]][j, ] <- (colSums(fe[[i]][j:1, ]^2)/sum(fe[[i]][j:1, ]^2)) * 100
    }
    fe2[[i]] <- fe2[[i]][seq_grid,]
  }
  return(fe2)
}

norm_irf <- function(irf_arr, 
                     norm_int, 
                     norm_pos=dim(irf_arr)[1], 
                     norm_scm = c("surp", "news"))
{
  
  if(norm_scm=="surp"){
    irf_arr <- irf_arr[,1,,drop=FALSE]/irf_arr[norm_pos,1,1]*norm_int
  } else if(norm_scm=="news"){
    irf_arr <- irf_arr[,1,,drop=FALSE]/max(irf_arr[norm_pos,1,])*norm_int
  }
  irf_arr
}

get_rest_irf <- function(tbl_slice, ...)
{
  nobs <- tbl_slice$nobs
  sd_mat <- readRDS("local_data/svarma_data_list.rds") %>% 
    filter(mp_type%in%tbl_slice$mp_type) %>%
    pull(std_dev) %>% .[[1]] %>%  diag
  
  irf_out <- sd_mat %r% tbl_slice$irf[[1]]
  
  opt_obj <- optim_zr(irf_out[,,1], ...)
  llf0 <- do.call(get_llf, tbl_slice %>% 
                    select(p, q, kappa, k, mp_type, sd) %>% 
                    rename(dtype = mp_type)
  )
  params0 <- tbl_slice$params_deep_final[[1]]
  
  params0[tbl_slice$params_deep_final[[1]] %in% tbl_slice$B_mat[[1]]] <- c(tbl0$B_mat[[1]]%*%opt_obj$rmat)
  pval <- 1-pchisq(2*nobs*(llf0(params0)-tbl0$value_final), opt_obj$n_rest)
  
  irf_out <- irf_whf(params0, tbl_slice$tmpl[[1]], 48)
  
  return(list(irf = irf_out, pval = pval, rmat = opt_obj$rmat))
}

choose_perm_sign <- function(target_mat, cand_mat, type = c("frob", "dg_abs"))
{
  nvar <- ncol(cand_mat)
  sign_ix <- nvar %>% replicate(list(c(-1,1))) %>% expand.grid
  perm_ix <- nvar %>% 
    replicate(list(1:nvar)) %>% 
    expand.grid %>% 
    filter(apply(., 1, n_distinct)==nvar)
  cr0 <- 1e25
  for(j in 1:nrow(sign_ix)){
    for(jj in 1:nrow(perm_ix)){
      rt_mat <- diag(sign_ix[j, ]) %*% diag(nvar)[, unlist(perm_ix[jj, ])]
      x1 <- cand_mat %*% rt_mat
      if(type=="frob"){
        cr1 <- norm(target_mat - x1, "F")
      } else if(type=="dg_abs"){
        cr1 <- if(all(diag(x1)>0)) -abs(prod(diag(x1))) else 1e25
      }
      if(cr1 < cr0){
        x_opt <- x1
        rt_opt <- rt_mat
        cr0 <- cr1
      }
    }
  }
  rt_opt
}

# RESULTS: FIGS & TBLS
plot_irf <- function(irf_arr, var_name, shock_name)
{
  n_var <- dim(irf_arr)[1]
  n_shock <- dim(irf_arr)[2]
  n_ahead <- dim(irf_arr)[3]-1
  tibble(irf = c(irf_arr),
         variable = var_name %>% 
           factor(levels = var_name) %>% 
           rep((n_ahead+1)*n_shock),
         shock = shock_name %>% 
           factor(levels=shock_name) %>% 
           rep(each=n_var) %>% 
           rep(n_ahead+1),
         months = rep(0:n_ahead, each = n_var*n_shock)) %>% 
    ggplot(aes(x=months, y = irf)) +
    geom_line() +
    geom_hline(yintercept = 0, size = 0.15) +
    facet_grid(variable ~ shock, scales = "free_y") +
    facet_rep_grid(variable ~ shock, scales = "free_y", repeat.tick.labels = 'left') +
    scale_x_continuous(breaks = seq(12, n_ahead, by = 12), expand = c(0,0)) +
    scale_y_continuous(n.breaks = 5) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 6, angle = 0),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=6),
          panel.grid.major = element_blank(),
          panel.spacing = unit(.5, "lines"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          axis.line = element_line(colour = "black"))
}

get_stats_tbl <- function(x)
{
  nvar <- ncol(x)
  res <- list()
  for(j in 1:nvar){
    sw_test <- shapiro.test(x[,j])
    jb_test <- tsoutliers::JarqueBera.test(x[,j])
    lb_lvl <- Box.test(x[,j], lag = 24, type = "Ljung-Box")
    lb_abs <- Box.test(abs(x[,j]), lag = 24, type = "Ljung-Box")
    lb_sq <- Box.test(x[,j]^2, lag = 24, type = "Ljung-Box")
    res[[j]] <- c(sw_test$statistic,
                  sw_test$p.value,
                  jb_test[[1]]$statistic,
                  jb_test[[1]]$p.value, 
                  lb_lvl$statistic,
                  lb_lvl$p.value, 
                  lb_abs$statistic, 
                  lb_abs$p.value,
                  lb_sq$statistic,
                  lb_sq$p.value) %>% round(3)
    res[[j]][!1:length(res[[j]])%%2] <- paste0("(", res[[j]][!1:length(res[[j]])%%2], ")")
  }
  res <- reduce(res, cbind)
  rownames(res) <- c("SW-test", "", "JB-test", "",
                     "LB 1", "", 
                     "LB 2", "",
                     "LB 3", "")
  res_ltx <- xtable::xtable(res)
  res_ltx
}

get_cor_tbl <- function(x)
{
  comb_i <- combn(ncol(x),2)
  l_corr <- cor(x, method = "pearson") %>% round(3)
  r_corr <- cor(x, method = "spearman") %>% round(3)
  for(j in 1:ncol(comb_i)){
    
    l_corr[comb_i[2, j], comb_i[1, j]] <- paste0("(", 
                                                 cor.test(x[, comb_i[1,j]], 
                                                          x[, comb_i[2,j]], 
                                                          method = "pearson")$p.value %>% round(3),
                                                 ")")
    r_corr[comb_i[2, j], comb_i[1, j]] <- paste0("(",
                                                 cor.test(x[, comb_i[1,j]], 
                                                          x[, comb_i[2,j]], 
                                                          method = "spearman")$p.value %>% round(3),
                                                 ")")
    
  }
  list(noquote(cbind(l_corr, r_corr)), xtable::xtable(cbind(l_corr, r_corr)))
}

get_qqplots <- function(x)
{
  plot_list <- list()
  for(j in 1:ncol(x)){
    df <- data.frame(y = x[,j])
    p <- ggplot(df, aes(sample = y))
    p <- p + 
      stat_qq() + 
      stat_qq_line() +
      ylab("") +
      xlab("") +
      if(j==1){
        ggtitle(bquote(epsilon[1]))
      } else if (j==2){
        ggtitle(bquote(epsilon[2]))
      } else if (j==3){
        ggtitle(bquote(epsilon[3]))
      } else if (j==4){
        ggtitle(bquote(epsilon[4]))
      }
    plot_list[[j]] <- p + theme(plot.title = element_text(hjust = 0.5))
  }
  plot_list
}

sim_comparison <- function(irf_arr, shock_ix, qntl, prms)
{
  n_ahead <- dim(irf_arr)[3]-1
  res_final <- expand_grid(prms, model = c("VAR", "VARMA"))
  res_final$irf <- t(expand_grid(1:nrow(prms), 1:2)) %>% 
    data.frame() %>% 
    as.list %>% 
    map(~apply(X = irf_arr[shock_ix[1], shock_ix[2],, .x[1],, .x[2]],
               MARGIN = 1,
               FUN = function(x) quantile(x, qntl)))
  res_final <- res_final %>%
    mutate(irf_lb = map(irf, ~.x[1,]),
           irf_ub = map(irf, ~.x[2,])) %>% 
    dplyr::select(-irf) %>%
    unnest_longer(starts_with("irf")) %>% 
    mutate(lag = rep(0:n_ahead, n()/(n_ahead+1)))
  
  res_final <- tibble(prms) %>% 
    mutate(irf_md = map(.x = beta, 
                        ~with(sim_news(beta = .x, rho = 0.5, no_sim = TRUE), 
                              pseries(sys, n_ahead)%r%sigma_L) %>%
                          unclass %>% 
                          .[shock_ix[1],shock_ix[2],]
    )
    ) %>% 
    unnest_longer("irf_md") %>% 
    mutate(lag = rep(0:n_ahead, n()/(n_ahead+1)),
           model = "MOD_0") %>% 
    bind_rows(res_final)
  
  plt_list <- list()
  for(jj in 1:nrow(prms))
  {
    dt_plt <- res_final %>% filter(beta==prms$beta[jj], nu==prms$nu[jj])
    dt_plt <- tibble(irf0 = dt_plt %>% filter(model=="MOD_0") %>% pull(irf_md),
                     varma_lb = dt_plt %>% filter(model=="VARMA") %>% pull(irf_lb),
                     varma_ub = dt_plt %>% filter(model=="VARMA") %>% pull(irf_ub),
                     var_lb = dt_plt %>% filter(model=="VAR") %>% pull(irf_lb),
                     var_ub = dt_plt %>% filter(model=="VAR") %>% pull(irf_ub),
                     lag = dt_plt %>% filter(model=="MOD_0") %>% pull(lag))
    plt_list[[jj]] <- dt_plt %>% ggplot(aes(x=lag)) + 
      geom_line(aes(y=irf0)) + 
      geom_line(aes(y=varma_lb), linetype="dashed", size = 0.25) +
      geom_line(aes(y=varma_ub), linetype="dashed", size = 0.25) +
      geom_line(aes(y=var_lb), linetype="dotted", size = 0.25) +
      geom_line(aes(y=var_ub), linetype="dotted", size = 0.25) +
      geom_ribbon(aes(ymin=varma_lb, ymax = varma_ub), alpha = 0.1, fill="blue") + 
      geom_ribbon(aes(ymin=var_lb, ymax = var_ub), alpha = 0.1, fill = "red") + 
      scale_x_continuous(expand = c(0,0)) +
      ggtitle(bquote(.(paste("(", letters[c(1,3,2,4)][jj], "): ", "")) ~
                       beta == .(paste(prms$beta[jj], ",", sep = "")) ~
                       m == .(prms$nu[jj])
      )
      ) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")) +
      ylab("") + 
      xlab("")
    
  }
  plt_list
}

# OTHER
pmap_tmpl_whf_rev = function(dim_out = DIM_OUT, p, q, kappa, k, shock_distr = "gaussian", ...)
{
  tmpl_whf_rev(dim_out = DIM_OUT, ARorder = p, MAorder = q, kappa = kappa, k = k, shock_distr = shock_distr)
}

pmap_get_residuals_once = function(params_deep_final, tmpl, data_list, ...)
{
  get_residuals_once(params_deep = params_deep_final, tmpl = tmpl, data_long = data_list)
}

hlp_parallel = function(list_input)
{
  return(create_results_list(theta_init = list_input[[1]], 
                             tmpl       = list_input[[2]],
                             params     = params, 
                             DATASET    = list_input[[3]],
                             shock_distr = list_input[[4]]))
}

k_kappa2nunst <- function(q, dim_out, k, kappa){
  
  kk <- c(k, kappa)
  hlp_tbl <- tibble(n_unst = map(q, ~0:(dim_out*.x))) %>% 
    unnest(n_unst) %>% 
    mutate(n_st = dim_out * q - n_unst) %>% 
    mutate(kappa = n_unst %/% dim_out,
           k = n_unst %% dim_out)
  hlp_tbl %>% 
    filter(k == kk[1], kappa == kk[2]) %>% 
    dplyr::select(n_unst) %>% 
    as.double
}

# Tools related to structrual models
get_struc_mat <- function(model_type = c("dynamic", "static", "r_smooth"), param_list = NULL){

  if(model_type == "dynamic"){
    
    if(is.null(param_list)) param_list <- prms <- list(beta = .99, kappa = .05, gamma = .5,
                                                       delta_x = .1, alfa = .5, 
                                                       tau_pi = 1.8, tau_x = .5, tau_r = .6, 
                                                       rho_pi = .5, rho_x = .5, rho_r = .5)
    param_vec <- c("beta", "alfa", "kappa", "gamma", "delta_x", 
                   "tau_x", "tau_pi", "tau_r", 
                   "rho_x", "rho_pi", "rho_r")
    check_vec <- param_vec%in%names(param_list)
    if(!all(check_vec)) stop("The following parameters are missing: ", 
                             paste(param_vec[!check_vec], collapse = ", "))
    
    Kmat <- with(param_list, matrix(c(1, -kappa, -tau_x*(1-tau_r),
                                      0, 1, -tau_pi*(1-tau_r),
                                      delta_x, 0, 1), 3, 3))
    
    Amat <- with(param_list, matrix(c(1-gamma, 0, 0,
                                      0, alfa/(1+alfa*beta), 0,
                                      0, 0, tau_r), 3, 3))
    Bmat <- with(param_list, matrix(c(gamma, 0, 0,
                                      delta_x, beta/(1+alfa*beta), 0,
                                      0, 0, 0), 3, 3))
    Hmat <- diag(3)
    
    Dmat <- with(param_list, diag(c(rho_x, rho_pi, rho_r)))
    
  } else if(model_type == "static"){
    
    if(is.null(param_list)) param_list <- list(beta = 0.995, kappa = 0.005,
                                                phi_pi = 1.5, phi_y = 0.1,
                                                sigma_d = 1.6, sigma_s = 0.95, sigma_m = 0.23)
    param_vec <- c("beta", "kappa",
                   "phi_pi", "phi_y",
                   "sigma_d", "sigma_s", "sigma_m")
    check_vec <- param_vec%in%names(param_list)
    if(!all(check_vec)) stop("The following parameters are missing: ", 
                             paste(param_vec[!check_vec], collapse = ", "))
    
    Kmat <- with(param_list, matrix(c(1, -kappa, -phi_y,
                                      0, 1, -phi_pi,
                                      1, 0, 1), 3, 3))
    
    Amat <- matrix(0, 3, 3)
    
    Bmat <- with(param_list, matrix(c(1, 0, 0,
                                      1, beta, 0,
                                      0, 0, 0), 3, 3))
    
    Hmat <- with(param_list, diag(c(sigma_d, sigma_s, sigma_m)))
    
    Dmat <- matrix(0, 3, 3)
  
  } else if(model_type == "r_smooth"){
    
    if(is.null(param_list)) param_list <- list(beta = .99, kappa = .75, tau = 2.08^-1,
                                               psi_pi = 2.19, psi_x = 0.3, rho_R = 0.84,
                                               sigma_g = .21, sigma_z = 1.16, sigma_R = 0.24)
    
    param_vec <- c("tau", "kappa", "rho_R", "beta", 
                   "psi_x", "psi_pi",
                   "sigma_g", "sigma_z", "sigma_R")
    check_vec <- param_vec%in%names(param_list)
    if(!all(check_vec)) stop("The following parameters are missing: ", 
                             paste(param_vec[!check_vec], collapse = ", "))
    
    Kmat <- with(param_list, matrix(c(1, -kappa, -(1-rho_R)*psi_x,
                                      0, 1, 0,
                                      tau, 0, 1), 3, 3))
    
    Amat <- with(param_list, matrix(c(rep(0,8), rho_R), 3, 3))
    
    Bmat <- with(param_list, matrix(c(1, 0, 0,
                                      tau, beta, (1-rho_R)*psi_pi, 
                                      0, 0, 0), 3, 3))
    
    Hmat <- with(param_list, matrix(c(sigma_g, 0, 0, 
                                      0, -kappa*sigma_z, -(1-rho_R)*psi_x*sigma_z,
                                      0, 0, sigma_R), 3, 3))
    Dmat <- matrix(0, 3, 3)
    
  } else{
    stop("Supply a valid model type.")
  }
  
  list(Kmat=Kmat, Amat=Amat, Bmat=Bmat, Hmat=Hmat, Dmat=Dmat)
}

# The system is of the form
# K * x[t] = A * x[t-1] + B * E(x[t+1]|I[t])
#                + H * w[t]
# w[t] = D * w[t-1] + v[t]

solve_re_mod_bp <- function(Kmat, Amat, Bmat, Hmat, Dmat, eps_val = 1e-9){
  
  # Transform System to Canonical Form:
  # x[t] = Ahat * x[t-1] + Bhat * E(x[t+1]|I[t]) 
  #        + Hhat * w[t]
  Ahat = solve(Kmat, Amat)
  Bhat = solve(Kmat, Bmat)
  Hhat = solve(Kmat, Hmat)
  
  dim1 = dim(Amat)[1]
  dim2 = dim(Amat)[2]
  
  # Compute Matrix C Using Brute-Force Iterative Procedure
  Cmat = diag(dim1)       # Initial Condition
  Fmat = diag(dim1)       # Initial Condition
  iter <- 1
  while(iter == 1 || crit1 >= eps_val || crit2 >= eps_val){
    Fi = solve(diag(dim1)-Bhat%*%Cmat, Bhat)
    Ci = solve(diag(dim1)-Bhat%*%Cmat, Ahat)
    crit1 = max(abs(Fi-Fmat))
    crit2 = max(abs(Ci-Cmat))
    Cmat = Ci 
    Fmat = Fi
    iter = iter+1
    if(iter > 500){ 
      stop("The brute-force iterative procedure did not converge after 500 iterations.")
    }
  }
  
  Gmat <- solve(diag(dim1)-Bhat%*%Cmat, Hhat)
  Pmat <- matrix(solve(diag(dim1^2) - t(Dmat)%x%Fmat, c(Gmat)), dim1, dim1)
  var_polm <- abind::abind(diag(dim1),
                           -Pmat%*%Dmat%*%solve(Pmat) - Cmat, 
                           Pmat%*%Dmat%*%solve(Pmat)%*%Cmat,
                           along = 3)
  armamod(lmfd(a = polm(var_polm), b = diag(dim1)), sigma_L = Pmat)
}

plot_irfs <- function(){
  for(prm_ix in 1:nrow(prms)){
    for(var_ix in 1:2){
      for(shock_ix in 1:2){
        irf_arr <- tt_full %>%
          filter(beta == prms[prm_ix,1],
                 nu == prms[prm_ix,2],
                 q==2) %>% 
          # mutate(value_aic = map2_dbl(.x = value_final, .y = tmpl, ~  .x + .y$n_par * 2/250)) %>%
          # group_by(q) %>%
          # slice_min(value_aic) %>%
          # ungroup() %>% 
          # filter(mb_length == prms[prm_ix]) %>%
          mutate(irf = map(.x=irf,~.x %>% unclass)) %>%
          pull(irf) %>%
          abind::abind(along=4)
        irf_qt <- apply(X = irf_arr[var_ix,shock_ix,,],
                        1, 
                        quantile, 
                        probs= c(.05, .14, .5, .86, .95))
        n_quantiles <- 1:dim(irf_qt)[1]
        n_ahead <- dim(irf_qt)[2]-1
        plot(0:n_ahead,
             irf_qt[median(n_quantiles),], 
             type = "l", 
             ylim = c(min(irf_qt),
                      max(irf_qt)),
             # main = paste0("mb_len = ", prms[prm_ix])
             main = paste("beta: ", prms[prm_ix, 1],
                          " nu: ", prms[prm_ix, 2])
        )
        lines(0:n_ahead, irf_qt[1,], col = "red")
        lines(0:n_ahead, irf_qt[max(n_quantiles),], col = "red")
        if(length(n_quantiles)==5){
          lines(0:n_ahead, irf_qt[2,], col = "blue")
          lines(0:n_ahead, irf_qt[4,], col = "blue")
        }
        abline(h=0, lty = "dashed")
      }
    }
  }
}

perm_init <- function(init_val, n_perms, ll_fun){
  
  init_list <- list(init_val)
  for(j in 1:n_perms){
    shrink_prm <- 1
    while(shrink_prm > 1e-9){
      aux_init <- rnorm(n = length(init_val),
                        mean = 0, 
                        sd = shrink_prm*init_val^2/(init_val^2+1))
      shrink_prm <- ifelse(test = abs(ll_fun(init_val + aux_init))>5*abs(ll_fun(init_val)), 
                           yes = shrink_prm/1.125,
                           no = 0) %>% as.double()
    }
    init_list[[j+1]] <- init_val + aux_init
  }
  init_list
}

# CHRIS SIMS' OPTIMIZATION ALGORITHMS

bfgsi <- function(H0, dg, dx) {
  ### dg is previous change in gradient; dx is previous change in x;
  ### 6/8/93 version that updates inverse hessian instead of hessian
  ### itself.
  ### Copyright by Christopher Sims 1996.  This material may be freely
  ### reproduced and modified.
  n <- length(dg)
  dim(dg) <- c(n,1)
  dim(dx) <- c(n,1)
  Hdg <- H0 %*% dg
  dgdx <- as.numeric(crossprod(dg,dx))
  dxHdg <- drop(dx %o% Hdg) # drops are needed to get rid of redundant unit-dimensions
  x1 <- as.numeric(1+crossprod(dg,Hdg)/dgdx)*drop(dx %o% dx)
  x2 <- dxHdg+t(dxHdg)
  ## code below causes problems when x1 and x2 have matching zeros, and I can't now (2005-12-22)
  ## figure out why I thought it was a good idea
  ##   if ( max(abs(x1-x2)/(abs(x1)+abs(x2))) <= 1e-12 ) {
  ##     cat("bfgs update failed.\n")
  ##     cat("x1, x2 = ",x1,x2,"\n")
  ##     H=H0
  ##     return(H)
  ##   }
  if (abs(dgdx) <= 1e-12)   {
    cat("bfgs update failed.\n")
    cat("|dg| =", sqrt(sum(dg^2)), "|dx| = ", sqrt(sum(dx^2)),"\n")
    cat("crossprod(dg,dx) =", dgdx,"\n")
    cat("|H*dg| =", crossprod(Hdg),"\n")
    H=H0
    return(H)
  }
  ## otherwise, 
  H <- H0 + x1/dgdx - x2/dgdx
  save(file="H.dat", H)
  return(H)
}

csminit <- function(fcn, x0, f0, g0, badg, H0,...){
  ### retcodes: 0, normal step.  5, largest step still improves too fast.
  ### 4,2 back and forth adjustment of stepsize didn't finish.  3, smallest
  ### stepsize still improves too slow.  6, no improvement found.  1, zero
  ### gradient.
  ###---------------------
  ### Fixed 7/17/93 to use inverse-hessian instead of hessian itself in bfgs
  ### update.
  ###
  ### Fixed 7/19/93 to flip eigenvalues of H to get better performance when
  ### it's not psd.
  ###
  ## ANGLE <- .0005
  ANGLE <- 1e-7
  THETA <- .01 #(0<THETA<.5) THETA near .5 makes long line searches, possibly fewer iterations.
  FCHANGE <- 1000
  MINLAMB <- 1e-9
  ### fixed 7/15/94
  ### MINDX <- .0001;
  ### MINDX <- 1e-6;
  MINDFAC <- .01
  fcount<-0
  lambda<-1
  xhat<-x0
  f<-f0
  fhat<-f0
  g <- g0
  gnorm <- sqrt(sum(g^2))
  ###
  if (!badg && (gnorm < 1.e-12)) {      # put !badg 8/4/94
    retcode <- 1
    dxnorm <- 0
    ## gradient convergence
  } else {
    ## with badg true, we don't try to match rate of improvement to directional
    ## derivative.  We're satisfied just to get some improvement in f.
    ##
    ##if(badg)
    ##   dx = -g*FCHANGE/(gnorm*gnorm);
    ##  dxnorm = norm(dx)
    ##  if dxnorm > 1e12
    ##     disp('Bad, small gradient problem.')
    ##     dx = dx*FCHANGE/dxnorm;
    ##   end
    ##else
    ## Gauss-Newton step;
    ##---------- Start of 7/19/93 mod ---------------
    ##[v d] = eig(H0);
    ##toc
    ##d=max(1e-10,abs(diag(d)));
    ##d=abs(diag(d));
    ##dx = -(v.*(ones(size(v,1),1)*d'))*(v'*g);
    dx <- -H0 %*% g
    dxnorm <- sqrt(sum(dx^2))
    if (dxnorm > 1e12){
      cat("Near-singular H problem.\n")
      dx <- dx*FCHANGE/dxnorm
    }
    dfhat <- crossprod(dx,g0)
    if(!badg){
      ## test for alignment of dx with gradient and fix if necessary
      a <- -dfhat/(gnorm*dxnorm)
      if(a<ANGLE){
        if (a < 0) {
          dx <- -g / gnorm^2
          dfhat <- -1
          cat("H unused\n")
          ## Don't bother with H.  It's not psd. Step length here appropriate for log LH,
          ## where 1.0 is a reasonable scale for changes.
        } else {
          dx <- dx - as.numeric(ANGLE*dxnorm/gnorm+dfhat/(gnorm*gnorm))*g
          dx <- dx * dxnorm / sqrt(sum(dx^2)) # This line added 2/18/2004
          dfhat <- crossprod(dx,g)
          ## dxnorm <- sqrt(sum(dx^2)) # No longer needed with 2/18/2004 change
          cat("Correct for low angle:" ,a,"\n")
        }
      }
    }
    cat(sprintf("Predicted improvement            = %18.9f\n",-dfhat/2))
    ## cat("Predicted improvement:", sprintf("%18.9f",-dfhat/2),"\n")
    ##
    ## Have OK dx, now adjust length of step (lambda) until min and
    ## max improvement rate criteria are met.
    done <- 0
    factor <- 3
    shrink <- 1
    lambdaMin <- 0
    lambdaMax <- Inf
    lambdaPeak <- 0
    fPeak <- f0
    lambdahat <- 0
    while(!done){
      ## argument of fcn retains its dim (or lack thereof), but g
      ## always emerges as n x 1
      ddx <- dx*lambda
      dim(ddx) <- dim(x0)
      dxtest <- x0+ddx
      f <- fcn(dxtest,...)
      cat(sprintf("lambda = %10.5g; f = %20.7f",lambda,f ),"\n")
      if(f < fhat){
        fhat <- f
        xhat <- dxtest
        lambdahat <- lambda
      }
      fcount <- fcount+1
      shrinkSignal <- (!badg && (f0-f < max(-THETA*dfhat*lambda,0))) || (badg && ((f0-f) < 0)) 
      growSignal <- (!badg && ( (lambda > 0)  &&  (f0-f > -(1-THETA)*dfhat*lambda) ))
      if(  shrinkSignal  &&   ( (lambda>lambdaPeak) || (lambda<0) )){
        if ((lambda>0) && ((!shrink || (lambda/factor <= lambdaPeak)))){
          shrink <- 1
          factor <- factor^.6
          while(lambda/factor <= lambdaPeak){
            factor <- factor^.6
          }
          if( abs(factor-1)<MINDFAC){
            if( abs(lambda) < 4){
              retcode <- 2
            }else{
              retcode <- 7
            }
            done <- 1
          }
        }
        if(  (lambda<lambdaMax) && (lambda>lambdaPeak) ){
          lambdaMax <- lambda
        }
        lambda <- lambda/factor
        if( abs(lambda) < MINLAMB ){
          if( (lambda > 0) & (f0 <= fhat) ){
            ## try going against gradient, which may be inaccurate
            lambda <- -lambda*factor^6
          }else{
            if( lambda < 0 ){
              retcode <- 6
            }else{
              retcode <- 3
            }
            done <- 1
          }
        }
      }else{
        if(  (growSignal && lambda>0) ||  (shrinkSignal && ((lambda <= lambdaPeak) && (lambda>0)))  ) {
          if( shrink ){
            shrink <- 0
            factor <- factor^.6
            if( abs(factor-1)<MINDFAC ) {
              if( abs(lambda)<4 ) {
                retcode <- 4
              }else{
                retcode <- 7
              }
              done <- 1
            }
          }
          if( ( f<fPeak ) && (lambda>0) ) {
            fPeak <- f
            lambdaPeak <- lambda
            if( lambdaMax<=lambdaPeak ) {
              lambdaMax <- lambdaPeak*factor*factor
            }
          }
          lambda <- lambda*factor
          if( abs(lambda) > 1e20 ) {
            retcode <- 5
            done <-1
          }
        } else {
          done <- 1
          if( factor < 1.2 ){
            retcode <- 7
          } else {
            retcode <- 0
          }
        }
      }
    }
  }
  cat(sprintf("Norm of dx %10.5g\n", dxnorm))
  return(list(fhat=fhat,xhat=xhat,fcount=fcount,retcode=retcode))
}

csminwelNew <- function(fcn, x0, H0, ..., grad=NULL, crit=1e-7, nit, Verbose=TRUE, Long=FALSE) {
  ### fcn:   the objective function to be minimized.  If it has a "gradient" attribute (like output of deriv), that is used
  ###        as analytic gradient.  If it has a "hessian" attribute, that is used as the hessian.
  ### fcn0:  Lean version of fcn, without grad or hessian attributes.  May save time to provide this. (Removed this option for now.)
  ### x0:    initial value of the parameter vector
  ### H0:    initial value for the inverse Hessian.  Must be positive definite, if used.  (Not used if attr(fcn,"hessian") exists.)
  ### grad:  If this is a numerical vector and attr(fcn,"gradient") is not present, then grad is used as an initial gradient vector.
  ###        Useful for restarting if numerical gradient calculation is slow.
  ### crit:  Convergence criterion.  Iteration will cease when it proves impossible to improve the
  ###        function value by more than crit.
  ### nit:   Maximum number of iterations.
  ### ...:   A list of optional length of additional parameters that get handed off to fcn each
  ###        time it is called.
  ###        Note that if the program ends abnormally, it is possible to retrieve the current x,
  ###        f, and H from the files g1.mat and H.mat that are written at each iteration and at each
  ###        hessian update, respectively.  (When the routine hits certain kinds of difficulty, it
  ###        writes g2.mat and g3.mat as well.  If all were written at about the same time, any of them
  ###        may be a decent starting point.  One can also start from the one with best function value.)
  dots <- list(...) # (need this to save these arguments in case of cliffs)
  nx <- length(x0)
  done <- 0
  itct <- 0
  fcount <- 0
  snit <- 100
  badg1 <- badg2 <- badg3 <- FALSE
  f0 <- fcn(x0,...)
  NumGrad <- is.null(attr(f0,"gradient"))
  NumHess <- is.null(attr(f0,"hessian"))
  badg <- FALSE
  ## browser()
  if( f0 > 1e50 ) {
    stop(paste("Bad initial parameter. f0 =", f0))
  }
  if( NumGrad ) {
    if( !is.numeric(grad) ) {
      gbadg <- numgrad(fcn,x0,...)
      g <- gbadg$g
      badg <- gbadg$badg
    } else {
      badg <- FALSE
      ## This is dangerous if you use a saved g file and it
      ## turns out to have been "bad".  We used to set badg to TRUE if
      ## grad contained any zeros.
      g <- grad
    }
  } else {
    g <- attr(f0,"gradient")
    badg <- FALSE
    gbadg <- list(g=g, badg=badg)
  }
  retcode3 <- 101
  x <- x0
  f <- f0
  if (is.null(attr(f0,"hessian"))) {
    H <- H0
  }else{
    H <- attr(f0,"hessian")
  }
  cliff <- 0
  while( !done ) {
    g1 <- NULL; g2 <- NULL; g3 <- NULL
    ##addition fj. 7/6/94 for control
    cat('-----------------\n')
    cat('-----------------\n')
    cat(sprintf('f at the beginning of new iteration, %20.10f',f),"\n")
    if (!Long && Verbose) { # set Long=TRUE if parameter vector printouts too long
      cat("x =\n")
      print(x,digits=12)
    }
    ##-------------------------
    itct <- itct+1
    itout <- csminit(fcn,x0=x, f0=f, g0=g,badg=badg, H0=H,...)
    f1 <- itout$fhat
    x1 <- itout$xhat
    fc <- itout$fcount
    retcode1 <- itout$retcode
    fcount <- fcount+fc
    if( retcode1 != 1 ) {         # Not gradient convergence
      ## if( retcode1==2 || retcode1==4) {
      ##   wall1 <- TRUE; badg1 <- TRUE
      ## } else {                          # not a back-forth wall, so check gradient
      if( NumGrad ) {
        gbadg <- numgrad(fcn, x1,...)
      } else {
        gbadg <- list(g=attr(f1,"gradient"),badg=FALSE)
      }
      g1 <- gbadg$g
      badg1 <- gbadg$badg
      wall1 <- (badg1 || retcode1==2 || retcode1 == 4) # A wall is back-forth line search close or bad gradient
      ## g1
      save(file="g1", g1, x1, f1, dots)
      ## }
      if( wall1 && dim(H)[1] > 1) {
        ## Bad gradient or back and forth on step length.  Possibly at
        ## cliff edge.  Try perturbing search direction, if problem is not unidimensional
        ##
        Hcliff <- H+diag(diag(H) * rnorm(nx))
        cat("Cliff.  Perturbing search direction. \n")
        itout <- csminit(fcn,x0=x,f0=f,g0=g, badg=badg, H0=Hcliff,...)
        f2 <- itout$fhat
        x2 <- itout$xhat
        fc <- itout$fcount
        retcode2 <- itout$retcode
        fcount <- fcount+fc   # put by Jinill
        ## if(  f2 < f ) {
        ##   if( retcode2 == 2 || retcode2==4 ){
        ##     wall2 <- 1
        ##     badg2 <- 1
        ##   } else {
        if( NumGrad ){
          gbadg <- numgrad(fcn, x2,...)
        }else{
          gbadg <- list(g=attr(f2,"gradient"),badg=FALSE)
        }
        g2 <- gbadg$g
        badg2 <- gbadg$badg
        wall2 <- (badg2 || retcode2==2 || retcode2==4)
        ## g2
        ## badg2
        save(file="g2", g2, x2, f2, dots)
        if( wall2 ){
          cat('Cliff again.  Try traversing\n')
          normdx <- sqrt(sum(x2-x1)^2)
          if( normdx < 1e-13 ) { # two trial x's too close, can't traverse
            f3 <- f; x3 <- x; badg3 <- NumGrad;retcode3 <- 101
          }else{
            ## as.numeric below for robustness against f's being 1x1 arrays
            gcliff <- (x2 - x1) * (as.numeric(f2 - f1)/(normdx^2))
            dim(gcliff) <- c(nx,1)
            itout <- csminit(fcn,x0=x,f0=f, g0=gcliff, badg=FALSE,H0=diag(nx),...)
            f3 <- itout$fhat
            x3 <- itout$xhat
            fc <- itout$fc
            retcode3 <- itout$retcode
            fcount <- fcount+fc  # put by Jinill
            ## if( retcode3==2 || retcode3==4 ) {
            ##   wall3 <- 1
            ##   badg3 <- 1
            ## } else {
            if( NumGrad ) {
              gbadg <- numgrad(fcn, x3,...) 
            }else{
              gbadg <- list(g=attr(f3,"gradient"),badg=FALSE)
            }
            g3 <- gbadg$g
            badg3 <- gbadg$badg
            wall3 <- (badg3 || retcode3==2 || retcode3==4)
            ## g3
            ## badg3
            save(file="g3", g3, x3, f3, dots)
          }
        } else { # wall1, but not wall2, so pack f3, etc with initial values
          f3 <- f
          x3 <- x
          badg3 <- NumGrad              #i.e., use the gradient if it's analytic
          g3 <- g
          retcode3 <- 101
        }
      } else {     # no walls, or one-dimensional, so no perturbations
        f3 <- f
        x3 <- x
        badg3 <- NumGrad
        g3 <- g
        retcode3 <- 101
        f2 <- f
        x2 <- x
        badg2 <- NumGrad
        g2 <- g
        retcode2 <- 101
      }
    } else { # gradient convergence --- csminit just looked at gradient and quit
      f1 <-  f2 <-  f3 <- f; g1 <- g <- FALSE; badg2 <-  badg3 <- TRUE; retcode1 <- 1; retcode2 <- 101; retcode3 <- 101
    }
    ## how to pick gh and xh:
    ## pick first fj that improves by crit and has badg==FALSE, starting with j=3
    ## may pick other than min fj to stay away from cliff
    if( f3 < f-crit & badg3==0 ) {
      ih <- 3
      fh <- f3;xh <- x3;gh <- g3;badgh <- badg3;retcodeh <- retcode3
    } else {
      if( f2 < f-crit && badg2==0 ) {
        ih <- 2
        fh <- f2;xh <- x2;gh <- g2;badgh <- badg2;retcodeh <- retcode2
      } else {
        if( f1 < f-crit && badg1==0 ) {
          ih <- 1
          fh <- f1;xh <- x1;gh <- g1;badgh <- badg1;retcodeh <- retcode1
        } else {
          ## none qualify, so pick the one that improves most.
          fh <- min(f1,f2,f3)
          ih <- which.min(c(f1,f2,f3))
          cat("ih =",ih,"\n")
          xh <- switch(ih,x1,x2,x3)
          retcodeh <- switch(ih,retcode1,retcode2,retcode3)
          nogh <- (!exists("gh",inherits=FALSE) || is.null(gh))
          if( nogh ) {
            if( NumGrad ) {
              gbadg <- numgrad(fcn,xh,...)
            } else {
              gbadg <- list(g=switch(ih,attr(f1,"gradient"),attr(f2,"gradient"),attr(f3,"gradient")),badg=FALSE)
            }
            gh <- gbadg$g
            badgh <- gbadg$badg
          }
          badgh <- NumGrad
        }
      }
    }
    ## end of picking
    ##ih
    ##fh
    ##xh
    ##gh
    ##badgh
    stuck <- (abs(fh-f) < crit)
    if ( (!badg) && (!badgh) && (!stuck) ) {
      if(NumHess){
        H <- bfgsi(H,gh-g,xh-x)
      } else {
        H <- attr(fh,"hessian")
      }
    } else {
      cat("skipped bfgsi\n")
    }
    ## if( Verbose ) {
    cat("----\n")
    cat(sprintf('Improvement on iteration %8.0f = %18.9f\n',itct,f-fh))
    ##}
    if( itct > nit ) {
      cat("iteration count termination\n")
      done <- 1
    } else {
      if( stuck ) {
        cat("improvement < crit termination\n")
        done <- 1
      }
      rc <- retcodeh
      switch( rc ,
              cat("zero gradient\n"),    #1
              cat("back and forth on step length never finished\n"), #2
              cat("smallest step still improving too slow\n"),       #3
              cat("back and forth on step length never finished\n"), #4
              cat("largest step still improving too fast\n"),        #5
              cat("smallest step still improving too slow, reversed gradient\n"), #6
              cat("warning: possible inaccuracy in H matrix\n"), #7
      )
    }
    f <- fh
    x <- xh
    g <- gh
    badg <- badgh
  }                                     # while !done
  return(list(fh=fh,xh=xh,gh=gh,H=H,itct=itct,fcount=fcount,retcodeh=retcodeh, match.call(), ...))
}

csminwel <- function(fcn, x0, H0, ..., grad=NULL, crit=1e-7, nit, Verbose=TRUE, Long=FALSE) {
  ### fcn:   the objective function to be minimized
  ### x0:    initial value of the parameter vector
  ### H0:    initial value for the inverse Hessian.  Must be positive definite.
  ### grad:  A function that calculates the gradient, or the null matrix, or a numerical
  ###        vector.  If it's not a function, numerical gradients are used.  If it's numerical, the
  ###        first iteration uses the numerical vector as the initial gradient.  This is useful when the algorithm
  ###        is restarted and the gradient calculation is slow, since the program stores the final gradient when it quits.
  ###        If it's a function, it must return a list, with elements named g and badg.  g is the gradient and badg
  ###        is logical, TRUE if the gradient routine ran into trouble so its return value should not be used.
  ### crit:  Convergence criterion.  Iteration will cease when it proves impossible to improve the
  ###        function value by more than crit.
  ### nit:   Maximum number of iterations.
  ### ...:   A list of optional length of additional parameters that get handed off to fcn each
  ###        time it is called.
  ###        Note that if the program ends abnormally, it is possible to retrieve the current x,
  ###        f, and H from the files g1.mat and H.mat that are written at each iteration and at each
  ###        hessian update, respectively.  (When the routine hits certain kinds of difficulty, it
  ###        write g2.mat and g3.mat as well.  If all were written at about the same time, any of them
  ###        may be a decent starting point.  One can also start from the one with best function value.)
  dots <- list(...)                     # (need this to save these arguments in case of cliffs)
  nx <- length(x0)
  NumGrad <- !is.function(grad)
  done <- 0
  itct <- 0
  fcount <- 0
  snit <- 100
  badg1 <- badg2 <- badg3 <- FALSE
  f0 <- fcn(x0,...)
  ## browser()
  if( f0 > 1e50 ) {
    stop("Bad initial parameter.")
  }
  if( NumGrad ) {
    if( !is.numeric(grad) ) {
      gbadg <- numgrad(fcn,x0,...)
      g <- gbadg$g
      badg <- gbadg$badg
    } else {
      badg <- any(grad==0)
      g <- grad
    }
  } else {
    gbadg <- grad(x0,...)
    badg <-gbadg$badg
    g <- gbadg$g
  }
  retcode3 <- 101
  x <- x0
  f <- f0
  H <- H0
  cliff <- 0
  while( !done ) {
    g1 <- NULL; g2 <- NULL; g3 <- NULL
    ##addition fj. 7/6/94 for control
    cat('-----------------\n')
    cat('-----------------\n')
    cat(sprintf('f at the beginning of new iteration, %20.10f',f),"\n")
    if (!Long && Verbose) {             # set Long=TRUE if parameter vector printouts too long
      cat("x =\n")
      print(x,digits=12)
    }
    ##-------------------------
    itct <- itct+1
    itout <- csminit(fcn,x,f,g,badg,H,...)
    f1 <- itout$fhat
    x1 <- itout$xhat
    fc <- itout$fcount
    retcode1 <- itout$retcode
    fcount <- fcount+fc
    if( retcode1 != 1 ) {               # Not gradient convergence
      ## if( retcode1==2 || retcode1==4) {
      ##   wall1 <- TRUE; badg1 <- TRUE
      ## } else {                          # not a back-forth wall, so check gradient
      if( NumGrad ) {
        gbadg <- numgrad(fcn, x1,...)
      } else {
        gbadg <- grad(x1,...)
      }
      g1 <- gbadg$g
      badg1 <- gbadg$badg
      wall1 <- (badg1 || retcode1==2 || retcode1 == 4) # A wall is back-forth line search close or bad gradient
      ## g1
      save(file="g1", g1, x1, f1, dots)
      ## }
      if( wall1 && dim(H)[1] > 1) {
        ## Bad gradient or back and forth on step length.  Possibly at
        ## cliff edge.  Try perturbing search direction, if problem is not unidimensional
        ##
        Hcliff <- H+diag(diag(H) * rnorm(nx))
        cat("Cliff.  Perturbing search direction. \n")
        itout <- csminit(fcn,x,f,g,badg,Hcliff,...)
        f2 <- itout$fhat
        x2 <- itout$xhat
        fc <- itout$fcount
        retcode2 <- itout$retcode
        fcount <- fcount+fc             # put by Jinill
        ## if(  f2 < f ) {
        ##   if( retcode2 == 2 || retcode2==4 ){
        ##     wall2 <- 1
        ##     badg2 <- 1
        ##   } else {
        if( NumGrad ){
          gbadg <- numgrad(fcn, x2,...)
        }else{
          gbadg <- grad(x2,...)
        }
        g2 <- gbadg$g
        badg2 <- gbadg$badg
        wall2 <- (badg2 || retcode2==2 || retcode2==4)
        ## g2
        ## badg2
        save(file="g2", g2, x2, f2, dots)
        if( wall2 ) {
          cat('Cliff again.  Try traversing\n')
          normdx <- sqrt(sum(x2-x1)^2)
          if( normdx < 1e-13 ) {        # two trial x's too close, can't traverse
            f3 <- f; x3 <- x; badg3 <- 1;retcode3 <- 101
          } else {
            ## as.numeric below for robustness against f's being 1x1 arrays
            gcliff <- (x2 - x1) * (as.numeric(f2 - f1)/(normdx^2))
            dim(gcliff) <- c(nx,1)
            itout <- csminit(fcn,x,f,gcliff,0,diag(nx),...)
            f3 <- itout$fhat
            x3 <- itout$xhat
            fc <- itout$fc
            retcode3 <- itout$retcode
            fcount <- fcount+fc         # put by Jinill
            ## if( retcode3==2 || retcode3==4 ) {
            ##   wall3 <- 1
            ##   badg3 <- 1
            ## } else {
            if( NumGrad ) {
              gbadg <- numgrad(fcn, x3,...)
            }else{
              gbadg <- grad(x3,...)
            }
            g3 <- gbadg$g
            badg3 <- gbadg$badg
            wall3 <- (badg3 || retcode3==2 || retcode3==4)
            ## g3
            ## badg3
            save(file="g3", g3, x3, f3, dots)
          } 
        } else {                        # wall1, but not wall2, so pack f3, etc with initial values
          f3 <- f
          x3 <- x
          badg3 <- 1
          retcode3 <- 101
        }
      } else {                          # no walls, or one-dimensional, so no perturbations
        f3 <- f
        x3 <- x
        badg3 <- TRUE
        retcode3 <- 101
        f2 <- f
        x2 <- f
        retcode2 <- 101
      }
    } else {                            # gradient convergence --- csminit just looked at gradient and quit
      f1 <-  f2 <-  f3 <- f; g1 <- g <- FALSE; badg2 <-  badg3 <- TRUE; retcode1 <- 1; retcode2 <- 101; retcode3 <- 101
    }
    ## how to pick gh and xh:
    ## pick first fj that improves by crit and has badg==FALSE, starting with j=3
    ## may pick other than min fj to stay away from cliff
    if( f3 < f-crit & badg3==0 ) {
      ih <- 3
      fh <- f3;xh <- x3;gh <- g3;badgh <- badg3;retcodeh <- retcode3
    } else {                              
      if( f2 < f-crit && badg2==0 ) {
        ih <- 2
        fh <- f2;xh <- x2;gh <- g2;badgh <- badg2;retcodeh <- retcode2
      } else {
        if( f1 < f-crit && badg1==0 ) {
          ih <- 1
          fh <- f1;xh <- x1;gh <- g1;badgh <- badg1;retcodeh <- retcode1
        } else {
          ## none qualify, so pick the one that improves most.
          fh <- min(f1,f2,f3)
          ih <- which.min(c(f1,f2,f3))
          cat("ih =",ih,"\n")
          xh <- switch(ih,x1,x2,x3)
          retcodeh <- switch(ih,retcode1,retcode2,retcode3)
          nogh <- (!exists("gh",inherits=FALSE) || is.null(gh))
          if( nogh ) {
            if( NumGrad ) {
              gbadg <- numgrad(fcn,xh,...)
            } else {
              gbadg <- grad( xh,...)
            }
            gh <- gbadg$gh
            badgh <- gbadg$badg
          }
          badgh <- 1
        }
      }
    }
    ## end of picking
    ##ih
    ##fh
    ##xh
    ##gh
    ##badgh
    stuck <- (abs(fh-f) < crit)
    if ( (!badg) && (!badgh) && (!stuck) ) {
      H <- bfgsi(H,gh-g,xh-x)
    }
    ## if( Verbose ) {
    cat("----\n")
    cat(sprintf('Improvement on iteration %8.0f = %18.9f\n',itct,f-fh))
    ##}
    if( itct > nit ) {
      cat("iteration count termination\n")
      done <- 1
    } else {
      if( stuck ) {
        cat("improvement < crit termination\n")
        done <- 1
      }
      rc <- retcodeh
      switch( rc ,
              cat("zero gradient\n"),    #1
              cat("back and forth on step length never finished\n"), #2
              cat("smallest step still improving too slow\n"), #3
              cat("back and forth on step length never finished\n"), #4
              cat("largest step still improving too fast\n"), #5
              cat("smallest step still improving too slow, reversed gradient\n"), #6
              cat("warning: possible inaccuracy in H matrix\n"), #7
      )
    }
    f <- fh
    x <- xh
    g <- gh
    badg <- badgh
  }                                     # while !done
  return(list(fh=fh,xh=xh,gh=gh,H=H,itct=itct,fcount=fcount,retcodeh=retcodeh,...))
}

csolve <- function(FUN, x,..., gradfun=NULL, crit=1e-7, itmax=20, verbose=TRUE, alpha=1e-3, delta=1e-6, long=FALSE) {
  ### FUN:      A function with vector argument x or (when numerical derivatives are used) with matrix argument x.
  ###           the number of rows in x matches the number of rows in the return value.  For numerical derivatives,
  ###           the number of columns in x is the number of distinct argument vectors at which FUN is evaluated.
  ###           **************************************
  ### x:        initial value for FUN's argument.
  ### gradfun:  the function called to evaluate the gradient matrix.  If this
  ###           is NULL (the default), a numerical gradient is used instead.  If it is identical to
  ###              FUN, then FUN returns a value v with attr(v,"grad") the gradient matrix.
  ### crit:     if the sum of absolute values that FUN returns is less than this,
  ###           the equation is solved.
  ### itmax:    the solver stops when this number of iterations is reached, with rc=4
  ### ...:      in this position the user can place any number of additional arguments, all
  ###           of which are passed on to FUN and gradfun (when it is non-empty) as a list of 
  ###           arguments following x.
  ### -------------------------------------------------------------------
  ###           Arguments below this usually can be left at defaults.
  ### verbose:  If set to FALSE, the amount of output during iterations is cut to zero.
  ### long:     Set to TRUE to suppress printout of x and f at each iteration. (No effect if verbose=FALSE)
  ###           Useful when x and FUN are long vectors.
  ### alpha:    Tolerance on rate of descent.  The algorithm may produce search directions nearly
  ###           orthogonal to the gradient, and hence nearly zero directional derivative.  Smaller
  ###           alpha allows closer approach to orthogonality.
  ### delta:    difference interval in (cheap, forward-difference) numerical derivative.  Ignored if gradfun non-NULL.
  ### rc:       0 means normal solution, 1 and 3 mean no solution despite extremely fine adjustments
  ###           in step length (very likely a numerical problem, or a discontinuity). 4 means itmax
  ###           termination.
  ### -------------------------------------------------------------------
  ### modified 12/23/05 to allow everything to be complex-valued.  Modifications only lightly tested.
  ###------------ analyticg --------------
  analyticg <- !is.null(gradfun) #if the grad argument is NULL, numerical derivatives are used.
  jointg <- identical(FUN, gradfun) # FUN returns both value and gradient
  ###-------------------------------------
  EPS <- .Machine$double.eps
  if (is.null(dim(x))) {
    vdn <- names(x)
    x <- matrix(x,length(x),1)
    dimnames(x)[[1]] <- vdn
  }
  nv <- dim(x)[1]
  vdn <- dimnames(x)
  tvec <- delta*diag(nv)
  done <- FALSE
  f0 <- FUN(x,...)
  af0 <- sum(abs(f0))
  af00 <- af0
  itct <- 0
  while(!done) {
    if((itct%%2)==1  && af00-af0 < crit*max(1,af0) && itct>3) {
      randomize <- TRUE
    } else {
      if(!analyticg) {
        grad <- (FUN(matrix(x,nv,nv,dimnames=vdn)+tvec,...)-matrix(f0,nv,nv))/delta
      } else {                          # use analytic gradient
        if (jointg) {
          grad <- attr(f0, "grad")
        } else {
          grad <- gradfun(x,...)
        }
      }
      ## if(is.real(grad) && is.finite(grad) && sum(abs(grad))> 4*nv^2*EPS) {
      if(is.finite(grad) && sum(abs(grad))> 4*nv^2*EPS) {
        svdg <- svd(grad)
        svdd <- svdg$d
        if(!(min(svdd)>0) || max(svdd)/min(svdd) > 1/(100*EPS)){
          svdd <- pmax(svdd,max(svdd)*1e-13)
          grad <- svdg$u%*% diag(svdd,ncol=length(svdd)) %*% t(svdg$v)
        }
        dx0 <- -solve(grad,f0)
        randomize <- FALSE
      } else {
        if(verbose){cat("gradient imaginary or infinite or zero\n")}
        randomize <- TRUE
      }
    }
    if(randomize) {
      if(verbose){cat("Random Search\n")}
      dx0 <- sqrt(sum(x*x))/matrix(rnorm(length(x)),length(x),1)
    }
    lambda <- 1
    lambdamin <- 1
    fmin <- f0
    xmin <- x
    afmin <- af0
    dxSize <- sqrt(sum(abs(dx0)^2))
    factor <- .6
    shrink <- TRUE
    subDone <- FALSE
    while(!subDone) {
      dx <- lambda*dx0
      f <- FUN(x+dx,...)
      af <- sum(abs(f))
      if(!is.nan(af) && af < afmin) {
        afmin <- af
        fmin <- f
        lambdamin <- lambda
        xmin <- x+dx
      }
      if (((lambda >0) && (is.nan(af) || (af0-af < alpha*lambda*af0))) || ((lambda<0) && (is.nan(af) || (af0-af < 0)) )) {
        if(!shrink) {
          factor <- factor^.6
          shrink <- TRUE
        }
        if(abs(lambda*(1-factor))*dxSize > .1*delta) {
          lambda <- factor*lambda
        } else {
          if( (lambda > 0) && (factor==.6) ) { #i.e., we've only been shrinking
            lambda <- -.3
          } else {
            subDone <- TRUE
            if(lambda > 0) {
              if(factor==.6) {
                rc <- 2
              } else {
                rc <- 1
              }
            } else {
              rc <- 3
            }
          }
        }
      } else {
        ## if( (lambda >0) && (af-af0 > (1-alpha)*lambda*af0) ) {
        ## I think af-af0 instead of af0-af in the line above is a long-standing error that
        ## meant this branch was never visited.  (cas 7/16/2009)
        if( (lambda >0) && (af0 - af > (1-alpha)*lambda*af0) ) {
          if(shrink) {
            factor <- factor^.6
            shrink <- FALSE
          }
          lambda <- lambda/factor
        } else {                        # good value found
          subDone <- TRUE
          rc <- 0
        }
      }
    }                                   # while ~subDone
    itct <- itct+1
    if(verbose){
      cat(paste("itct ",itct,", af ",afmin,", lambda ", lambdamin,", rc ",rc),"\n")
      if ( !long ) {
        if ( !is.complex(xmin) && !is.complex(fmin) ) {
          cat("x: \n",formatC(xmin,width=10,digits=6),"\n")
          cat("fmin:\n",formatC(fmin,width=10,digits=6),"\n")
        } else {
          cat("x: \n")
          print(as.complex(x))
          cat("fmin:\n")
          print(as.complex(fmin))
        }
      }
    }
    x <- xmin
    f0 <- fmin
    af00 <- af0
    af0 <- afmin
    if(itct >= itmax) {
      done <- TRUE
      rc <- 4
    } else {
      if(af0<crit) {
        done <- TRUE
        rc <- 0
      }
    }
  }                                     #while not done
  return(list(x=xmin,f=fmin,itcount=itct,retcode=rc))
}

derivVec <- function(ex, x, param) {
  ##    ex:    a vector of expressions
  ##     x:    a character vector of variable names w.r.t. which the expressions will be differentiated
  ## param:    a vector of parameter values that will be constant in repeated calls to fret (the returned function)
  ##-------------------
  ##  fret:    returned value; a function that when evaluated at a numerical x, 
  ##           returns the vector of expression values, but also, as attr(value, "grad"), the gradient matrix
  ##           It uses the fixed param values set in the call to derivVec.
  nq <- length(ex)
  nv <- length(x)
  param <- param                        #to put param into this namespace, so it does not need to be in every call.
  outf <- vector("expression",nq)
  for (iq in 1:nq) {
    f <- deriv(ex[iq], x)
    outf[iq] <- f
  }
  fret <- function(z) {
    fval <- vector("numeric", nq)
    gval <- matrix(0, nq, nv)
    zv <- as.vector(z)
    names(zv) <- if (is.null(names(z))) dimnames(z)[[1]] else names(z)
    for (iq in 1:nq) {
      fg <- eval(outf[iq], as.list(c(zv, param)))
      fval[iq] <- fg
      gval[iq, ] <- attr(fg, "grad")
    }
    fval <- c(fval)
    names(fval) <- names(ex)
    dimnames(gval) <- list(eq=names(ex), vbl=names(zv))
    attr(fval, "grad") <- gval
    return(fval)
  }
  return(fret)
}

numgrad <- function(fcn, x, ...) {
  ## fcn can return a vector, in which case numgrad returns a matrix.
  delta <- 1e-6
  ## delta <- 1e-8
  n <- length(x)
  ## we tolerate x's that may be n x 1, 1 x n, or R vectors (with no dim),
  ## but note that g comes out as n x k matrix regardless. 
  tvec <- delta*diag(n)
  f0 <- fcn(x,...)
  k <- length(f0)
  g <- matrix(0,n,k)
  badg <- FALSE
  for (i in 1:n){
    scale <- 1
    tvecv <- tvec[,i]
    if(is.null(dim(x))){
      tvecv <- as.vector(tvecv)
    }else{
      dim(tvecv) <- dim(x)
    }
    g0 <- (fcn(x+scale*tvecv,...) - f0)/(scale*delta)
    if (max(abs(g0))< 1e15){
      g[i, ] <- as.vector(g0)
    }else{
      cat("bad gradient ------------------------\n")
      badg <- TRUE
    }
  }
  return(list(g=g,badg=badg))
}

numHess <- function(fcn, x, ...) {
  f1 <- fcn
  n <- length(x)
  h <- matrix(0, n, n)
  f2 <- function(z, ...) { numgrad(fcn=f1, z, ...)$g}
  h <- numgrad(fcn=f2, x=x, ...)$g
  return(h)
}

volcano <- function(x) {
  rho <- sqrt(sum(x))
  f <- rho^2*exp(-rho)+.02*x[1]/(1+.001*x[1]^2)
  return(-f)
}
