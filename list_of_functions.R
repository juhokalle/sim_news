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

sim_news <- function(beta, rho, no_sim = FALSE, nobs = 250, nu = 3)
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
                     rand.gen =  function(x) stats::rt(x, nu),
                     n.burnin = 1000,
                     n.obs = nobs)
  return(list(y = data_out, 
              mod = re_mod, 
              prms = c("beta" = beta, "rho" = rho, "nobs" = nobs, "nu" = nu))
         )
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

id_policy_shox <- function(irf_arr, policy_var)
{
  
  dim_out <- dim(irf_arr)[1]
  n_ahead <- dim(irf_arr)[3]
  fevd_obj <- get_fevd(irf_arr, int_var = policy_var)[[1]]
  combs <- combn(dim_out, 2)
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