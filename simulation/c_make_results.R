# ----------------------------------------------------------------------- #
# This script creates model output: diagnostics, model selection, IRFs... #
# ----------------------------------------------------------------------- #

pkgs = c("tidyverse", "svarmawhf")
void = lapply(pkgs, library, character.only = TRUE)
select <- dplyr::select
params <- list(PATH = "local_data/jobid_",
               JOBID = "20221123")
n_ahead <- 8
# functions for the analysi
get_llf <- function(p, q, kappa, k, dtype)
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
  ll_whf_factory(data_wide = t(data_i), tmpl = tmpl_i, shock_distr = "sgt")
  
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
                     n.burnin = 2*nobs,
                     n.obs = nobs)
  return(list(y = data_out, 
              mod = re_mod, 
              prms = c("beta" = beta, "rho" = rho, "nobs" = nobs, "nu" = nu))
  )
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

get_rest_irf <- function(tbl_slice, rest_ix)
{
  nobs <- tbl_slice$nobs
  sd_mat <- readRDS("local_data/svarma_data_list.rds") %>% 
    filter(mp_type%in%tbl_slice$mp_type) %>%
    pull(std_dev) %>% .[[1]] %>%  diag
  
  irf_out <- sd_mat %r% tbl_slice$irf[[1]]
  
  rmat <- optim_zr(irf_out[,,1], rest_ix)
  llf0 <- do.call(get_llf, tbl_slice %>% select(p, q, kappa, k, mp_type) %>% rename(dtype = mp_type))
  params0 <- tbl_slice$params_deep_final[[1]]
  
  params0[tbl_slice$params_deep_final[[1]] %in% tbl_slice$B_mat[[1]]] <- c(tbl0$B_mat[[1]]%*%rmat)
  pval <- 1-pchisq(2*nobs*(llf0(params0)-tbl0$value_final), length(rest_ix))
  
  irf_out <- irf_whf(params0, tbl_slice$tmpl[[1]], 48)
  
  return(list(irf = irf_out, pval = pval, rmat = rmat))
}

pmap_tmpl_whf_rev = function(dim_out = DIM_OUT, p, q, kappa, k, shock_distr = "sgt", ...)
{
  tmpl_whf_rev(dim_out = DIM_OUT, ARorder = p, MAorder = q, kappa = kappa, k = k, shock_distr = shock_distr)
}

pmap_get_residuals_once = function(params_deep_final, tmpl, data_list, ...)
{
  get_residuals_once(params_deep = params_deep_final, tmpl = tmpl, data_long = data_list)
}

get_fevd <- function (irf_arr) 
{
  nvar <- dim(irf_arr)[1]
  n.ahead <- dim(irf_arr)[3]
  fe <- list()
  for (i in 1:nvar) {
    fe[[i]] <- as.data.frame(t(irf_arr[i, , ]))
  }
  fe2 <- fe
  for (i in 1:length(fe)) {
    for (j in 1:n.ahead) {
      fe2[[i]][j, ] <- (colSums(fe[[i]][j:1, ]^2)/sum(fe[[i]][j:1, ]^2)) * 100
    }
  }
  return(fe2)
}

choose_perm_sign <- function(target_mat, cand_mat)
{
  nvar <- ncol(target_mat)
  sign_ix <- replicate(nvar, list(c(-1,1))) %>% expand.grid
  perm_ix <- replicate(nvar, list(1:nvar)) %>% 
    expand.grid %>% 
    filter(apply(., 1, n_distinct)==nvar)
  cr0 <- 1e25
  for(j in 1:nrow(sign_ix)){
    for(jj in 1:nrow(perm_ix)){
      rt_mat <- diag(sign_ix[j, ]) %*% diag(nvar)[, unlist(perm_ix[jj, ])]
      x1 <- cand_mat %*% rt_mat
      cr1 <- norm(target_mat - x1, "F")
      if(cr0 > cr1){
        x_opt <- x1
        rt_opt <- rt_mat
        cr0 <- cr1
      }
    }
  }
  return(list(x_opt, rt_opt))
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

vec_files = list.files(paste0(params$PATH, params$JOBID))
vec_files = vec_files[grepl("arrayjob", vec_files)]
SCRIPT_PARAMS = readRDS(paste0(params$PATH, params$JOBID, "/", vec_files[1]))[[1]]$results_list$script_params
DIM_OUT = SCRIPT_PARAMS$DIM_OUT

tt_full <- tibble()
for (ix_file in seq_along(vec_files))
{
  file_this = readRDS(paste0(params$PATH, params$JOBID, "/",
                             vec_files[ix_file]))
  
  SCRIPT_PARAMS_this = file_this[[1]]$results_list$script_params
  
  IX_ARRAY_JOB_this = SCRIPT_PARAMS_this$IX_ARRAY_JOB
  N_MODS_this = with(SCRIPT_PARAMS_this, N_MODS_PER_CORE * N_CORES)
  
  tt_full =  
    enframe(file_this) %>% 
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>% 
    unnest_wider(results_list) %>% 
    select(nr, params_deep_final, value_final, input_integerparams) %>% 
    mutate(n_params = map_int(params_deep_final, length)) %>% 
    unnest_wider(input_integerparams) %>% 
    mutate(readRDS("local_data/total_data_sim.rds") %>% slice(nr)) %>% 
    mutate(punish_aic = n_params * 2/nobs) %>% 
    mutate(punish_bic = n_params * log(nobs)/nobs) %>% 
    mutate(value_aic = value_final + punish_aic) %>% 
    mutate(value_bic = value_final + punish_bic) %>% 
    mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
    mutate(res = pmap(., pmap_get_residuals_once)) %>% 
    mutate(B_mat = map2(params_deep_final, tmpl, 
                        ~fill_tmpl_whf_rev(theta = .x, 
                                           tmpl = .y)$B)) %>% 
    mutate(shocks = map2(res, B_mat, ~ solve(.y, t(.x)) %>% t())) %>%
    select(nr, p, q, kappa, k, n_st, n_unst, value_final, value_aic, value_bic, nobs, beta, rho, nu, mc_ix, B_mat, shocks, params_deep_final) %>% 
    bind_rows(tt_full)
  #mutate(cov_shocks = map(shocks, function(x){y = abs(cov(x) - diag(DIM_OUT)); names(y) = paste0("cov_el_", letters[1:(DIM_OUT^2)]); y})) %>% 
  #unnest_wider(cov_shocks) %>% 
  #mutate(cov_el_sum = rowSums(across(contains("cov_el")))) # %>% select(-tmpl, -starts_with("punish"), -res, -B_mat)
}

tt_opt <- tt_full %>% 
  group_by(mc_ix, n_unst, beta, nu) %>% 
  slice_min(value_aic) %>% 
  ungroup() %>% 
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~irf_whf(.x, .y, n_ahead)))

mc_n <- unique(tt_opt$mc_ix)
prms <- expand.grid(beta = unique(tt_opt$beta), nu = unique(tt_opt$nu))
# dim_out x dim_out x n_ahead x {(beta_1,nu_1), ...,(beta_2,nu_2)} x {mc_1, mc_2, ..., mc_n} x {VAR, VARMA}
irf_arr <- array(NA, c(DIM_OUT, DIM_OUT, n_ahead+1, nrow(prms), length(mc_n), 2))

for(prm_ix in 1:nrow(prms)){
  for(mc_i in mc_n){
    for(mod_ix in 0:1){
      # Extract IRF
      irf_i <- tt_opt %>% 
        filter(beta == prms[prm_ix,1], 
               nu == prms[prm_ix,2], 
               mc_ix == mc_i, 
               n_unst == mod_ix) %>% 
        pull(irf) %>% .[[1]]
      # Order news shock last by identifying it from FEVD
      fevd_i <- irf_i %>% unclass %>% get_fevd
      # Contribution per horizon normalized w.r.t. initial period
      max_i <- lapply(fevd_i, function(x)  t(t(x)/unlist(x[1,])) %>% apply(2, max))
      # Which horizon accounts for the most contribution w.r.t. initial period?
      whichMax_i <- sapply(max_i, which.max)
      # Choose the column of FEVD which shows a) the largest jump w.r.t. init. 
      # period across variables or b) the largest delayed jump in every variable,
      # which should be the case with pure news shock
      if(length(unique(whichMax_i))>1){
        news_ix <- max_i %>% sapply(max) %>% which.max()
      } else{
        news_ix <- unique(whichMax_i)
      }
      irf_i <- irf_i%r%diag(DIM_OUT)[, c((1:DIM_OUT)[-news_ix], news_ix)]
      # Impose zero restriction 
      rt_zero <- optim_zr(input_mat = irf_i[,,1],
                          zr_ix = c(1,2),
                          opt_it = FALSE)
      irf_i <- irf_i%r%rt_zero
      # Impose positive signs of the shocks: 1) conventional shock; 2) news shock
      if(all(irf_i[,1,1]<0)) irf_i <- irf_i%r%diag(c(-1, 1))
      max_ix <- apply(abs(irf_i[,DIM_OUT,]), 1, which.max)
      if(all(apply(irf_i[,DIM_OUT,max_ix], 1, max) < 0)) irf_i <- irf_i%r%diag(c(1,-1))
      irf_arr[,,,prm_ix,mc_i,mod_ix+1] <- unclass(irf_i)
    }
  }
}

margs_plus_font <- theme(plot.title = element_text(size = 8, hjust = 0.5),
                         plot.margin=grid::unit(c(1, 1, 0, 0), "mm"),
                         axis.text.x = element_text(size = 6),
                         axis.text.y = element_text(size = 6))

plot1 <- map2(.x = sim_comparison(irf_arr, c(2, 2), c(.14, .86), prms),
              .y = list(c(-0.2, 1.5), c(-0.25, 2.2),
                        c(-0.2, 1.5), c(-0.25, 2.2)),
                 ~ .x + margs_plus_font + ylim(.y)
                )

plot2 <- map2(.x = sim_comparison(irf_arr, c(1, 2), c(.14,.86), prms),
              .y = list(c(-0.3,1.2), c(-0.25,1.2),
                        c(-0.3,1.2), c(-0.25,1.2)),
                 ~ .x + margs_plus_font + ylim(.y)
              )

ggsave(filename = "./paper_output/sim_plt.pdf",
       gridExtra::marrangeGrob(plot1, nrow = 2, ncol = 2, top = NULL),
       width = 7.5,
       height = 5)
ggsave(filename = "./paper_output/sim_plt_a.pdf",
       plot1[[1]],
       width = 7.5,
       height = 5)
ggsave(filename = "./paper_output/sim_plt_b.pdf",
       plot1[[3]],
       width = 7.5,
       height = 5)
ggsave(filename = "./paper_output/sim_plt_d.pdf",
       plot1[[4]],
       width = 7.5,
       height = 5)

ggsave(filename = "./paper_output/appendix1.pdf",
       gridExtra::marrangeGrob(plot2, nrow = 2, ncol = 2, top = NULL),
       width = 7.5,
       height = 5)
