# Packages ####
pkgs <- c("lubridate", "xts", "parallel", "svarmawhf", "svars", "fitdistrplus", "sgt", "tidyverse")
void = lapply(pkgs, library, character.only = TRUE)
DIM_OUT <- 4
pmap_tmpl_whf_rev = function(dim_out = DIM_OUT, p, q, kappa, k, shock_distr = "sgt", ...){
  tmpl_whf_rev(dim_out = DIM_OUT, ARorder = p, MAorder = q, kappa = kappa, k = k, shock_distr = shock_distr)
}
pmap_get_residuals_once <- function(params_deep_final, tmpl, data_list, ...){
  get_residuals_once(params_deep = params_deep_final, tmpl = tmpl, data_long = data_list)
}

# preliminaries ####
main_path <- "local_data/jobid_20221021/"
vec_files <- list.files(main_path)
tibble_list = vector("list", length(vec_files))
aux_data <- readRDS(paste0("local_data/total_data.rds"))

for (ix_file in seq_along(vec_files)){
  file_this = readRDS(paste0(main_path, vec_files[ix_file]))
  
  SCRIPT_PARAMS_this = file_this[[1]]$results_list$script_params
  
  IX_ARRAY_JOB_this = SCRIPT_PARAMS_this$IX_ARRAY_JOB
  N_MODS_this = with(SCRIPT_PARAMS_this, N_MODS_PER_CORE * N_CORES)
  
  tibble_list[[ix_file]] =  
    enframe(file_this) %>% 
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>% 
    unnest_wider(results_list) %>% 
    select(nr, params_deep_final, value_final, input_integerparams) %>% 
    mutate(n_params = map_int(params_deep_final, length)) %>% 
    unnest_wider(input_integerparams)
  if(ix_file==length(vec_files)) tt_full <- reduce(tibble_list, bind_rows)
}

tt <- tt_full %>% bind_cols(aux_data %>% select(-intersect(names(tt_full), names(aux_data)))) %>% 
  mutate(n_obs = map_dbl(data_list, ~nrow(.x))) %>% 
  mutate(punish_aic = n_params * 2/n_obs) %>% 
  mutate(punish_bic = n_params * log(n_obs)/n_obs) %>% 
  mutate(value_aic = value_final + punish_aic) %>% 
  mutate(value_bic = value_final + punish_bic) %>% 
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  mutate(res = pmap(., pmap_get_residuals_once)) %>% 
  mutate(B_mat = map2(params_deep_final, tmpl, 
                      ~fill_tmpl_whf_rev(theta = .x, 
                                         tmpl = .y)$B)) %>% 
  mutate(shocks = map2(res, B_mat, ~ solve(.y, t(.x)) %>% t()))


THRESHOLD_SW <- 0.05
tt <- tt %>% 
  mutate(sw = map(shocks, ~apply(.x, 2, FUN = function(x){shapiro.test(x)$p.value}))) %>% 
  mutate(sw_flag = map_int(sw, ~sum(.x > THRESHOLD_SW))) %>% # one component may be Gaussian
  mutate(sw_pval_sum = map_dbl(sw, sum)) %>% 
  unnest_wider(sw, names_sep = "_pval") %>% 
  arrange(desc(sw_pval_sum))

tt %>% group_by(length, prst) %>% summarise_at(vars(sw_flag), max)

THRESHOLD_JB = 0.05 
# filter good moodels by flag == 0
# H_0: Normality -> good models have small p-values
tt = tt %>% 
  mutate(jb = map(shocks, ~apply(.x, 2, FUN = function(x){tsoutliers::JarqueBera.test(x)[[1]]$p.value}))) %>% 
  mutate(jb_flag = map_int(jb, ~sum(.x > THRESHOLD_JB))) %>% # one component may be Gaussian
  mutate(jb_pval_sum = map_dbl(jb, sum)) %>% 
  unnest_wider(jb, names_sep = "_pval") %>% 
  arrange(desc(jb_pval_sum))

tt %>% group_by(length, prst) %>% summarise_at(vars(jb_flag), mean)

tt = tt %>% mutate(normality_flag = sw_flag + jb_flag)

tt %>% group_by(length, prst) %>% summarise(across(normality_flag, ~ n_distinct(.x)))

THRESHOLD_LB = 0.05
  # filter good moodels by flag == 0
# H_0: No autocorrelation of (transformation of) residuals
# -> good models have high p-values
tt = tt %>% 
  mutate(lb = map(shocks, ~ apply(.x, 2, FUN = function(x){Box.test(x, lag = 24, type = "Ljung-Box")$p.value }))) %>% 
  mutate(lb_flag = map_lgl(lb, ~any(.x < THRESHOLD_LB))) %>% 
  mutate(lb_pval_sum = map_dbl(lb, sum)) %>%
  unnest_wider(lb, names_sep = "_pval") %>% 
  
  mutate(lb_abs = map(shocks, ~ apply(.x, 2, FUN = function(x){Box.test(abs(x), lag = 24, type = "Ljung-Box")$p.value }))) %>% 
  mutate(lb_abs_flag = map_lgl(lb_abs, ~any(.x < THRESHOLD_LB))) %>% 
  mutate(lb_abs_pval_sum = map_dbl(lb_abs, sum)) %>%
  unnest_wider(lb_abs, names_sep = "_pval") %>% 
  
  mutate(lb_sq = map(shocks, ~ apply(.x, 2, FUN = function(x){Box.test(x^2, lag = 24, type = "Ljung-Box")$p.value }))) %>% 
  mutate(lb_sq_flag = map_lgl(lb_sq, ~any(.x < THRESHOLD_LB))) %>% 
  mutate(lb_sq_pval_sum = map_dbl(lb_sq, sum)) %>%
  unnest_wider(lb_sq, names_sep = "_pval") %>% 
  
  mutate(lb_all_pval_sum = lb_pval_sum + lb_abs_pval_sum + lb_sq_pval_sum) %>% 
  arrange(lb_all_pval_sum)

tt = tt %>% mutate(indep_flag = lb_flag + lb_abs_flag + lb_sq_flag)

tt %>% group_by(length, prst) %>% summarise_at(vars(lb_flag), mean)
tt %>% group_by(length, prst) %>% summarise_at(vars(lb_abs_flag), mean)
tt %>% group_by(length, prst) %>% summarise_at(vars(lb_sq_flag), mean)
tt %>% group_by(length, prst) %>% summarise_at(vars(indep_flag), min)
tt %>% select(contains("pval4")) %>% apply(2, mean)


# helper functions ####
permute_chgsign = function(irf_array, 
                           perm = rep(1, dim(irf_array)[1]), 
                           sign = rep(1, dim(irf_array)[1])){
  
  dim_out = dim(irf_array)[1]
  
  perm_mat = diag(dim_out)
  perm_mat = perm_mat[,perm]
  
  sign = diag(sign)
  
  ll = map(1:dim(irf_array)[3], ~ irf_array[,,.x] %*% perm_mat %*% sign)
  
  irf_array = array(0, dim = c(dim_out, dim_out, dim(irf_array)[3]))
  for (ix in 1:dim(irf_array)[3]){
    irf_array[,,ix] = ll[[ix]] 
  }
  
  return(irf_array)
}
rotmat <- function(x, n){
  
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

ff <- function(x, zero_ix, input_mat){
  
  row_ix <- if(is.null(dim(zero_ix))) zero_ix[1] else zero_ix[,1]
  col_ix <- if(is.null(dim(zero_ix))) zero_ix[2] else zero_ix[,2]
  sum(sqrt(diag(input_mat[row_ix,] %*% rotmat(x, ncol(input_mat))[, col_ix])^2))
}

# for imposing the zero restriction, take the resulting VARMA model,
# generate the corresponding ll-function, conditional on the model
# parameters optimize the likelihood function wrt the rotation to
# obtain rotation which makes the desired element equal to zero
# and maximizes the likelihood function.

for (ix_file in seq_along(vec_files)){
  file_this = readRDS(paste0(main_path, vec_files[ix_file]))
  
  SCRIPT_PARAMS_this = file_this[[1]]$results_list$script_params
  N_NOISE_PARAMS_SGT = DIM_OUT^2 + DIM_OUT*3
  
  IX_ARRAY_JOB_this = SCRIPT_PARAMS_this$IX_ARRAY_JOB
  N_MODS_this = with(SCRIPT_PARAMS_this, N_MODS_PER_CORE * N_CORES)
  
  tibble_list[[ix_file]] =  
    enframe(file_this) %>% 
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>% 
    unnest_wider(results_list) %>% 
    select(nr, params_deep_final, value_final, input_integerparams) %>% 
    mutate(n_params = map_int(params_deep_final, length)) %>% 
    mutate(n_params_sys = n_params - N_NOISE_PARAMS_SGT) %>% # SGT has (too) many noise parameters, so I want to check what happens when one uses only the number of system parameters as punishment
    unnest_wider(input_integerparams) %>% 
    mutate(punish_aic = n_params * 2/SCRIPT_PARAMS_this$N_OBS) %>% 
    mutate(punish_bic = n_params * log(SCRIPT_PARAMS_this$N_OBS)/SCRIPT_PARAMS_this$N_OBS) %>% 
    mutate(value_aic = value_final + punish_aic) %>% 
    mutate(value_bic = value_final + punish_bic) %>% 
    select(-value_final, -starts_with("punish"))
}

tt_full = reduce(tibble_list, bind_rows) %>% 
  bind_cols(tt %>% select(setdiff(names(tt), names(tt_full)))) %>% 
  group_by(prm_ix, mc_ix, n_unst) %>% 
  slice_min(order_by = value_bic) %>% 
  ungroup %>%
  mutate(irf = map2(params_deep_final, tmpl, ~irf_whf(.x, .y, n_lags = n_ahead)))

irf_svarma <- array(0, c(2, 2, n_ahead+1, max(tt_full$mc_ix), max(tt_full$prm_ix) + max(tt_full$n_unst)))

# sign and permute irfs and impose the zero restriction
for(para_ix in 1:max(tt_full$prm_ix)){
  for(root_ix in 1:max(tt_full$n_unst)){
    for(mc_i in 1:dim(irf_svarma)[4]){
      irf0 <- tt_full %>% 
        filter(prm_ix==para_ix) %>% 
        select(irf) %>%
        slice(mc_i) %>% 
        deframe %>% .[[1]] %>% 
        unclass
      perm_ix <- order(apply(abs(sapply(1:dim(irf0)[3], function(x) diag(irf0[,,x]))/diag(irf0[,,1])), 1, max))
      sign_ix <- sapply(1:dim(irf0)[2], function(x) ifelse(abs(min(irf0[,x,]))>abs(max(irf0[,x,])), -1, 1))
      irf0 <- permute_chgsign(irf0, perm = perm_ix, sign = sign_ix[perm_ix])
      irf_svarma[,,,mc_i,(para_ix-1)+root_ix] <- (pseries(irf0, n_ahead)%r%rotmat(atan(-irf0[1,2,1]/irf0[1,1,1]),2)) %>% unclass
    }
  }
}
saveRDS(tt_full, file = "./local_data/data_list.rds")
saveRDS(irf_svarma, file = "./local_data/irf_svarma.rds")