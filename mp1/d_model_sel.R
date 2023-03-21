# ----------------------------------------------------------------------- #
# This script creates model output: diagnostics, model selection, IRFs... #
# ----------------------------------------------------------------------- #

source("list_of_functions.R")
pkgs = c("tidyverse", "svarmawhf")
void = lapply(pkgs, library, character.only = TRUE)
select <- dplyr::select
params <- list(PATH = "local_data/jobid_",
               JOBID = "20230223")

# sftp::sftp_connect(server = "turso.cs.helsinki.fi",
#                    folder = "/proj/juhokois/sim_news/local_data/",
#                    username = "juhokois",
#                    password = "***") -> scnx
# sftp::sftp_download(file = "jobid_20230201.zip",
#                     tofolder = "/local_data/",
#                     sftp_connection = scnx)
vec_files = list.files(paste0(params$PATH, params$JOBID))
vec_files = vec_files[grepl("arrayjob", vec_files)]
SCRIPT_PARAMS = readRDS(paste0(params$PATH, params$JOBID, "/", vec_files[1]))[[1]]$results_list$script_params
DIM_OUT = SCRIPT_PARAMS$DIM_OUT
  
tibble_list <- vector("list", length(vec_files))
TOTAL_DATA <- readRDS(paste0(params$PATH, params$JOBID, "/total_data.rds"))

for (ix_file in seq_along(vec_files)){
  
  file_this = readRDS(paste0(params$PATH, params$JOBID, "/",
                             vec_files[ix_file]))
  
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
    unnest_wider(input_integerparams) %>% 
    mutate(TOTAL_DATA %>% slice(nr)) %>% 
    mutate(nobs = map_dbl(data_list, ~nrow(.x))) %>% 
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
    select(nr, p, q, kappa, k, n_st, n_unst,
           value_final, value_aic, value_bic, nobs,
           long_sample, short_sample, fgr_sample,
           log_diff, pi, log_level,
           # mp_type, mpr_lvl, log_level,
           shock_distr, 
           B_mat, shocks, res, params_deep_final, tmpl)
    #mutate(cov_shocks = map(shocks, function(x){y = abs(cov(x) - diag(DIM_OUT)); names(y) = paste0("cov_el_", letters[1:(DIM_OUT^2)]); y})) %>% 
    #unnest_wider(cov_shocks) %>% 
    #mutate(cov_el_sum = rowSums(across(contains("cov_el")))) # %>% select(-tmpl, -starts_with("punish"), -res, -B_mat)
}

tt_full = reduce(tibble_list, bind_rows)
tt = tt_full %>% 
  mutate(rk_aic = rank(value_aic),
         rk_bic = rank(value_bic),
         rk_mle = rank(value_final)) %>% 
  arrange(value_aic) %>% 
  mutate(p_plus_q = p+q)

THRESHOLD_SW = 0.05 

# filter good models by flag == 0
# H_0: Normality -> good models have small p-values
tt = tt %>% 
  mutate(sw = map(shocks, ~apply(.x, 2, FUN = function(x){shapiro.test(x)$p.value}))) %>% 
  mutate(sw_flag = map_int(sw, ~sum(.x > THRESHOLD_SW))) %>% # one component may be Gaussian
  mutate(sw_pval_sum = map_dbl(sw, sum)) %>% 
  unnest_wider(sw, names_sep = "_pval") %>% 
  arrange(desc(sw_pval_sum))

tt %>% 
  pull(sw_flag) %>% table()

THRESHOLD_JB = 0.05 

# filter good moodels by flag == 0
# H_0: Normality -> good models have small p-values
tt = tt %>% 
  mutate(jb = map(shocks, ~apply(.x, 2, FUN = function(x){tsoutliers::JarqueBera.test(x)[[1]]$p.value}))) %>% 
  mutate(jb_flag = map_int(jb, ~sum(.x > THRESHOLD_JB))) %>% # one component may be Gaussian
  mutate(jb_pval_sum = map_dbl(jb, sum)) %>% 
  unnest_wider(jb, names_sep = "_pval") %>% 
  arrange(desc(jb_pval_sum))

tt %>% 
  pull(jb_flag) %>% table()

tt = tt %>% mutate(normality_flag = sw_flag + jb_flag)
tt %>% pull(normality_flag) %>% table()

THRESHOLD_LB = 0.05

# filter good moodels by flag == 0
# H_0: No autocorrelation of (transformation of) residuals
# -> good models have high p-values
tt = tt %>% 
  mutate(lb = map2(.x = shocks, .y = p_plus_q, ~ apply(.x, 2, FUN = function(x){ Box.test(x, lag = 24, type = "Ljung-Box")$p.value }))) %>% 
  mutate(lb_flag = map_lgl(lb, ~ any(.x < THRESHOLD_LB))) %>% 
  mutate(lb_pval_sum = map_dbl(lb, sum)) %>%
  unnest_wider(lb, names_sep = "_pval") %>% 
  
  mutate(lb_abs = map2(.x = shocks, .y = p_plus_q, ~ apply(.x, 2, FUN = function(x){ Box.test(abs(x), lag = 24, type = "Ljung-Box")$p.value }))) %>% 
  mutate(lb_abs_flag = map_lgl(lb_abs, ~any(.x < THRESHOLD_LB))) %>% 
  mutate(lb_abs_pval_sum = map_dbl(lb_abs, sum)) %>%
  unnest_wider(lb_abs, names_sep = "_pval") %>% 
  
  mutate(lb_sq = map2(.x = shocks, .y = p_plus_q, ~ apply(.x, 2, FUN = function(x){ Box.test(x^2, lag = 24, type = "Ljung-Box")$p.value }))) %>% 
  mutate(lb_sq_flag = map_lgl(lb_sq, ~any(.x < THRESHOLD_LB))) %>% 
  mutate(lb_sq_pval_sum = map_dbl(lb_sq, sum)) %>%
  unnest_wider(lb_sq, names_sep = "_pval") %>% 
  
  mutate(lb_all_pval_sum = lb_pval_sum + lb_abs_pval_sum + lb_sq_pval_sum) %>% 
  arrange(lb_all_pval_sum)

tt = tt %>% mutate(indep_flag = lb_flag + lb_abs_flag + lb_sq_flag)
tt %>% pull(indep_flag) %>% table()

tt <- tt %>% mutate(norm_indep_flag = indep_flag+normality_flag)

tt %>% pull(norm_indep_flag) %>% table

tt %>%
  #mutate(n_params = map_int(params_deep_final, length)) %>% 
  filter(norm_indep_flag==0) %>%
  group_by(log_level, pi, log_diff,
           FFR, GS1,
           short_sample, long_sample, fgr_sample) %>%
  summarise(n=n()) %>%
  arrange(n)
  pivot_wider(names_from = shock_distr, values_from = n)

irf_arr <- tt %>%
  # Filter models according to some criteria
  filter(norm_indep_flag==0,
         value_final != 1e25,
         log_level,
         FFR,
         fgr_sample) %>% 
  arrange(value_bic) %>% 
  # Merge data
  mutate(TOTAL_DATA %>% slice(nr) %>% dplyr::select(-sd)) %>%
  # Save position of structural impact matrix in prm vector
  mutate(aux_ix = map2(.x = params_deep_final, .y = B_mat, ~ .x%in%c(.y))) %>%
  # Calculate unique rotation for the B matrix and replace old
  mutate(B_mat = map(.x = B_mat, ~ choose_perm_sign(cand_mat = .x, type = "dg_abs")[[1]])) %>%
  # Replace old vector of B matrix values with the rotated ones
  mutate(params_deep_final = pmap(list(x = params_deep_final, y = aux_ix, z = B_mat), function(x,y,z) replace(x, y, c(z)))) %>% 
  # Calculate unique irf
  mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~ irf_whf(.x, .y, n_lags = 48))) %>%
  # Shocks, also rotated now appropriately
  mutate(shocks = map2(res, B_mat, ~ solve(.y, t(.x)) %>% t())) %>%
  # Shock covariance
  mutate(Sigma = map(.x = shocks, ~ cov(.x))) %>% 
  # Normalize to one std dev shocks
  mutate(irf = map2(.x = irf, .y = Sigma, ~ .x%r%diag(diag(.y)^-.5))) %>% 
  # Finally multiply the irf by the resp. std devs to get units in the original scale
  mutate(irf = map2(.x = std_dev, .y = irf, ~  diag(.x)%r%.y)) %>%
  # save FEVD for obtaining boot-CI's, this can be time consuming
  mutate(ffr_fevd = map(.x = irf, ~ get_fevd(.x, int_var = 3, by_arg = 12))) %>% 
  arrange(value_bic)

irf_arr$irf %>% 
  map(~ .x %>% unclass) %>% 
  abind::abind(along=4) %>% 
  apply(c(1,2,3), mean) %>%
  #rmfd4dfm:::cum_irf(c(5,5,1)) %>% 
  pseries(lag.max=48) %>% 
  plot()


plot(pseries(irf_md, lag.max=48))


for(j in 1:6){
  irf_arr %>% 
    slice(j) %>% 
    pull(irf) %>% .[[1]] %>%
    unclass %>% 
    rmfd4dfm:::cum_irf(trans_ix = c(5,5,1)) %>% 
    pseries(lag.max = 48) %>% 
    plot
}

saveRDS(irf_arr %>% slice(1), "local_data/target_model.rds")

irf_tot <- irf_arr %>% 
  mutate(irf = map(.x = irf, ~ unclass(.x))) %>% 
  pull(irf) %>% 
  abind::abind(along=4)

irf_md <- apply(irf_tot, c(1,2,3), mean)

n_ahead <- 48
tbl0 <- tt %>% 
  filter(nr %in% 10817) %>% 
  mutate(TOTAL_DATA %>% slice(nr) %>% dplyr::select(-sd)) %>% 
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~ irf_whf(.x, .y, lag.max = n_ahead))) %>% 
  mutate(irf = map2(.x = irf, .y = std_dev, ~ diag(.y)%r%.x)) #%>% 
  #mutate(irf = map2(.x = irf, .y = shocks, ~ .x%r%diag(diag(var(.y)))))

irf_out[,c(3,4),,drop=FALSE] %>%
  plot_irf(var_name = c("IP", "CPI", "FFR", "MPR"),
           shock_name = c("Monetary policy", "Forward guidance")) -> p1

ggsave(filename = "paper_output/IRF1.pdf", plot = p1)
