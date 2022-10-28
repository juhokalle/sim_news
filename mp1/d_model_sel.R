pkgs = c("tidyverse", "svarmawhf")
select <- dplyr::select
void = lapply(pkgs, library, character.only = TRUE)
params <- list(PATH = "local_data/jobid_",
               JOBID = "20221024")

# sftp::sftp_connect(server = "turso.cs.helsinki.fi",
#                    folder = "/proj/juhokois/sim_news/local_data/",
#                    username = "juhokois",
#                    password = "***") -> scnx
# file_ix <- 1
# file_dl <- NULL
# while(!inherits(file_dl, 'try-error')){
#   
#   sftp::sftp_download(paste0("arrayjob_", if(file_ix<10) "0", file_ix, ".rds"),
#                       tofolder = "/local_data/jobid_20221024/",
#                       sftp_connection = scnx) %>% try() -> file_dl
#   file_ix <- file_ix + 1 
# }
# sftp::sftp_download(file = "total_data.rds",
#                     tofolder = "/local_data/",
#                     sftp_connection = scnx)
vec_files = list.files(paste0(params$PATH, params$JOBID))
vec_files = vec_files[grepl("arrayjob", vec_files)]
SCRIPT_PARAMS = readRDS(paste0(params$PATH, params$JOBID, "/", vec_files[1]))[[1]]$results_list$script_params
SCRIPT_PARAMS$FILE_NAME_INPUT <-  "local_data/svarma_data.rds"
DIM_OUT = SCRIPT_PARAMS$DIM_OUT
total_data <- readRDS("local_data/total_data.rds")
  
pmap_tmpl_whf_rev = function(dim_out = DIM_OUT, p, q, kappa, k, shock_distr = "sgt", ...){
  tmpl_whf_rev(dim_out = DIM_OUT, ARorder = p, MAorder = q, kappa = kappa, k = k, shock_distr = shock_distr)
}

pmap_get_residuals_once = function(params_deep_final, tmpl, data_list, ...){
  get_residuals_once(params_deep = params_deep_final, tmpl = tmpl, data_long = data_list)
}

tibble_list = vector("list", length(vec_files))

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
    mutate(total_data %>% slice(nr)) %>% 
    mutate(punish_aic = map2_dbl(.x = data_list, .y = n_params, ~ .y * 2/nrow(.x))) %>% 
    mutate(punish_bic = map2_dbl(.x = data_list, .y = n_params, ~ .y * log(nrow(.x))/nrow(.x))) %>% 
    mutate(value_aic = value_final + punish_aic) %>% 
    mutate(value_bic = value_final + punish_bic) %>% 
    mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
    mutate(res = pmap(., pmap_get_residuals_once)) %>% 
    mutate(B_mat = map2(params_deep_final, tmpl, 
                        ~fill_tmpl_whf_rev(theta = .x, 
                                           tmpl = .y)$B)) %>% 
    mutate(shocks = map2(res, B_mat, ~ solve(.y, t(.x)) %>% t())) #%>% 
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
  mutate(p_plus_q = p+q,
         n_unstable = kappa * DIM_OUT + k) %>% 
  select(-params_deep_final, -kappa, -k, -n_params, -contains("punish"), -tmpl, -res, -B_mat) %>% 
  select(nr, p_plus_q, p, q, n_unstable, rk_aic, rk_bic, rk_mle, everything())

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
  mutate(lb = map2(.x = shocks, .y = p_plus_q, ~ apply(.x, 2, FUN = function(x){ Box.test(x, lag = 24)$p.value }))) %>% 
  mutate(lb_flag = map_lgl(lb, ~ any(.x < THRESHOLD_LB))) %>% 
  mutate(lb_pval_sum = map_dbl(lb, sum)) %>%
  unnest_wider(lb, names_sep = "_pval") %>% 
  
  mutate(lb_abs = map2(.x = shocks, .y = p_plus_q, ~ apply(.x, 2, FUN = function(x){ Box.test(abs(x), lag = 24)$p.value }))) %>% 
  mutate(lb_abs_flag = map_lgl(lb_abs, ~any(.x < THRESHOLD_LB))) %>% 
  mutate(lb_abs_pval_sum = map_dbl(lb_abs, sum)) %>%
  unnest_wider(lb_abs, names_sep = "_pval") %>% 
  
  mutate(lb_sq = map2(.x = shocks, .y = p_plus_q, ~ apply(.x, 2, FUN = function(x){ Box.test(x^2, lag = 24)$p.value }))) %>% 
  mutate(lb_sq_flag = map_lgl(lb_sq, ~any(.x < THRESHOLD_LB))) %>% 
  mutate(lb_sq_pval_sum = map_dbl(lb_sq, sum)) %>%
  unnest_wider(lb_sq, names_sep = "_pval") %>% 
  
  mutate(lb_all_pval_sum = lb_pval_sum + lb_abs_pval_sum + lb_sq_pval_sum) %>% 
  arrange(lb_all_pval_sum)

tt %>% 
  mutate(indep_flag = lb_flag + lb_abs_flag + lb_sq_flag) %>% 
  filter(length=="short") %>%
  group_by(p) %>% 
  summarise_at(vars(contains("lb_pval")), median)
  #group_by(length, prst) %>% 
  #slice_min(value_aic)

tt = tt %>% mutate(indep_flag = lb_flag + lb_abs_flag + lb_sq_flag)
tt %>% pull(indep_flag) %>% table()

tt %>% group_by(length, prst) %>% 
  summarise_at(vars(contains("lb_sq_pval")), min)

tt %>% 
  filter(sw_flag == 0) %>% 
  filter(jb_flag == 0) %>% 
  filter(lb_flag == 0) %>% 
  filter(lb_sq_flag == 0) %>% 
  filter(lb_abs_flag == 0) %>% 
  arrange(value_aic)

tt_full %>% filter(nr==736) %>% 
  mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~irf_whf(.x, .y, n_lags = 48))) %>% 
  .$irf %>% .[[1]] %>% plot
