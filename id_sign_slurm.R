# -------------------------------------------------------------- #
# This script generates the output regarding the main simulation #
# exercise. The simulation results are loaded from a data frame  #
# containing the necessary components for this script to run and #
# make the figures and tables in the paper. -------------------- #
# -------------------------------------------------------------- #

source("/home/juhokois/proj/sim_news/list_of_functions.R")
pkgs <- c("svarmawhf", "tidyverse")
void = lapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE)))

# Arguments from Rscript call: Parameters from SLURM script ####
args = commandArgs(trailingOnly=TRUE)
params = list()

# Parameters from Rmarkdown
params$NEW_DIR = args[1]
  
THRESHOLD_SW = 0.05
THRESHOLD_JB = 0.05
THRESHOLD_LB = 0.05
THRESHOLD_LB_AUX = 0.01
N_LB_LAG = 24

tt_empex <- readRDS("/home/juhokois/proj/sim_news/local_data/tt_empex_20230907.rds")
tt_irf <- tt_empex %>%
  filter(abs(value_final)<100) %>%
  # Extract residuals
  mutate(res = pmap(., pmap_get_residuals_once)) %>%
  # Test of normality pt. 1: Shapiro-Wilk
  mutate(sw = map(.x = res,
                  ~ apply(X = .x, 
                          MARGIN = 2, 
                          FUN = function(x){shapiro.test(x)$p.value}
                          )
                  )
         ) %>% 
  mutate(sw_flag = map_lgl(sw, ~any(.x > THRESHOLD_SW))) %>% 
  # Test of normality pt. 2: Jarque-Bera
  mutate(jb = map(.x = res,
                  ~apply(X = .x,
                         MARGIN = 2, 
                         FUN = function(x){tsoutliers::JarqueBera.test(x)[[1]]$p.value}
                         )
                  )
         ) %>% 
  mutate(jb_flag = map_lgl(.x = jb, ~ any(.x > THRESHOLD_JB))) %>%
  mutate(norm_flag = rowSums(across(contains("flag")))) %>% 
  # Test of serial correlation: Level
  mutate(lb = map(.x = res,
                  ~ apply(X = .x,
                          MARGIN = 2, 
                          FUN = function(x){Box.test(x = x, 
                                                     lag = N_LB_LAG, 
                                                     type = "Ljung-Box")$p.value}
                          )
                  )
         ) %>% 
  mutate(lb_flag = map_lgl(lb, ~ any(.x < THRESHOLD_LB))) %>% 
  # Test of serial correlation: Shocks squared
  mutate(lb_sq = map(.x = res,
                     ~ apply(X = .x, 
                             MARGIN = 2, 
                             FUN = function(x){ Box.test(x = x^2,
                                                         lag = N_LB_LAG, 
                                                         type = "Ljung-Box")$p.value}
                             )
                     )
         ) %>% 
  mutate(lb_sq_flag = map_lgl(lb_sq, ~any(.x < THRESHOLD_LB_AUX))) %>% 
  mutate(indep_flag = rowSums(across(contains(c("lb_flag", "lb_sq_flag"))))) %>% 
  mutate(norm_indep_flag = norm_flag+indep_flag) %>% 
  mutate(nobs = map_dbl(data_list, ~nrow(.x))) %>%
  mutate(n_params = map_dbl(.x = params_deep_final,
                            ~length(.x))) %>% 
  mutate(punish_aic = n_params * 2/nobs) %>% 
  mutate(punish_bic = n_params * log(nobs)/nobs) %>% 
  mutate(value_aic = value_final + punish_aic) %>% 
  mutate(value_bic = value_final + punish_bic) %>%
  mutate(rk_aic = rank(value_aic),
         rk_bic = rank(value_bic),
         rk_mle = rank(value_final)) %>% 
  mutate(armamod = map2(.x = params_deep_final, .y = tmpl, ~armamod_whf(.x, .y))) %>% 
  mutate(chol_fac = map(.x = res, ~ t(chol(cov(.x))))) %>% 
  mutate(irf = map(.x = armamod, ~pseries(lmfd(.x$polm_ar, .x$polm_ma), 48))) %>%
  mutate(chol_irf = map2(.x = irf, .y = chol_fac, ~.x%r%.y))

tt_tmp <- tt_irf %>%
  filter(norm_indep_flag==0) %>% 
  mutate(lbl = names(data_list)) %>% 
  group_by(lbl, shock_distr) %>% 
  slice_min(rk_aic) %>% 
  ungroup()
  
rest_hor <- 3 # sign restrictions holding up to rest_hor - 1 months
sgn_mat <- array(matrix(NA, 5, 5), c(5, 5, rest_hor))
sgn_mat[2,4:5,] <- -1 # CPI
sgn_mat[3,4:5,] <- -1 # EBP
sgn_mat[4,4,] <- 1 # FFR
sgn_mat[5,4,1] <- 1 # MPR
sgn_mat[5,5,1] <- 0

map(.x = tt_tmp$chol_irf, ~ list(chol_irf = unclass(.x),
                                 sign_mat = sgn_mat,
                                 ndraws = 10,
                                 max_draws = 1e4,
                                 verbose = FALSE)) -> param_list

id_obj <- lapply(param_list, function(x) do.call(id_mixed_new, x))

tibble_id <- paste0("/tibble_",
                    paste(sample(0:9, 5, replace = TRUE), collapse = ""), 
                    paste(sample(letters, 5, replace = TRUE), collapse = ""),
                    ".rds")
tibble_out <- tt_tmp %>% 
  dplyr::select(p, q, kappa, k, n_unst, irf, value_final, value_aic, value_bic) %>% 
  mutate(id_obj)

#Create new directory from slurm and save id_objects there
saveRDS(tibble_out, file = paste0(params$NEW_DIR, tibble_id))