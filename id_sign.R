# -------------------------------------------------------------- #
# This script generates the output regarding the main simulation #
# exercise. The simulation results are loaded from a data frame  #
# containing the necessary components for this script to run and #
# make the figures and tables in the paper. -------------------- #
# -------------------------------------------------------------- #

pkgs <- c("svarmawhf", "tidyverse", "lemon")
void = lapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE)))
source("list_of_functions.R")
THRESHOLD_SW = 0.05
THRESHOLD_JB = 0.05
THRESHOLD_LB = 0.01
THRESHOLD_LB_AUX = 0.01
THRESHOLD_C = 0.05
N_LB_LAG = 24

tt_empex <- readRDS("~/Documents/Rscripts/sim_news/local_data/tt_empex_20231024.rds")
tt_irf <- tt_empex %>%
  mutate(lbl = names(data_list)) %>% 
  # Calculate shocks
  mutate(res = pmap(., pmap_get_residuals_once)) %>% 
  mutate(B_mat = map2(.x = params_deep_final,
                      .y = tmpl,
                      ~fill_tmpl_whf_rev(theta = .x,
                                         tmpl = .y)$B
                      )
         ) %>%
  # mutate(B_mat = map(.x = res, ~t(chol(cov(.x))))) %>%
  mutate(shock = map2(.x = res,
                      .y = B_mat,
                      ~ t(solve(.y, t(.x)))
                      )
         ) %>%
  # Test of normality pt. 1: Shapiro-Wilk
  mutate(sw = map(.x = shock,
                  ~ apply(X = .x, 
                          MARGIN = 2, 
                          FUN = function(x){shapiro.test(x)$p.value}
                          )
                  )
         ) %>% 
  mutate(sw_flag = map_lgl(sw, ~any(.x > THRESHOLD_SW))) %>% 
  # Test of normality pt. 2: Jarque-Bera
  mutate(jb = map(.x = shock,
                  ~apply(X = .x,
                         MARGIN = 2, 
                         FUN = function(x){tsoutliers::JarqueBera.test(x)[[1]]$p.value}
                         )
                  )
         ) %>% 
  mutate(jb_flag = map_lgl(.x = jb, ~ any(.x > THRESHOLD_JB))) %>%
  # Test of serial correlation: Level
  mutate(lb = map(.x = shock,
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
  mutate(lb_sq = map(.x = shock,
                     ~ apply(X = .x, 
                             MARGIN = 2, 
                             FUN = function(x){ Box.test(x = x^2,
                                                         lag = N_LB_LAG, 
                                                         type = "Ljung-Box")$p.value}
                             )
                     )
         ) %>% 
  mutate(lb_sq_flag = map_lgl(lb_sq, ~any(.x < THRESHOLD_LB_AUX))) %>% 
  # Test of mutual rank correlation in X
  mutate(cor_p = map(.x = shock,
                     ~ apply(
                       X = combn(1:ncol(.x), 2),
                       MARGIN = 2,
                       FUN =
                         function(x)
                           cor.test(x = .x[,x[1]],
                                    y = .x[,x[2]],
                                    method = "pearson")$p.value
                       )
                     )
         ) %>%
  mutate(cor_p_flag = map_lgl(cor_p, ~ any(.x < THRESHOLD_C))) %>%
  # # Test of mutual linear correlation in X
  mutate(cor_s = map(.x = shock,
                     ~ apply(
                       X = combn(1:ncol(.x), 2),
                       MARGIN = 2,
                       FUN =
                         function(x)
                           cor.test(x = .x[,x[1]],
                                    y = .x[,x[2]],
                                    method = "spearman",
                                    exact = FALSE)$p.value
                       )
                     )
         ) %>%
  mutate(cor_s_flag = map_lgl(cor_s, ~ any(.x < THRESHOLD_C))) %>%
  # # Test of mutual linear correlation in X^2
  mutate(cor_sq_p = map(.x = shock,
                        ~ apply(
                          X = combn(1:ncol(.x), 2),
                          MARGIN = 2,
                          FUN =
                            function(x)
                              cor.test(x = .x[,x[1]]^2,
                                       y = .x[,x[2]]^2,
                                       method = "pearson")$p.value
                          )
                        )
         ) %>%
  mutate(cor_sq_p_flag = map_lgl(cor_sq_p, ~ any(.x < THRESHOLD_C))) %>%
  # # Test of mutual rank correlation in X^2
  mutate(cor_sq_s = map(.x = shock,
                        ~ apply(
                          X = combn(1:ncol(.x), 2),
                          MARGIN = 2,
                          FUN = function(x)
                            cor.test(x = .x[,x[1]]^2,
                                     y = .x[,x[2]]^2,
                                     method = "spearman",
                                     exact = FALSE)$p.value
                          )
                        )
         ) %>%
  mutate(cor_sq_s_flag = map_lgl(cor_sq_s, ~ any(.x < THRESHOLD_C))) %>%
  # Model diagnostics checks
  mutate(norm_flag = rowSums(across(contains("jb") & contains("sw") & contains("flag")))) %>%
  mutate(indep_flag = rowSums(across(contains("lb") & contains("flag")))) %>%
  mutate(cor_flag = rowSums(across(contains("cor") & contains("flag")))) %>%
  mutate(norm_indep_flag = norm_flag + indep_flag + cor_flag)

tt_irf <- tt_irf %>% 
  filter(norm_indep_flag==0) %>%
  mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~irf_whf(.x,.y,48))) %>% 
  # filter gaussian processes out
  # choose best model per shock and distr
  # aic, bic
  # mutate(nobs = map_dbl(data_list, ~nrow(.x))) %>%
  # mutate(n_params = map_dbl(.x = params_deep_final,
  #                           ~length(.x))) %>% 
  # mutate(punish_aic = n_params * 2/nobs) %>% 
  # mutate(punish_bic = n_params * log(nobs)/nobs) %>% 
  # mutate(value_aic = value_final + punish_aic) %>% 
  # mutate(value_bic = value_final + punish_bic) %>%
  # group_by(shock_distr, lbl) %>%
  # mutate(rk_aic = rank(value_aic),
  #        rk_bic = rank(value_bic),
  #        rk_mle = rank(value_final)) %>% 
  # slice_min(rk_aic) %>%
  # ungroup() %>%
  # mutate(armamod = map2(.x = params_deep_final, .y = tmpl, ~armamod_whf(.x, .y))) %>% 
  # mutate(irf = map(.x = armamod, ~pseries(lmfd(.x$polm_ar, .x$polm_ma), 48))) %>%
  # mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~ irf_whf(.x, .y, 48))) %>% 
  # mutate(irf = map(.x = irf, ~ .x[1:2,,] %>% 
  #                    apply(c(1,2), cumsum) %>%
  #                    aperm(c(2,3,1)) %>%  
  #                    abind::abind(.x[3:4,,], 
  #                                 along = 1)
  #                  )
  #        ) %>% 
  mutate(irf = map2(.x = sd_vec, .y = irf, ~ diag(.x)%r%.y))

rest_hor <- 6
sgn_mat <- array(matrix(NA, 4, 4), c(4, 4, rest_hor))
sgn_mat[,1,1] <- c(1,1,1,NA) # Demand shock 
sgn_mat[,2,1] <- c(-1,1,1,NA) # Supply shock
sgn_mat[2,3,] <- -1 # MP -> CPI
sgn_mat[3,3,] <- 1 # MP -> FFR
sgn_mat[4,3,1] <- 1 # MP -> MPR
#sgn_mat[2,4,] <- -1 # FG -> CPI
sgn_mat[3,4,1] <- 0 # FG -> FFR

# param_list <- list(chol_irf = unclass(tt_tmp$chol_irf[[7]]),
#                    rest_mat = sgn_mat,
#                    news_rest = list(1, c(3,4)),
#                    ndraws = 1e3,
#                    verbose = TRUE)

param_list <- pmap(list(x = tt_irf$irf,
                        y = tt_irf$res,
                        z = tt_irf$B_mat),
                   function(x, y, z) 
                     list(irf_arr = unclass(x),
                          emp_innov = y,
                          w_mat = z,
                          rest_mat = sgn_mat,
                          news_rest = c(1, 3, 4),
                          rest_mom = list(type = "kurt",
                                          test_shocks = 3:4,
                                          threshold = 1.5),
                          ndraws = 1e2,
                          max_draws = 1e4,
                          irf_cb = c(0.16, 0.84),
                          replace_md = TRUE,
                          verbose = TRUE)
                   )

ncores <- parallel::detectCores()-1
par_fun <- function(param_list){
  source("list_of_functions.R")
  do.call(id_mixed_new, param_list)
}
# param_list$ndraws <- ceiling(1e2/3)
# param_list$max_draws <- ceiling(1e4/3)
cl <- parallel::makeCluster(ncores)
parallel::clusterExport(cl, ls())
parallel::clusterEvalQ(cl, library("tidyverse"))
parallel::clusterEvalQ(cl, library("svarmawhf"))
test <- parallel::parLapplyLB(cl = cl, X = param_list, fun = par_fun)
parallel::stopCluster(cl)

irf_arr <- lapply(test[!ok], function(x) x$irf) %>% abind::abind(along=4)

