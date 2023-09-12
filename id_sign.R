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
THRESHOLD_LB = 0.05
THRESHOLD_LB_AUX = 0.01
N_LB_LAG = 24

tt_empex <- readRDS("~/Documents/Rscripts/sim_news/local_data/tt_empex_20230907.rds")
tt_irf <- tt_empex %>%
  # Calculate shocks
  filter(abs(value_final)<100) %>% 
  mutate(res = pmap(., pmap_get_residuals_once)) %>% 
  # mutate(B_mat = map2(.x = params_deep_final, 
  #                     .y = tmpl, 
  #                     ~fill_tmpl_whf_rev(theta = .x, 
  #                                        tmpl = .y)$B
  #                     )
  #        ) %>% 
  # mutate(res = map2(.x = res, 
  #                      .y = B_mat, 
  #                      ~ solve(.y, t(.x)) %>% t()
  #                      )
  #        ) %>%
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
  filter(norm_indep_flag==0) %>% 
  # filter gaussian processes out
  # choose best model per shock and distr
  # aic, bic
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
  group_by(shock_distr, data_list) %>%
  slice_min(rk_aic) %>%
  ungroup() %>%
  mutate(chol_fac = map(.x = res, ~ t(chol(cov(.x))))) %>% 
  mutate(armamod = map2(.x = params_deep_final, .y = tmpl, ~armamod_whf(.x, .y))) %>% 
  mutate(irf = map(.x = armamod, ~pseries(lmfd(.x$polm_ar, .x$polm_ma), 48))) %>%
  mutate(chol_irf = map2(.x = irf, .y = chol_fac, ~.x%r%.y))

tt_tmp <- tt_irf %>%
  filter(norm_indep_flag==0) %>% 
  mutate(lbl = names(data_list)) %>% 
  group_by(lbl, shock_distr) %>% 
  slice_min(rk_aic) %>% 
  ungroup()
  
rest_hor <- 6
sgn_mat <- array(matrix(NA, 5, 5), c(5, 5, rest_hor))
sgn_mat[2,4,] <- -1 # MP,FG -> CPI
sgn_mat[3,4,] <- -1 # MP,FG -> EBP
sgn_mat[4,4,] <- 1 # MP -> FFR
sgn_mat[5,4,1] <- 1 # MP -> MPR
sgn_mat[5,5,1] <- 0 # FG -> MPR

ncores <- parallel::detectCores()-1
cl <- parallel::makeCluster(ncores)
par_fun <- function(param_list){
  do.call(id_mixed_new, param_list)
}
param_list <- list(chol_irf = unclass(tt_tmp$chol_irf[[2]]),
                   rest_mat = sgn_mat,
                   ndraws = ceiling(1e3/ncores/6),
                   max_draws = ceiling(1e6/ncores/6),
                   verbose = FALSE)
# map(.x = tt_tmp$chol_irf, ~ list(chol_irf = unclass(.x),
#                                  ndraws = 1e2,
#                                  max_draws = 1e5,
#                                  sign_mat = sgn_mat,
#                                  policy_var = 3,
#                                  mp_hor = 24,
#                                  zr_ix = c(3,4),
#                                  verbose = TRUE)) -> param_list

parallel::clusterExport(cl, ls())
parallel::clusterEvalQ(cl, library("svarmawhf"))
parallel::clusterEvalQ(cl, library("tidyverse"))
system.time(test <- parallel::parLapplyLB(cl = cl, X = rep(list(param_list), ncores*6), fun = par_fun))
parallel::stopCluster(cl)

res_list <- list()
res_list$tot_ok <- sapply(test, function(x) x$draws_ok)
res_list$draws_total <- sapply(test, function(x) x$draws_total/ceiling(1e6/42))
res_list$struc_irf <- abind::abind(lapply(test[as.logical(res_list$tot_ok)], function(x) x$irf), 
                                   along = 4)
aux_vec <- apply(res_list$struc_irf[4,2,,], 2, function(x) ifelse(max(abs(x))==max(x), 1, -1))
# 
for(j in 1:dim(res_list$struc_irf)[4]){
  res_list$struc_irf[, , , j] <- apply(res_list$struc_irf[, , , j], 3, function(x) x%*%diag(c(1, aux_vec[j])))
  # test$bmat[,,j] <- test$bmat[,,j]%*%aa[[j]]
}

apply(res_list$struc_irf, c(1,2,3), quantile, probs = c(.05, .14, .5, .86, .95)) %>% 
  aperm(c(2,3,4,1)) -> irf_qt

plot_irf(irf_arr = irf_qt, 
         var_name = c("IP", "CPI", "EBP", "FFR", "MPR"), 
         shock_name = c("MP shock", "FG shock"))
print(paste0("Share of draws retained: ", 
             round(
               100*dim(struc_irf)[4]/sum(sapply(test, function(x) x$draws_total)), 2
               ),
             "%"
             )
      )
