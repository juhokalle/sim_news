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
THRESHOLD_LB_AUX = 0.05
THRESHOLD_C = 0.05
N_LB_LAG = 24

tt_empex <- readRDS("~/Documents/Rscripts/sim_news/local_data/tt_empex_20231128.rds")
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
                             FUN = function(x){Box.test(x = x^2,
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
  mutate(norm_flag = sw_flag + jb_flag) %>%
  mutate(indep_flag = lb_flag + lb_sq_flag) %>%
  mutate(cor_flag = rowSums(across(contains("cor") & contains("flag")))) %>%
  mutate(norm_indep_flag = norm_flag + indep_flag + cor_flag)
  
tt_est <- tt_irf %>% 
    # filter(norm_flag == 0,
    #        indep_flag == 0,
    #        cor_p_flag == 0,
    #        cor_s_flag == 0,
    #        lbl == "GSS22",
    #        shock_distr == "skewed_t") %>%
    filter(n_unst==1, q==1, shock_distr =="skewed_t") %>% 
    # filter gaussian processes out
    # choose best model per shock and distr: aic, bic
    # ----------------------- #
    mutate(nobs = map_dbl(data_list, ~nrow(.x))) %>%
    mutate(n_params = map_dbl(.x = params_deep_final, ~length(.x))) %>%
    mutate(punish_aic = n_params * 2/nobs) %>%
    mutate(punish_bic = n_params * log(nobs)/nobs) %>%
    mutate(value_aic = value_final + punish_aic) %>%
    mutate(value_bic = value_final + punish_bic) %>%
    group_by(shock_distr, lbl) %>%
    mutate(rk_aic = rank(value_aic),
           rk_bic = rank(value_bic),
           rk_mle = rank(value_final)) %>%
    slice_min(value_bic) %>%
    ungroup() %>%
    # ----------------------- #
    mutate(armamod = map2(.x = params_deep_final, .y = tmpl, ~armamod_whf(.x, .y))) %>%
    mutate(irf = map(.x = armamod, ~pseries(lmfd(.x$polm_ar, .x$polm_ma), 48))) %>%
    mutate(irf = map2(.x = sd_vec, .y = irf, ~ diag(.x)%r%.y))

rest_hor <- 6
dim_out <- 5
sgn_mat <- array(NA, c(dim_out, dim_out, rest_hor))
sgn_mat[,1,] <- c(1, 1, NA, 1, NA) # Demand shock 
sgn_mat[,2,] <- c(-1, 1, NA, 1, NA) # Supply shock
sgn_mat[,3,1] <- c(NA, NA, 1, NA, NA) # Demand shock
sgn_mat[,4,1] <- c(NA, -1, -1, 1, 1) # MP impact
sgn_mat[,4,-1] <- c(NA, -1, -1, 1, NA) # MP dynamic
sgn_mat[,5,1] <- c(NA, -1, -1,  NA, 0)  # FG impact
sgn_mat[,5,-1] <- c(NA, -1, -1, NA, NA)  # FG dynamic

param_list <- pmap(with(tt_est, list(x = irf, y = res, z = B_mat)),
                   function(x, y, z) 
                     list(irf_arr = unclass(x)[,,1:37],
                          emp_innov = y,
                          w_mat = z,
                          rest_mat = sgn_mat,
                          news_rest = c(1, 4, 5),
                          rest_mom = list(type = "kurt",
                                          test_shocks = 4:5,
                                          threshold = 1.5),
                          ndraws = 50, #ceiling(1000/14),
                          sub_max = 100,
                          max_draws = 1e6, #ceiling(1e6/14), 
                          irf_cb = NULL, #c(.05,.95),
                          replace_md = FALSE,
                          verbose = TRUE)
                   )

system.time( test1 <- do.call(id_mixed_new, param_list[[1]]) )

ncores <- parallel::detectCores()-1
par_fun <- function(param_list){
  source("list_of_functions.R")
  do.call(id_mixed_new, param_list)
}

# param_list$ndraws <- ceiling(1e2/3)
# param_list$max_draws <- ceiling(1e4/3)
cl <- parallel::makeCluster(ncores)
#parallel::clusterExport(cl, ls())
parallel::clusterEvalQ(cl, library("tidyverse"))
parallel::clusterEvalQ(cl, library("svarmawhf"))
test <- parallel::parLapplyLB(cl = cl, X = rep(param_list, 14), fun = par_fun)
parallel::stopCluster(cl)

# create and save irfs
irf_fnl <- abind::abind(lapply(test, function(x) x$irf),
                        along = 4)
shocks_fnl <- abind::abind(lapply(test, function(x) x$shocks),
                           along = 3)
for(j in 1:dim(irf_fnl)[4]){
  shock_sd <- diag(apply(shocks_fnl[,,j], 2, sd)^-1)
  irf_fnl[,,,j] <- apply(irf_fnl[,,,j], 3, function(x) x%*%shock_sd)
}


irf_fnl <- irf_fnl[,,1:37,]
irf_qt <- apply(irf_fnl, c(1,2,3), quantile, probs = c(.05,.5,.95)) %>% 
  aperm(c(2,3,4,1))
irf_qt[,,,2] <- fry_pagan_mt(irf_arr = irf_fnl)

tt_final <- expand_grid(type = c("lb", "true", "ub"),
              horizon = 0:36,
              shock = factor(c("MP", "FG"), levels = c("MP", "FG")), 
              var = factor(c("IP", "CPI", "SP500", "FFR", "MPR"), 
                           levels= c("IP", "CPI", "SP500", "FFR", "MPR"))
              ) %>% 
  bind_cols(irf_value = c(irf_qt[,4:5,,]))



plt1 <- tt_final %>% 
  pivot_wider(names_from = type, 
              values_from = irf_value) %>%
  ggplot(aes(x=horizon)) + 
  geom_line(aes(y = true)) +
  geom_line(aes(y = lb), linetype = "dashed", size = 0.25) +
  geom_line(aes(y = ub), linetype = "dashed", size = 0.25) +
  geom_hline(yintercept = 0, size = .25) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.1, fill = "blue") + 
  theme(panel.grid.minor = element_blank()) +
  ylab("") + 
  xlab("") +
  facet_grid(var ~ shock, 
             scales = "free_y"
             ) +
  scale_x_continuous(breaks = seq(0, 36, by = 6), expand = c(0.025,0.025))

ggsave("paper_output/MAIN_IRF2.pdf", plt1)