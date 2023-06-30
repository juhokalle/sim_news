# ----------------------------------------------------------------------- #
# This script creates model output: diagnostics, model selection, IRFs... #
# ----------------------------------------------------------------------- #

source("list_of_functions.R")
pkgs = c("tidyverse", "svarmawhf")
void = lapply(pkgs, library, character.only = TRUE)
select <- dplyr::select
tt_simu <- readRDS("~/Documents/Rscripts/sim_news/local_data/tt_simu_main.rds")
tt <- tt_simu %>% 
  mutate(rmat = map(.x = irf, 
                    ~optim_zr(input_mat = unclass(.x)[,,1],
                              zr_ix = c(1,2),
                              opt_it = FALSE)
                    )
         ) %>% 
  mutate(irf = map2(.x = irf, 
                    .y = rmat, 
                    ~ .x%r%.y
                    )
         ) %>% 
  mutate(irf = map(.x = irf,
                   ~ .x %>% unclass
                   )
         ) %>%
  group_by(q, beta, nu) %>%
  # slice_sample(n=2000) %>% 
  reframe(irf_arr = list(abind::abind(irf, along=4))) %>% 
  mutate(irf_qt = map(.x = irf_arr,
                      ~ apply(X = .x, 
                              MARGIN = c(1,2,3), 
                              FUN = quantile, 
                              probs = c(.14, .86)
                              ) %>% 
                        aperm(c(2,3,4,1))
                      )
         ) %>% 
  mutate(true_irf = map(.x = beta,
                        ~with(sim_news(beta = .x, 
                                       rho = 0.5, 
                                       no_sim = TRUE),
                              pseries(sys, 8)%r%sigma_L
                              ) %>% unclass
                        )
         ) %>% 
  mutate(irf_qt = map2(.x = irf_qt, 
                       .y = true_irf,
                       ~abind::abind(.x, .y, along = 4)[,,,c(1,3,2)]
                       )
         )


THRESHOLD_SW = 0.05
THRESHOLD_JB = 0.05
THRESHOLD_LB = 0.01
THRESHOLD_C = 0.05
tt_simu <- tt_simu %>%
  # Calculate shocks
  mutate(res = pmap(., pmap_get_residuals_once)) %>% 
  mutate(B_mat = map2(.x = params_deep_final, 
                      .y = tmpl, 
                      ~fill_tmpl_whf_rev(theta = .x, 
                                         tmpl = .y)$B
  )
  ) %>% 
  mutate(shocks = map2(.x = res, 
                       .y = B_mat, 
                       ~ solve(.y, t(.x)) %>% t()
  )
  ) %>%
  # Test of normality pt. 1: Shapiro-Wilk
  mutate(sw = map(.x = shocks,
                  ~ apply(X = .x, 
                          MARGIN = 2, 
                          FUN = function(x){shapiro.test(x)$p.value}
                  )
  )
  ) %>% 
  mutate(sw_flag = map_lgl(sw, ~any(.x > THRESHOLD_SW))) %>% 
  # Test of normality pt. 2: Jarque-Bera
  mutate(jb = map(.x = shocks,
                  ~apply(X = .x,
                         MARGIN = 2, 
                         FUN = function(x){tsoutliers::JarqueBera.test(x)[[1]]$p.value}
                  )
  )
  ) %>% 
  mutate(jb_flag = map_lgl(.x = jb, ~ any(.x > THRESHOLD_JB))) %>%
  # Test of serial correlation: Level
  mutate(lb = map(.x = shocks,
                  ~ apply(X = .x,
                          MARGIN = 2, 
                          FUN = function(x){Box.test(x = x, 
                                                     lag = 24, 
                                                     type = "Ljung-Box")$p.value}
                  )
  )
  ) %>% 
  mutate(lb_flag = map_lgl(lb, ~ any(.x < THRESHOLD_LB))) %>% 
  # Test of serial correlation: Absolute values
  mutate(lb_abs = map(.x = shocks,
                      ~ apply(X = .x, 
                              MARGIN = 2, 
                              FUN = function(x){Box.test(x = abs(x), 
                                                         lag = 24, 
                                                         type = "Ljung-Box")$p.value}
                      )
  )
  ) %>% 
  mutate(lb_abs_flag = map_lgl(lb_abs, ~any(.x < THRESHOLD_LB))) %>% 
  # Test of serial correlation: Shocks squared
  mutate(lb_sq = map(.x = shocks, 
                     .y = policy_id, 
                     ~ apply(X = .x, 
                             MARGIN = 2, 
                             FUN = function(x){ Box.test(x = x^2,
                                                         lag = 24, 
                                                         type = "Ljung-Box")$p.value}
                     )
  )
  ) %>% 
  mutate(lb_sq_flag = map_lgl(lb_sq, ~any(.x < THRESHOLD_LB))) %>% 
  # Test of mutual rank correlation in X
  mutate(cor_p = map_dbl(.x = shocks,
                         ~ cor.test(x = .x[,1],
                                    y = .x[,2],
                                    method = "pearson")$p.value
  )
  ) %>%  
  mutate(cor_p_flag = map_lgl(cor_p, ~ .x < THRESHOLD_C)) %>% 
  # Test of mutual linear correlation in X
  mutate(cor_s = map_dbl(.x = shocks,
                         ~ cor.test(x = .x[,1],
                                    y = .x[,2],
                                    method = "spearman",
                                    exact = FALSE)$p.value
  )
  ) %>%
  mutate(cor_s_flag = map_lgl(cor_s, ~ .x < THRESHOLD_C)) %>% 
  # Test of mutual linear correlation in X^2
  mutate(cor_sq_p = map_dbl(.x = shocks,
                            .y = policy_id,
                            ~ cor.test(x = .x[,1]^2,
                                       y = .x[,2]^2,
                                       method = "pearson")$p.value
  )
  ) %>%  
  mutate(cor_sq_p_flag = map_lgl(cor_sq_p, ~ .x < THRESHOLD_C)) %>% 
  # Test of mutual rank correlation in X^2
  mutate(cor_sq_s = map_dbl(.x = shocks,
                            ~ cor.test(x = .x[,1]^2,
                                       y = .x[,2]^2,
                                       method = "spearman",
                                       exact = FALSE)$p.value
  )
  ) %>%
  mutate(cor_sq_s_flag = map_lgl(cor_sq_s, ~ .x < THRESHOLD_C)) %>% 
  mutate(norm_indep_flag = rowSums(across(contains("flag"))))

tt_simu %>% 
  group_by(p, beta, nu) %>% 
  summarise_at(vars(contains("flag")),
               ~100*mean(.)) %>%
  ungroup()


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