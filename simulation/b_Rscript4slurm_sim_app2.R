# ------------------------------------------------------------- #
# Script to be called via SLURM ------------------------------- #
# ------------------------------------------------------------- #
# This is the main simulation script -------------------------- #
# 1) Simulate data from the underlying structural model ------- #
# 2) Estimate model with the pre-specified model specifications #
# 3) Calculate the IRF and other objects of interest ---------- #
# 4) Save the resulting tibble with random named file name ---- #
# ------------------------------------------------------------- #

# PREAMBLE ####
source("/proj/juhokois/sim_news/list_of_functions.R")
.libPaths(c("/proj/juhokois/R/", .libPaths()))
pkgs <- c("svarmawhf", "fitdistrplus", "sgt", "tidyverse")
void = lapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE)))

# Arguments from Rscript call: Parameters from SLURM script ####
args = commandArgs(trailingOnly=TRUE)
params = list()

# Parameters from Rmarkdown
params$IX_ARRAY_JOB = as.integer(args[1]) # index of array-job. Number of array-jobs is determined from number of rows of dataframe containing all integer-values parameters
params$SLURM_JOB_ID = as.integer(args[2])
params$MANUALLY_ASSIGNED_ID = as.integer(args[3])
params$SLURM_ARRAY_TASK_MAX = as.integer(args[4])
params$NEW_DIR = args[5]

# OPTIMIZATION PARAMS

## general
params$RESTART_W_NOISE = 1
params$FIX_INIT = FALSE
params$IC <- TRUE
params$penalty_prm = 100

## gaussian density
params$IT_OPTIM_GAUSS = 1
params$USE_BFGS_GAUSS = TRUE
params$USE_NM_GAUSS = TRUE
params$USE_CS_GAUSS = FALSE
params$MAXIT_BFGS_GAUSS = 200
params$MAXIT_NM_GAUSS = 3000
params$MAXIT_CS_GAUSS = 500

## laplacian density
params$IT_OPTIM_LAPLACE = 3
params$USE_BFGS_LAPLACE = TRUE
params$USE_NM_LAPLACE = TRUE
params$USE_CS_LAPLACE = FALSE
params$MAXIT_BFGS_LAPLACE = 200 # default for derivative based methods
params$MAXIT_NM_LAPLACE = 3000 # default for NM is 500
params$MAXIT_CS_LAPLACE = 500

## sgt density
params$IT_OPTIM_SGT = 3
params$USE_BFGS_SGT = TRUE
params$USE_NM_SGT = TRUE
params$USE_CS_SGT = FALSE
params$MAXIT_BFGS_SGT = 200 # default for derivative based methods
params$MAXIT_NM_SGT = 3000 # default for NM is 500
params$MAXIT_CS_SGT = 500

# SIMULATION SPECS: MODEL PARAMS
sim_prm <- expand_grid(beta = c(0.5, 0.9), 
                       rho = 0.5,
                       nobs = 250,
                       nu = c(3, 20))
sim_prm <- sim_prm[(params$IX_ARRAY_JOB-1)%/%(params$SLURM_ARRAY_TASK_MAX/nrow(sim_prm))+1, ]
DIM_OUT <- 2
params$DIM_OUT = DIM_OUT

# GENERATE DATA
sim_obj <- do.call(what = sim_news, args = sim_prm)

# Tibble with integer-valued parameters
tt =
  # orders (p,q)
  expand_grid(p = 1:2, q = 0:2) %>% 
  # number of unstable zeros
  mutate(n_unst = map(q, ~0:(DIM_OUT*.x))) %>% 
  unnest(n_unst) %>% 
  mutate(n_st = DIM_OUT * q - n_unst) %>% 
  mutate(kappa = n_unst %/% DIM_OUT,
         k = n_unst %% DIM_OUT) %>% 
  # assume correct specification w.r.t. MA polm
  # VAR benchmark
  bind_rows(c(p = 12, q = 0, n_unst = 0, n_st = 0, kappa = 0, k = 0)) %>% 
  # this way of including data makes it convenient for slicing
  expand_grid(sim_prm) %>%
  mutate(data_list = list(sim_obj$y$y)) %>% 
  mutate(data_list = map(.x = data_list, ~ apply(X = .x, 
                                                 MARGIN = 2, 
                                                 FUN = function(x) x-mean(x)
                                                 )
                         )
         ) %>% 
  # template
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  # generate initial values and likelihood functions (we can use the same template for initial values and likelihood fct bc both have no parameters for density)
  mutate(theta_init = map2(.x = data_list, .y = tmpl, ~get_init_armamod_whf_random(.x, .y))) %>%
  mutate(shock_distr = "tdist") %>%
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev))

# Parallel setup ####
tt_optim_parallel = tt %>% 
  dplyr::select(theta_init, tmpl, data_list, shock_distr)

params_parallel = lapply(1:nrow(tt_optim_parallel),
                         function(i) t(tt_optim_parallel)[,i])
mods_parallel_list <- lapply(params_parallel, hlp_parallel)

tibble_out =  
  enframe(mods_parallel_list) %>% 
  unnest_wider(value) %>% 
  unnest_wider(results_list) %>%
  unnest_wider(input_integerparams) %>% 
  mutate(tt) %>%
  mutate(n_params = map_dbl(.x = params_deep_final,
                            ~ length(.x)
                            )
         ) %>% 
  mutate(punish_aic = n_params * 2/nobs) %>% 
  mutate(punish_bic = n_params * log(nobs)/nobs) %>% 
  mutate(value_bic = value_final + punish_bic) %>%
  mutate(svar_bool = map_lgl(.x = p, ~p==12)) %>% 
  group_by(svar_bool) %>%
  slice_min(value_bic) %>%
  ungroup() %>%
  mutate(irf = map2(.x = params_deep_final, 
                    .y = tmpl, 
                    ~ irf_whf(.x, .y, n_lags = 8)
                    )
         ) %>% 
  mutate(true_irf = map(.x = beta,
                        ~with(sim_news(beta = .x, 
                                       rho = 0.5, 
                                       no_sim = TRUE),
                              pseries(sys, 8)%r%sigma_L
                              )
                        )
         ) %>% 
  mutate(rmat = map2(.x = true_irf, 
                     .y = irf, 
                     ~ choose_perm_sign(target_mat = .x,
                                        cand_mat = .y,
                                        type = "min_rmse")
                     )
         ) %>% 
  mutate(rmat = map2(.x = irf, 
                     .y = rmat, 
                     ~ .y%*%optim_zr(input_mat = unclass(.x)[,,1]%*%.y,
                                     zr_ix = c(1,2),
                                     opt_it = FALSE)
                     )
         ) %>%
  mutate(irf = map2(.x = irf,
                    .y = rmat,
                    ~ .x%r%.y
                    )
         ) %>%
  mutate(mad = map2_dbl(.x = irf,
                        .y = true_irf,
                        ~ mean(abs(unclass(.x-.y)))
                        )
         ) %>% 
  dplyr::select(beta, nu, p, nobs, mad, irf)

tibble_id <- paste0("/tibble_",
                    paste(sample(0:9, 5, replace = TRUE), collapse = ""), 
                    paste(sample(letters, 5), collapse = ""),
                    ".rds")

saveRDS(tibble_out, file = paste0(params$NEW_DIR, tibble_id))