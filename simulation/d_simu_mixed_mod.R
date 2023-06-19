# ----------------------------------------------------- #
# Script to be called via SLURM ----------------------- #
# ----------------------------------------------------- #
# This script carries out the simulation exercise ----- #
# showing that NG-SVARMA partially identifies --------- #
# a subset of shocks that are sufficiently non-Gaussian #
# 1) Simulate data from 3-eq. NK model with 2 different #
# distributional settings, purely NG and partially NG - #
# 2) Estimate under both setups ----------------------- #
# 3) Calculate IRFs and other objects of interest ----- #
# 4) Save the resulting tibble ------------------------ #
# ----------------------------------------------------- #

# PREAMBLE ####
source("/proj/juhokois/sim_news/list_of_functions.R")
.libPaths(c("/proj/juhokois/R/", .libPaths()))
pkgs <- c("svarmawhf", "fitdistrplus", "sgt", "tidyverse")
void = lapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE)))
set.seed(20230619)

# SIMU MODEL ####
DIM_OUT <- 3
n_obs <- 1e5
mod_mat <- get_struc_mat(model_type = "dynamic")
mod_mat$eps_val <- 10^-12
mod_str <- do.call(solve_re_mod_bp, mod_mat)
ma_polm <- abind::abind(diag(DIM_OUT),
                        with(mod_str, sigma_L%*%diag(c(0.5, 0.5, 2))%*%solve(sigma_L)), 
                        along = 3) %>% polm
dgp_mod <- armamod(sys = lmfd(mod_str$sys$a, ma_polm), # reduced-from varma 
                   sigma_L = mod_str$sigma_L) # with m0, sigma_L makes impact mat align with Bmat

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
params$PERM_INIT = 5
params$FIX_INIT = FALSE
params$IC <- TRUE
params$penalty_prm = 100
params$DIM_OUT = DIM_OUT

## gaussian density
params$IT_OPTIM_GAUSS = 1
params$USE_BFGS_GAUSS = FALSE
params$USE_NM_GAUSS = FALSE
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
params$MAXIT_NM_LAPLACE = 2000 # default for NM is 500
params$MAXIT_CS_LAPLACE = 500

## sgt density
params$IT_OPTIM_SGT = 3
params$USE_BFGS_SGT = TRUE
params$USE_NM_SGT = TRUE
params$USE_CS_SGT = FALSE
params$MAXIT_BFGS_SGT = 200 # default for derivative based methods
params$MAXIT_NM_SGT = 2000 # default for NM is 500
params$MAXIT_CS_SGT = 500

# GENERATE DATA
rg_fun <- 
  if(params$IX_ARRAY_JOB%%2){
    list(fun =  mixed_marg_dists(DIM_OUT, 3),
         lbl = "mixed")
  } else {
    list(fun = function(x) stats::rt(x, 3),
         lbl = "tdist")
  }

DATASET = do.call(what = simu_y, 
                  args = list(model = dgp_mod, 
                              n.obs = n_obs,
                              rand.gen = rg_fun$fun,
                              n.burnin = n_obs*2))$y

# MODEL SPECIFICATION
tt =
  # orders (p,q)
  expand_grid(p = 2, q = 1) %>% 
  # number of unstable zeros
  mutate(n_unst = map(q, ~0:(DIM_OUT*.x))) %>% 
  unnest(n_unst) %>% 
  mutate(n_st = DIM_OUT * q - n_unst) %>% 
  mutate(kappa = n_unst %/% DIM_OUT,
         k = n_unst %% DIM_OUT) %>% 
  # assume correct specification w.r.t. MA polm
  filter(n_unst==1) %>% 
  mutate(data_list = list(DATASET)) %>% 
  mutate(sd_vec = map(.x = data_list, ~ apply(.x, 2, sd))) %>% 
  mutate(data_list = map(.x = data_list, ~ apply(.x, 2, function(x) (x-mean(x))/sd(x)))) %>% 
  # template
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  # generate initial values and likelihood functions (we can use the same template for initial values and likelihood fct bc both have no parameters for density)
  mutate(theta_init = map2(.x = data_list, .y = tmpl, ~get_init_armamod_whf_random(.x, .y))) %>%
  mutate(theta_init = map2(.x = theta_init, .y = tmpl, ~ perm_init(.x, params$PERM_INIT, .y))) %>% 
  unnest_longer(theta_init) %>% 
  mutate(init_ix = rep(1:(params$PERM_INIT+1), n()/(params$PERM_INIT+1))) %>% 
  mutate(shock_distr = "tdist")

# Parallel setup ####
tt_optim_parallel = tt %>%
  slice((params$IX_ARRAY_JOB+1)%/%2) %>% 
  dplyr::select(theta_init, tmpl, data_list, shock_distr)

params_parallel = lapply(1:nrow(tt_optim_parallel),
                         function(i) t(tt_optim_parallel)[,i])

mods_parallel_list = lapply(params_parallel, FUN = hlp_parallel)

tibble_out =  
  enframe(mods_parallel_list) %>% 
  unnest_wider(value) %>% 
  unnest_wider(results_list) %>%
  unnest_wider(input_integerparams) %>% 
  mutate(tt %>% slice((params$IX_ARRAY_JOB+1)%/%2)) %>% 
  mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~ irf_whf(.x, .y, n_lags = 12))) %>% 
  mutate(irf = map2(.x = sd_vec, .y = irf, ~ diag(.x)%r%.y)) %>% 
  mutate(true_irf = list(with(dgp_mod, pseries(sys, 12)%r%sigma_L))) %>% 
  mutate(rmat = map2(.x = true_irf, .y = irf, ~ choose_perm_sign(target_mat = .x, 
                                                                 cand_mat = .y, 
                                                                 type = "min_rmse"))) %>%
  mutate(irf = map2(.x = irf, .y = rmat, ~ .x%r%.y)) %>% 
  expand_grid(rg = rg_fun$lbl)

tibble_id <- paste0("/tibble_",
                    paste(sample(0:9, 5, replace = TRUE), collapse = ""), 
                    paste(sample(letters, 5), collapse = ""),
                    ".rds")

saveRDS(tibble_out, file = paste0(params$NEW_DIR, tibble_id))
