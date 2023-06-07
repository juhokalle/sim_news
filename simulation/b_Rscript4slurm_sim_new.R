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
params$SRUN_CPUS_PER_TASK = as.integer(args[6])

# OPTIMIZATION PARAMS

## general
params$RESTART_W_NOISE = 1
params$PERM_INIT = 5
params$FIX_INIT = FALSE
params$IC <- TRUE
params$penalty_prm = 1000
params$PATH_RESULTS_HELPER = "/proj/juhokois/sim_news/local_data/"

## gaussian density
params$IT_OPTIM_GAUSS = 2
params$USE_BFGS_GAUSS = TRUE
params$USE_NM_GAUSS = TRUE
params$USE_CS_GAUSS = FALSE
params$MAXIT_BFGS_GAUSS = 100
params$MAXIT_NM_GAUSS = 500
params$MAXIT_CS_GAUSS = 500

## laplacian density
params$IT_OPTIM_LAPLACE = 2
params$USE_BFGS_LAPLACE = TRUE
params$USE_NM_LAPLACE = TRUE
params$USE_CS_LAPLACE = FALSE
params$MAXIT_BFGS_LAPLACE = 100 # default for derivative based methods
params$MAXIT_NM_LAPLACE = 500 # default for NM is 500
params$MAXIT_CS_LAPLACE = 1000

## sgt density
params$IT_OPTIM_SGT = 2
params$USE_BFGS_SGT = TRUE
params$USE_NM_SGT = TRUE
params$USE_CS_SGT = FALSE
params$MAXIT_BFGS_SGT = 100 # default for derivative based methods
params$MAXIT_NM_SGT = 500 # default for NM is 500
params$MAXIT_CS_SGT = 500

# SIMULATION SPECS: MODEL PARAMS
sim_prm <- expand_grid(beta = c(0.5, 0.9), 
                       rho = 0.5, 
                       nobs = 250, 
                       nu = c(3, 20))
sim_prm <- sim_prm[(params$IX_ARRAY_JOB-1)%/%(params$SLURM_ARRAY_TASK_MAX/nrow(sim_prm))+1, ]
DIM_OUT <- 2
params$DIM_OUT = DIM_OUT

# SIMULATION SPECS: TARGET PATH FOR SAVING RESULTS
pap = pap_factory(params$PATH_RESULTS_HELPER)

# New directory for saving all tibbles
new_dir_path = pap(paste0("jobid_", params$MANUALLY_ASSIGNED_ID))

if(params$IX_ARRAY_JOB==1){
  if (!dir.exists(new_dir_path)){
    dir.create(new_dir_path)
  }
}

# GENERATE DATA
sim_obj <- map(1:params$SRUN_CPUS_PER_TASK, ~ list(sim_ix = .x, data_list = do.call(what = sim_news, args = sim_prm)$y$y))

# Tibble with integer-valued parameters
tt =
  # orders (p,q)
  expand_grid(p = 1, q = 2) %>% 
  # number of unstable zeros
  mutate(n_unst = map(q, ~0:(DIM_OUT*.x))) %>% 
  unnest(n_unst) %>% 
  mutate(n_st = DIM_OUT * q - n_unst) %>% 
  mutate(kappa = n_unst %/% DIM_OUT,
         k = n_unst %% DIM_OUT) %>% 
  # assume correct specification w.r.t. MA polm
  filter(n_unst==1) %>% 
  # VAR benchmark
  bind_rows(c(p = 12, q = 0, n_unst = 0, n_st = 0, kappa = 0, k = 0)) %>% 
  # this way of including data makes it convenient for slicing
  # expand_grid(sim_prm, data_list = sim_obj) %>%
  expand_grid(sim_prm) %>%
  mutate(data_list = list(sim_obj)) %>% 
  unnest_longer(data_list) %>% 
  unnest_wider(data_list) %>% 
  mutate(sd_vec = map(.x = data_list, ~ apply(.x, 2, sd))) %>% 
  mutate(data_list = map(.x = data_list, ~ apply(.x, 2, function(x) (x-mean(x))/sd(x)))) %>% 
  # template
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  # generate initial values and likelihood functions (we can use the same template for initial values and likelihood fct bc both have no parameters for density)
  mutate(ll_fun = map2(.x = data_list, .y = tmpl, ~ ll_whf_factory(t(unclass(.x)), .y, shock_distr = "gaussian"))) %>% 
  mutate(theta_init = map2(.x = data_list, .y = tmpl, ~get_init_armamod_whf_random(.x, .y))) %>%
  mutate(theta_init = map2(.x = theta_init, .y = ll_fun, ~ perm_init(.x, params$PERM_INIT, .y))) %>% 
  unnest_longer(theta_init) %>% 
  mutate(init_ix = rep(1:(params$PERM_INIT+1), n()/(params$PERM_INIT+1))) %>% 
  mutate(shock_distr = "tdist")

# Parallel setup ####
tt_optim_parallel = tt %>% 
  dplyr::select(theta_init, tmpl, data_list, shock_distr)

params_parallel = lapply(1:nrow(tt_optim_parallel),
                         function(i) t(tt_optim_parallel)[,i])
cl <- parallel::makeCluster(params$SRUN_CPUS_PER_TASK, "FORK", methods = FALSE)
mods_parallel_list = parallel::clusterApply(cl, params_parallel, hlp_parallel)
parallel::stopCluster(cl)

tibble_out =  
  enframe(mods_parallel_list) %>% 
  unnest_wider(value) %>% 
  unnest_wider(results_list) %>%
  unnest_wider(input_integerparams) %>% 
  mutate(tt) %>%
  group_by(q, sim_ix) %>%
  slice_min(value_final) %>%
  ungroup() %>%
  # mutate(res = pmap(., pmap_get_residuals_once)) %>% 
  # mutate(B_mat = map2(params_deep_final, tmpl, 
  #                     ~fill_tmpl_whf_rev(theta = .x, 
  #                                        tmpl = .y)$B)) %>% 
  #mutate(shocks = map2(res, B_mat, ~ solve(.y, t(.x)) %>% t())) %>%
  mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~ irf_whf(.x, .y, n_lags = 8))) %>% 
  mutate(irf = map2(.x = sd_vec, .y = irf, ~ diag(.x)%r%.y)) %>% 
  mutate(rmat = map(.x = irf, ~ id_news_shox(irf_arr = .x, policy_var = 1))) %>% 
  mutate(rmat = map2(.x = irf, .y = rmat, ~ .y%*%optim_zr(input_mat = unclass(.x)[,,1]%*%.y,
                                                          zr_ix = c(1,2),
                                                          opt_it = FALSE))) %>%
  mutate(irf = map2(.x = irf, .y = rmat, ~ .x%r%.y))

tibble_id <- paste0("/tibble_",
                    paste(sample(0:9, 5, replace = TRUE), collapse = ""), 
                    paste(sample(letters, 5), collapse = ""),
                    ".rds")

saveRDS(tibble_out, file = paste0(new_dir_path, tibble_id))
