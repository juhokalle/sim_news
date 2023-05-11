# ----------------------------- #
# Script to be called via SLURM #
# ----------------------------- #

# PREAMBLE ####
source("/proj/juhokois/sim_news/list_of_functions.R")
.libPaths(c("/proj/juhokois/R/", .libPaths()))
pkgs <- c("svarmawhf", "fitdistrplus", "sgt", "tidyverse")
void = lapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE)))
nrep_est <- 20

# Arguments from Rscript call: Parameters from SLURM script ####
args = commandArgs(trailingOnly=TRUE)
params = list()

# Parameters from Rmarkdown
params$IX_ARRAY_JOB = as.integer(args[1]) # index of array-job. Number of array-jobs is determined from number of rows of dataframe containing all integer-values parameters
params$SLURM_JOB_ID = as.integer(args[2])
params$MANUALLY_ASSIGNED_ID = as.integer(args[3])
params$SLURM_ARRAY_TASK_MAX = as.integer(args[4])

# OPTIMIZATION PARAMS
params$RESTART_W_NOISE = 0
params$FIX_INIT = FALSE
params$IC <- TRUE
params$penalty_prm = 100

params$IT_OPTIM_GAUSS = 3
params$USE_BFGS_GAUSS = TRUE
params$USE_NM_GAUSS = TRUE
params$MAXIT_BFGS_GAUSS = 80
params$MAXIT_NM_GAUSS = 1000

params$IT_OPTIM_LAPLACE = 3
params$USE_BFGS_LAPLACE = TRUE
params$USE_NM_LAPLACE = TRUE
params$MAXIT_BFGS_LAPLACE = 100 # default for derivative based methods
params$MAXIT_NM_LAPLACE = 2000 # default for NM is 500

params$IT_OPTIM_SGT = 4
params$USE_BFGS_SGT = TRUE
params$USE_NM_SGT = TRUE
params$MAXIT_BFGS_SGT = 100 # default for derivative based methods
params$MAXIT_NM_SGT = 3000 # default for NM is 500

params$PATH_RESULTS_HELPER = "/proj/juhokois/sim_news/local_data/"

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

for(i in 1:nrep_est){
  
  # GENERATE DATA
  sim_obj <- do.call(what = sim_news, 
                     args = sim_prm)
  DATASET = apply(X = sim_obj$y$y,
                  MARGIN = 2,
                  FUN = function(x) x-mean(x))
  
  # Tibble with integer-valued parameters
  tt =
    # orders (p,q)
    expand_grid(p = 1, q = 1:2) %>% 
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
    expand_grid(sim_prm, data_list = list(DATASET)) %>% 
    mutate(sd_vec = map(.x = data_list, ~ apply(.x, 2, sd))) %>% 
    mutate(data_list = map(.x = data_list, ~ apply(.x, 2, function(x) x/sd(x))))
  
  # Parallel setup ####
  tt_optim_parallel = tt %>% 
    # template
    mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
    # generate initial values and likelihood functions (we can use the same template for initial values and likelihood fct bc both have no parameters for density)
    mutate(theta_init = map2(tmpl, data_list, ~get_init_armamod_whf_random(.y, .x))) %>% 
    mutate(shock_distr = "tdist") %>% 
    select(theta_init, tmpl, data_list, shock_distr)
  
  params_parallel = lapply(1:nrow(tt_optim_parallel),
                           function(i) t(tt_optim_parallel)[,i])
  
  mods_parallel_list = lapply(params_parallel, FUN = hlp_parallel)
  
  tibble_out =  
    enframe(mods_parallel_list) %>% 
    unnest_wider(value) %>% 
    unnest_wider(results_list) %>%
    unnest_wider(input_integerparams) %>% 
    mutate(tt) %>% 
    mutate(shock_distr = "tdist") %>% 
    mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>%
    select(p, q, kappa, k, n_st, n_unst, beta, rho, nu, value_final,
           params_deep_final,
           tmpl) %>% 
    mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~ irf_whf(.x, .y, n_lags = 8))) %>% 
    mutate(irf = map2(.x = sd_vec, .y = irf, ~ diag(.x)%r%.y)) %>% 
    mutate(rmat = map(.x = irf, ~ choose_perm_sign(cand_mat = unclass(.x)[,,1],
                                                   target_mat = sim_obj$mod$sigma_L,
                                                   type = "frob"))) %>%
    mutate(rmat = map2(.x = irf, .y = rmat, ~ .y%*%optim_zr(unclass(.x)[,,1]%*%.y, c(1,2), opt_it = FALSE))) %>%
    mutate(irf = map2(.x = irf, .y = rmat, ~ .x%r%.y))
  
  tibble_id <- paste0("/tibble_",
                      paste(sample(0:9, 5, replace = TRUE), collapse = ""), 
                      paste(sample(letters, 5), collapse = ""),
                      ".rds")
  
  saveRDS(tibble_out, file = paste0(new_dir_path, tibble_id))
}