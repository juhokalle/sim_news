# ----------------------------- #
# Script to be called via SLURM #
# ----------------------------- #

# PREAMBLE ####
source("/home/juhokois/proj/sim_news/list_of_functions.R")
.libPaths(c("/home/juhokois/proj/R/", .libPaths()))
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
params$RESTART_W_NOISE = 3
params$FIX_INIT = FALSE
params$IC <- TRUE
params$penalty_prm = 25
params$AR_ORDER_MAX = 12
params$MA_ORDER_MAX = 3

## gaussian density
params$IT_OPTIM_GAUSS = 1
params$USE_BFGS_GAUSS = FALSE
params$USE_NM_GAUSS = FALSE
params$USE_CS_GAUSS = FALSE
params$MAXIT_BFGS_GAUSS = 200
params$MAXIT_NM_GAUSS = 3000
params$MAXIT_CS_GAUSS = 500

## laplacian density
params$IT_OPTIM_LAPLACE = 1
params$USE_BFGS_LAPLACE = FALSE
params$USE_NM_LAPLACE = FALSE
params$USE_CS_LAPLACE = TRUE
params$MAXIT_BFGS_LAPLACE = 200 # default for derivative based methods
params$MAXIT_NM_LAPLACE = 3000 # default for NM is 500
params$MAXIT_CS_LAPLACE = 1e5

## sgt density
params$IT_OPTIM_SGT = 3
params$USE_BFGS_SGT = TRUE
params$USE_NM_SGT = FALSE
params$USE_CS_SGT = FALSE
params$MAXIT_BFGS_SGT = 1e5 # default for derivative based methods
params$MAXIT_NM_SGT = 3000 # default for NM is 500
params$MAXIT_CS_SGT = 500
params$FTOL_REL = 1e-9

params$FILE_NAME_INPUT = "/home/juhokois/proj/sim_news/local_data/svarma_data_list.rds"
params$USE_PARALLEL = FALSE

# Data and derived PARAMETERS ####
DATASET = readRDS(params$FILE_NAME_INPUT)
DIM_OUT = dim(DATASET[[1]])[2]
params$DIM_OUT = DIM_OUT

# Tibble with integer-valued parameters
set.seed(123)
tt = 
  # orders (p,q)
  expand_grid(p = 1:params$AR_ORDER_MAX,
              q = 1:params$MA_ORDER_MAX) %>% 
  # number of unstable zeros
  mutate(n_unst = map(q, ~0:(DIM_OUT*.x))) %>% 
  unnest(n_unst) %>% 
  mutate(n_st = DIM_OUT * q - n_unst) %>% 
  mutate(kappa = n_unst %/% DIM_OUT,
         k = n_unst %% DIM_OUT) %>% 
  # filter(n_unst<=params$MA_ORDER_MAX*DIM_OUT/2) %>% 
  # Estimate SVAR(12) for comparison 
  # Join data sets and use two distributions for the estimation
  expand_grid(data_list = DATASET) %>% 
  # Select a subset of models to be estimated according to slurm task 
  slice(split(x = 1:n(),
              f = cut(x = sample(1:n(), replace = FALSE),
                      breaks = params$SLURM_ARRAY_TASK_MAX))[[params$IX_ARRAY_JOB]]
        )
rm(.Random.seed, envir=globalenv())

tt <- tt %>%
  # Standardize data and save sd's for IRFs
  mutate(sd_vec = map(.x = data_list, ~ apply(.x, 2, sd))) %>%
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  mutate(data_list = map(.x = data_list, 
                         ~ apply(.x, 2, function(x) (x - mean(x))/sd(x))
                         )
         ) %>%
  # generate initial values and likelihood functions (we can use the same template for initial values and likelihood fct bc both have no parameters for density)
  mutate(theta_init = map2(.x = data_list, 
                           .y = tmpl, 
                           ~ get_init_armamod_whf_random(.x, .y)
                           )
         ) %>%
  # estimate with a set of initial values
  mutate(theta_init = map2(.x = theta_init, 
                           .y = tmpl, 
                           ~ perm_init(.x, 500, .y, max_dist = TRUE)
                           )
         ) %>%
  # unnest_longer(theta_init) %>% 
  # mutate(init_ix = rep(1:(params$PERM_INIT+1), n()/(params$PERM_INIT+1))) %>% 
  # update template
  expand_grid(shock_distr = c("skewed_t", "sgt")) %>% 
  # template
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev))
  # filter(!(q==0 & init_ix>1))

# Parallel setup ####
tt_optim_parallel = tt %>% 
  dplyr::select(theta_init, tmpl, data_list, shock_distr)

params_parallel = lapply(1:nrow(tt_optim_parallel),
                         function(i) t(tt_optim_parallel)[,i])

if(params$USE_PARALLEL){
  
  # Parallel computations
  cl = try(makeCluster(params$N_CORES, type = "FORK"))
  if(inherits(cl, 'try-error')){
    mods_parallel_list = lapply(params_parallel, FUN = hlp_parallel)
  } else{
    mods_parallel_list <- clusterApply(cl, params_parallel, fun = hlp_parallel)
    stopCluster(cl)
  }
  
  cat("Parallel finished \n")
} else {
  mods_parallel_list = lapply(params_parallel, FUN = hlp_parallel)
}

tibble_out =  
  enframe(mods_parallel_list) %>% 
  unnest_wider(value) %>% 
  unnest_wider(results_list) %>%
  unnest_wider(input_integerparams) %>% 
  mutate(tt)
  # group_by(p, q, kappa, k, shock_distr, data_list) %>%
  # slice_min(value_final) %>%
  # ungroup()

tibble_id <- paste0("/tibble_",
                    paste(sample(0:9, 5, replace = TRUE), collapse = ""), 
                    paste(sample(letters, 5, replace = TRUE), collapse = ""),
                    ".rds")

saveRDS(tibble_out, file = paste0(params$NEW_DIR, tibble_id))
