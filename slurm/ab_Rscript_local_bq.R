
params = list()

# Parameters from Rmarkdown
params$USE_PARALLEL = TRUE
params$N_CORES = parallel::detectCores() - 1
params$MANUALLY_ASSIGNED_ID = "2021_11_15"

params$FILE_NAME_INPUT = "local_data/g_gmr_bq/data_xts.rds"

params$AR_ORDER_MAX = 5
params$MA_ORDER_MAX = 5

params$IT_OPTIM_GAUSS = 2
params$USE_BFGS_GAUSS = TRUE
params$USE_NM_GAUSS = TRUE
params$MAXIT_BFGS_GAUSS = 80
params$MAXIT_NM_GAUSS = 1000

params$IT_OPTIM_LAPLACE = 2
params$USE_BFGS_LAPLACE = TRUE
params$USE_NM_LAPLACE = TRUE
params$MAXIT_BFGS_LAPLACE = 100 # default for derivative based methods
params$MAXIT_NM_LAPLACE = 1000 # default for NM is 500

params$IT_OPTIM_SGT = 3
params$USE_BFGS_SGT = TRUE
params$USE_NM_SGT = TRUE
params$MAXIT_BFGS_SGT = 100 # default for derivative based methods
params$MAXIT_NM_SGT = 1000 # default for NM is 500

params$PATH_RESULTS_HELPER = "local_data/p_whf/local_bq_"

# Packages ####
pkgs <- c("lubridate", "xts", "RLDM", "parallel", "svarmawhf", "fitdistrplus", "sgt", "tidyverse")
void = lapply(pkgs, library, character.only = TRUE)

# Data and derived PARAMETERS ####
DATASET = readRDS(params$FILE_NAME_INPUT)

DATASET = DATASET %>% as.matrix()
DIM_OUT = dim(DATASET)[2]
params$DIM_OUT = DIM_OUT
N_OBS = dim(DATASET)[1]
params$N_OBS = N_OBS


# Tibble with integer-valued parameters
tt = 
  # orders (p,q)
  expand_grid(p = 0:params$AR_ORDER_MAX,
              q = 0:params$MA_ORDER_MAX) %>% 
  filter(!(p == 0 & q == 0)) %>% 
  # number of unstable zeros
  mutate(n_unst = map(q, ~0:(DIM_OUT*.x))) %>% 
  unnest(n_unst) %>% 
  mutate(n_st = DIM_OUT * q - n_unst) %>% 
  mutate(kappa = n_unst %/% DIM_OUT,
         k = n_unst %% DIM_OUT)

# Select rows ####
pmap_tmpl_whf_rev = function(dim_out = DIM_OUT, p, q, kappa, k, shock_distr = "gaussian", ...){
  tmpl_whf_rev(dim_out = DIM_OUT, ARorder = p, MAorder = q, kappa = kappa, k = k, shock_distr = shock_distr)
}

tt = tt %>% 
  # template
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  # generate initial values and likelihood functions (we can use the same template for initial values and likelihood fct bc both have no parameters for density)
  mutate(theta_init = map(tmpl, ~get_init_armamod_whf_random(DATASET, .x))) 



# Parallel setup ####
tt_optim_parallel = tt %>% 
  select(theta_init, tmpl)

params_parallel = lapply(1:nrow(tt_optim_parallel),
                         function(i) t(tt_optim_parallel)[,i])

hlp_parallel = function(list_input){
  return(create_results_list(theta_init = list_input[[1]], 
                             tmpl = list_input[[2]],
                             params = params,
                             DATASET = DATASET))
}



if (params$USE_PARALLEL){
  
  # Parallel computations
  cl = makeCluster(params$N_CORES, type = "FORK")
  mods_parallel_list <- clusterApply(cl, params_parallel, fun = hlp_parallel)
  stopCluster(cl)

  cat("Parallel finished \n")
} else {
  mods_parallel_list = lapply(params_parallel, FUN = hlp_parallel)
}

pap = pap_factory(params$PATH_RESULTS_HELPER)

# New directory for saving all
new_dir_path = pap(paste0("jobid_", params$MANUALLY_ASSIGNED_ID))
  
if (!dir.exists(new_dir_path)){
  dir.create(new_dir_path)
}

saveRDS(mods_parallel_list, paste0(new_dir_path, "/local_results.rds"))

