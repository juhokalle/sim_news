# ----------------------------- #
# Script to be called via SLURM #
# ----------------------------- #

# PREAMBLE ####
source("/proj/juhokois/sim_news/list_of_functions.R")
.libPaths(c("/proj/juhokois/R/", .libPaths()))
pkgs <- c("lubridate", "xts", "parallel", "svarmawhf", 
          "mixtools", "fitdistrplus", "sgt", "tidyverse")
void = lapply(pkgs, library, character.only = TRUE)

# Arguments from Rscript call: Parameters from SLURM script ####
args = commandArgs(trailingOnly=TRUE)
params = list()

# Parameters from Rmarkdown
params$USE_PARALLEL = TRUE
if (params$USE_PARALLEL){
  params$N_CORES = as.integer(args[6])
} else {
  params$N_CORES = 1
}

params$RESTART_W_NOISE = 0
params$FIX_INIT = FALSE

params$N_MODS_PER_CORE = as.integer(args[1]) # important param: specifies how many models are estimated by each array-job
params$IX_ARRAY_JOB = as.integer(args[2]) # index of array-job. Number of array-jobs is determined from number of rows of dataframe containing all integer-values parameters
params$SLURM_JOB_ID = as.integer(args[3])
params$MANUALLY_ASSIGNED_ID = as.integer(args[4])

params$FILE_NAME_INPUT = "/proj/juhokois/sim_news/local_data/svarma_data_list.rds"

params$AR_ORDER_MAX = 4
params$MA_ORDER_MAX = 4

params$IT_OPTIM_GAUSS = 3
params$USE_BFGS_GAUSS = TRUE
params$USE_NM_GAUSS = TRUE
params$MAXIT_BFGS_GAUSS = 500
params$MAXIT_NM_GAUSS = 1000

params$IT_OPTIM_LAPLACE = 3
params$USE_BFGS_LAPLACE = TRUE
params$USE_NM_LAPLACE = TRUE
params$MAXIT_BFGS_LAPLACE = 500 # default for derivative based methods
params$MAXIT_NM_LAPLACE = 1000 # default for NM is 500

params$IT_OPTIM_SGT = 4
params$USE_BFGS_SGT = TRUE
params$USE_NM_SGT = TRUE
params$MAXIT_BFGS_SGT = 1000 # default for derivative based methods
params$MAXIT_NM_SGT = 2500 # default for NM is 500

params$PATH_RESULTS_HELPER = "/proj/juhokois/sim_news/local_data/"

cat("\n--------------------------------------------------\n")
cat(paste0("This is array task ", params$IX_ARRAY_JOB, "\n"))
cat(paste0("The job ID is ", params$SLURM_JOB_ID, "\n"))
cat(paste0("The manually assigned ID is ", params$MANUALLY_ASSIGNED_ID, "\n\n"))
cat(paste0("The number of cores of this job is ", params$N_CORES, 
           "and there are ", params$N_MODS_PER_CORE,
           " models estimated per core, for a total of ", params$N_CORES*params$N_MODS_PER_CORE, "\n\n"))
cat(paste0("Rows ", 1 + (params$IX_ARRAY_JOB-1) * params$N_CORES * params$N_MODS_PER_CORE,
           " to ", params$IX_ARRAY_JOB * params$N_CORES * params$N_MODS_PER_CORE,
           " from the tibble containing integer-valued parameters are chosen. \n\n"))
cat("\n--------------------------------------------------\n")


# Data and derived PARAMETERS ####
DATASET = readRDS(params$FILE_NAME_INPUT)

DIM_OUT = dim(DATASET$data_list[[1]])[2]
params$DIM_OUT = DIM_OUT

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
         k = n_unst %% DIM_OUT) %>% 
  # Estimate SVAR for comparison 
  bind_rows(tibble(p = (params$AR_ORDER_MAX+1):12, q = 0, n_unst = 0, n_st = 0, kappa = 0, k = 0)) %>% 
  expand_grid(DATASET, sd = c("tdist", "gt", "skewed_ged", "sgt", "skewed_t"))

pap = pap_factory(params$PATH_RESULTS_HELPER)

# New directory for saving all
new_dir_path = pap(paste0("jobid_", params$MANUALLY_ASSIGNED_ID))

if (!dir.exists(new_dir_path)){
  dir.create(new_dir_path)
}

if(params$IX_ARRAY_JOB==1){
  saveRDS(tt, file = paste0(new_dir_path, "/total_data.rds"))
}

# Parallel setup ####
tt_optim_parallel = tt %>% 
  slice((1 + (params$IX_ARRAY_JOB-1) * params$N_CORES * params$N_MODS_PER_CORE):(params$IX_ARRAY_JOB * params$N_CORES * params$N_MODS_PER_CORE)) %>% 
  # template
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>%
  # generate initial values and likelihood functions (we can use the same template for initial values and likelihood fct bc both have no parameters for density)
  mutate(theta_init = map2(tmpl, data_list, ~get_init_armamod_whf_random(.y, .x))) %>% 
  select(theta_init, tmpl, data_list, sd) %>% 
  rename(shock_distr = sd)

params_parallel = lapply(1:nrow(tt_optim_parallel),
                         function(i) t(tt_optim_parallel)[,i])

if(params$USE_PARALLEL){
  
  # Parallel computations
  cl = makeCluster(params$N_CORES, type = "FORK")
  mods_parallel_list <- clusterApply(cl, params_parallel, fun = hlp_parallel)
  stopCluster(cl)
  
  cat("Parallel finished \n")
} else {
  mods_parallel_list = lapply(params_parallel, FUN = hlp_parallel)
}

saveRDS(mods_parallel_list, paste0(new_dir_path, "/arrayjob_", 
                                   if(params$IX_ARRAY_JOB<10) "00" else if(params$IX_ARRAY_JOB<100) "0",
                                   params$IX_ARRAY_JOB,".rds"), 
        version = 3)

