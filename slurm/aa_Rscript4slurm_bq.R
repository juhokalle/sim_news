# Script to be called via SLURM

# Arguments from Rscript call: Parameters from SLURM script ####
args = commandArgs(trailingOnly=TRUE)
params = list()
.libPaths(c(.libPaths(), "/proj/juhokois/R/"))

# Parameters from Rmarkdown
params$USE_PARALLEL = FALSE
if (params$USE_PARALLEL){
  params$N_CORES = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
} else {
  params$N_CORES = 1
} 
params$N_MODS_PER_CORE = as.integer(args[1]) # important param: specifies how many models are estimated by each array-job
params$IX_ARRAY_JOB = as.integer(args[2]) # index of array-job. Number of array-jobs is determined from number of rows of dataframe containing all integer-values parameters
params$SLURM_JOB_ID = as.integer(args[3])
params$MANUALLY_ASSIGNED_ID = as.integer(args[4])

params$FILE_NAME_INPUT = "./local_data/data_svarma.rds"

params$AR_ORDER_MAX = 5
params$MA_ORDER_MAX = 5

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

params$PATH_RESULTS_HELPER = "~/sim_news/local_data/"


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

# Packages ####
pkgs <- c("lubridate", "xts", "parallel", "svarmawhf", "fitdistrplus", "sgt", "tidyverse")
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
  slice((1 + (params$IX_ARRAY_JOB-1) * params$N_CORES * params$N_MODS_PER_CORE):(params$IX_ARRAY_JOB * params$N_CORES * params$N_MODS_PER_CORE)) %>% 
  # template
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  # generate initial values and likelihood functions (we can use the same template for initial values and likelihood fct bc both have no parameters for density)
  mutate(theta_init = map(tmpl, ~get_init_armamod_whf_random(DATASET, .x))) #%>% 
  #mutate(ll_fun_gaussian = map(tmpl, ~ll_whf_factory(data_wide = t(DATASET), tmpl = .x, shock_distr = "gaussian", use_cpp = TRUE)))

# Parallel setup ####
tt_optim_parallel = tt %>% 
  select(theta_init, tmpl)

params_parallel = lapply(1:nrow(tt_optim_parallel),
                         function(i) t(tt_optim_parallel)[,i])

hlp_parallel = function(list_input){
  return(create_results_list(theta_init = list_input[[1]], 
                             tmpl       = list_input[[2]],
                             params     = params, 
                             DATASET    = DATASET))
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

saveRDS(mods_parallel_list, paste0(new_dir_path, "/arrayjob_", params$IX_ARRAY_JOB,".rds"), version = 3)
