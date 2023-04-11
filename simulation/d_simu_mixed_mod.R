# ----------------------------- #
# Script to be called via SLURM #
# ----------------------------- #

# PREAMBLE ####
source("/proj/juhokois/sim_news/list_of_functions.R")
.libPaths(c("/proj/juhokois/R/", .libPaths()))
pkgs <- c("svarmawhf", "fitdistrplus", "sgt", "tidyverse")
void = lapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE)))
set.seed(1303)
nrep_est <- 10

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
params$penalty_prm = 100

params$IT_OPTIM_GAUSS = 3
params$USE_BFGS_GAUSS = TRUE
params$USE_NM_GAUSS = TRUE
params$MAXIT_BFGS_GAUSS = 500
params$MAXIT_NM_GAUSS = 3000

params$IT_OPTIM_LAPLACE = 3
params$USE_BFGS_LAPLACE = TRUE
params$USE_NM_LAPLACE = TRUE
params$MAXIT_BFGS_LAPLACE = 500
params$MAXIT_NM_LAPLACE = 3000

params$IT_OPTIM_SGT = 4
params$USE_BFGS_SGT = TRUE
params$USE_NM_SGT = TRUE
params$MAXIT_BFGS_SGT = 1000
params$MAXIT_NM_SGT = 3000

params$PATH_RESULTS_HELPER = "/proj/juhokois/sim_news/local_data/"
params$DIM_OUT <- 4
DIM_OUT <- params$DIM_OUT
tmpl0 <- tmpl_whf_rev(dim_out = DIM_OUT,
                      ARorder = 1,
                      MAorder = 2,
                      kappa = 0,
                      k = 1,
                      shock_distr = "gaussian")
ar_flag <- imp_flag <- TRUE
while(ar_flag || imp_flag)
{
  imp_mat <- matrix(rnorm(DIM_OUT^2, 0, sqrt(2)), DIM_OUT, DIM_OUT)
  imp_mat <- imp_mat*sign(imp_mat)
  sign_mat <- matrix(1, DIM_OUT, DIM_OUT)
  sign_mat[upper.tri(sign_mat)] <- -1
  imp_mat <- imp_mat*sign_mat
  dgp_mod <- tmpl_whf_rev(DIM_OUT, 1, 2, 0, 1, "gaussian") %>%  
    get_simu_model(imp_mat)
  ar_flag <- min(abs(zeroes(dgp_mod$sys$a)))>1.2
  imp_flag <- max(abs(imp_mat))/min(abs(imp_mat))>4
}

# SIMULATION SPECS: TARGET PATH FOR SAVING RESULTS
pap = pap_factory(params$PATH_RESULTS_HELPER)

# New directory for saving all tibbles
new_dir_path = pap(paste0("jobid_", params$MANUALLY_ASSIGNED_ID))

if (!dir.exists(new_dir_path)){
  dir.create(new_dir_path)
}

tibble_out <- tibble()
for(i in 1:nrep_est){
  
  # GENERATE DATA
  DATASET = apply(X = do.call(what = simu_y, 
                              args = list(model = dgp_mod, 
                                          n.obs = 250,
                                          rand.gen = mixed_marg_dists(dim_out, 3),
                                          n.burnin = 500))$y,
                  MARGIN = 2,
                  FUN = function(x) x-mean(x))
  
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
    # this way of including data makes it convenient for slicing
    expand_grid(data_list = list(DATASET))
  
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
    mutate(res = pmap(., pmap_get_residuals_once)) %>% 
    select(p, q, kappa, k, n_st, n_unst, beta, rho, nu, value_final,
           params_deep_final,
           tmpl) %>% 
    mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~ irf_whf(.x, .y, n_lags = 8))) %>% 
    mutate(rmat = map(.x = irf, ~ choose_perm_sign(target_mat = sign_mat, 
                                                   cand_mat = unclass(.x)[,,1], 
                                                   type = "frob"))) %>%
    mutate(irf = map2(.x = irf, .y = rmat, ~ .x%r%.y)) %>% 
    mutate(B_mat = map2(.x = params_deep_final,.y = tmpl,
                        ~fill_tmpl_whf_rev(theta = .x,
                                           tmpl = .y)$B)) %>%
    mutate(B_mat = map2(.x = B_mat, .y = rmat, ~ .x%*%.y)) %>% 
    select(p, q, kappa, k, n_st, n_unst, beta, rho, nu, value_final, irf)
  
  tibble_id <- paste0("/tibble_",
                      paste(sample(0:9, 5, replace = TRUE), collapse = ""), 
                      paste(sample(letters, 5), collapse = ""),
                      ".rds")
  
  saveRDS(tibble_out, file = paste0(new_dir_path, tibble_id))
}