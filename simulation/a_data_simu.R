# --------------------------------------------------------------- #
# purpose and the output of the script: ------------------------- #
# 1) simulates data from the model with news shocks ------------- #
# with different degrees of non-Gaussianity and non-invertibility #
# 2) saving the dataset to be used in the estimation ------------ #
# --------------------------------------------------------------- #

# Packages ####
source("list_of_functions.R")
pkgs <- c("svarmawhf", "tidyverse")
void = lapply(pkgs, library, character.only = TRUE)

# Simulation params ####
mc_n <- 5000
sim_prm <- expand_grid(beta=c(0.5, 0.9), rho = 0.5, nobs = 250, nu = c(3, 20))
dl <- vector("list", mc_n*nrow(sim_prm))
# Simulation and data save ####
for(prm_ix in 1:nrow(sim_prm)){
  for(mc_ix in 1:mc_n){
    dl[[(prm_ix-1)*mc_n+mc_ix]] <- do.call(sim_news, sim_prm[prm_ix,])$y
  }
}
data_tbl <- tibble(data_list = map(dl, ~ apply(.x$y, 2, function(x) x-mean(x))), 
                   expand_grid(sim_prm, mc_ix = 1:mc_n))
#data_tbl$shock_list <- map(data_list, ~ .x$u)
saveRDS(data_tbl, file = "./local_data/data_list.rds")
