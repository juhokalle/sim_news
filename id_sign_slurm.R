# -------------------------------------------------------------- #
# This script generates the output regarding the main simulation #
# exercise. The simulation results are loaded from a data frame  #
# containing the necessary components for this script to run and #
# make the figures and tables in the paper. -------------------- #
# -------------------------------------------------------------- #

pkgs <- c("svarmawhf", "tidyverse")
void = lapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE)))
source("/home/juhokois/proj/sim_news/list_of_functions.R")

args = commandArgs(trailingOnly=TRUE)
path2results <- args[[1]]
## IRF sec. 7.2.
tt_est <- readRDS("/home/juhokois/proj/sim_news/local_data/aux_est_mod.rds")

# CHOL light restrictions
rest_hor <- 6
dim_out <- 5
sgn_mat <- array(NA, c(dim_out, dim_out, rest_hor))
sgn_mat[,1,] <- c(1, 1, NA, 1, NA) # Demand shock 
sgn_mat[,2,] <- c(-1, 1, NA, 1, NA) # Supply shock
sgn_mat[,3,] <- c(NA, NA, 1, NA, NA) # Financial shock
sgn_mat[,4,1] <- c(NA, -1, -1, 1, 1) # MP impact
sgn_mat[,4,-1] <- c(NA, -1, -1, 1, NA) # MP dynamic
sgn_mat[,5,1] <- c(NA, -1, -1,  NA, 0)  # FG impact
sgn_mat[,5,-1] <- c(NA, -1, -1, NA, NA)  # FG dynamic

sr_list <- list(sgn_mat)
sgn_mat[,3,] <- sgn_mat[,3,] <- c(1, 1, 1, 1, NA)
sr_list[[2]] <- sgn_mat

param_list <- expand_grid(w_mat = list(tt_est$B_mat[[1]], "chol"),
                          rest_mat = sr_list)
param_list$irf_arr <- list(unclass(tt_est$irf[[1]])[,,1:37])
param_list$emp_innov <- tt_est$res
param_list$news_rest <- list(c(1, 4, 5))
param_list$rest_mom <- list(list(type = "kurt", 
                                 test_shocks = 4:5, 
                                 threshold = .5))
param_list$ndraws <- 1
param_list$sub_max <- 100
param_list$max_draws <- 10
param_list$replace_md <- FALSE
param_list = lapply(1:nrow(param_list),
                    function(i) t(param_list)[,i])
est_obj <- lapply(param_list, function(x) do.call(id_mixed_new, x))
res_tbl <- expand_grid(bmat = c("ML", "chol"), sr = c("few", "many"))
res_tbl$results <- est_obj
tibble_id <- paste0(sample(0:9, 5, replace = TRUE),
                    sample(letters, 5, replace = TRUE))
tibble_id <- paste0("/tibble_",
                    paste(sample(0:9, 5, replace = TRUE), collapse = ""), 
                    paste(sample(letters, 5, replace = TRUE), collapse = ""),
                    ".rds")
saveRDS(res_tbl, paste0(path2results, tibble_id))