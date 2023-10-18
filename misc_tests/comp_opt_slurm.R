# ----------------------------- #
# Script to be called via SLURM #
# ----------------------------- #

# PREAMBLE ####
.libPaths(c("/home/juhokois/proj/R/", .libPaths()))
pkgs <- c("svarmawhf", "tidyverse")
void = lapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE)))

# Arguments from Rscript call: Parameters from SLURM script ####
args = commandArgs(trailingOnly=TRUE)

tbl_out <- expand_grid(nobs = c(250, 500),
                       dim_out = c(2,3,5),
                       ar_ord = 2,
                       ma_ord = c(2,3),
                       k_val = 1,
                       kappa_val = 1,
                       max_p_q = 100,
                       init_opt = c("min", "max"),
                       path2f = "/home/juhokois/proj/sim_news/list_of_functions.R")

params_parallel <- split(tbl_out, seq(nrow(tbl_out)))
tbl_out$res <- lapply(params_parallel, function(x) do.call(comp_fn, x))

tibble_id <- paste0("/tibble_",
                    paste(sample(letters, 5, replace = TRUE), collapse = ""),
                    paste(sample(0:9, 5, replace = TRUE), collapse = ""), 
                    ".rds")

saveRDS(tbl_out, file = paste0(args[1], tibble_id))