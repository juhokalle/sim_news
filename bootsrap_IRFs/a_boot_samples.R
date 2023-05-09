# Create bootstrap samples for the IRF estimation
# Packages ####
source("list_of_functions.R")
pkgs <- c("svarmawhf", "tidyverse")
void = lapply(pkgs, library, character.only = TRUE)
nboot <- 2000

tbl0 <- readRDS("./local_data/target_model.rds")
ds <- readRDS("./local_data/jobid_20230303/total_data_sim.rds") %>% 
  slice(tbl0$nr) %>%  
  pull(data_list) %>% 
  .[[1]]
  
bl_vec <- c(5, 10, 20, 50)

arg_list <- map(bl_vec, ~
                  list(y = ds, 
                       prms = tbl0 %>% pull(params_deep_final) %>% .[[1]],
                       tmpl = tbl0 %>% pull(tmpl) %>% .[[1]],
                       b.length = .x,
                       nboot = nboot)
                )

dl <- map(arg_list, ~ do.call(mb_boot, .x)) %>% unlist(recursive = FALSE)

# standardize data for estimation
data_list <- tibble(data_list = lapply(dl, function(x) apply(x, 2, function(xx) (xx-mean(xx))/sd(xx))),
                    std_dev = lapply(dl, function(x) apply(x, 2, sd)),
                    mb_length = rep(bl_vec, each = nboot))

saveRDS(data_list, file = "./local_data/data_list_boot.rds")
