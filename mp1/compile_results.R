# Collect results
library(tidyverse)
local_data <- "/home/juhokois/proj/sim_news/local_data/"
tibble_path <- local_data %>% 
  paste0(list.files(local_data)[grepl("jobid_", list.files(local_data))], "/")
tibble_out <- paste0(tibble_path, list.files(tibble_path)) %>%
  map_dfr(readRDS) %>% 
  group_by(p, q, kappa, k, shock_distr, data_list) %>%
  slice_min(value_final) %>%
  ungroup()
saveRDS(tibble_out, paste0(local_data, "tt_simu.rds"))