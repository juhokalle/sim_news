# Collect results
library(tidyverse)
local_data <- "/home/juhokois/proj/sim_news/misc_tests/"
tibble_path <- local_data %>% 
  paste0(list.files(local_data)[grepl("jobid_", list.files(local_data))], "/")
tibble_out <- paste0(tibble_path, list.files(tibble_path)) %>%
  map_dfr(readRDS)
saveRDS(tibble_out, paste0(local_data, "opt_comp.rds"))


