# Collect results
library(tidyverse)
local_data <- "/proj/juhokois/sim_news/local_data/"
tibble_path <- local_data %>% 
  paste0(list.files(local_data)[grepl("jobid_", list.files(local_data))], "/")
tibble_out <- paste0(tibble_path, list.files(tibble_path)) %>% map_dfr(readRDS)
saveRDS(tibble_out, paste0(local_data, "tt_boot.rds"))