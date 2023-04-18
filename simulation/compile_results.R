# Collect results
library(tidyverse)
local_data <- "/proj/juhokois/sim_news/local_data/"
tibble_path <- local_data %>% 
  paste0(list.files(local_data)[grepl("jobid_", list.files(local_data))], "/")
tibble_list <- vector(mode = "list", length = length(list.files(tibble_path)))
for(i in seq_along(tibble_list)){ 
  file_path <- tibble_path %>% 
    paste0(paste0(tibble_path) %>% list.files() %>% .[i])
  tibble_list[[i]] <- file_path %>% readRDS()
  file_path %>% file.remove()
}
tibble_out <- do.call(bind_rows, tibble_list)
saveRDS(tibble_out, paste0(local_data, "tt_full.rds"))
