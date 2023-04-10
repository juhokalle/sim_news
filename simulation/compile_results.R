# Collect results
library(tidyverse)
local_data <- "/proj/juhokois/sim_news/local_data/"
tibble_path <- local_data %>% 
  paste0(list.files(local_data)[grepl("jobid_", list.files(local_data))], "/")
tibble_list <- vector(mode = "list", legnth = length(list.files(tibble_path)))
for(i in seq_along(tibble_list)){ 
  tibble_list[[i]] <- paste0(tibble_path, 
                             paste0(tibble_path) %>% list.files() %>% .[i]) %>% 
    readRDS()
}
tibble_out <- reduce(tibble_list, bind_rows)
saveRDS(tibble_out, paste0(local_data, "tt_full.rds"))
