library(tidyverse)
library(lubridate)
library(fbi)
fred_md <- fredmd("http://files.stlouisfed.org/files/htdocs/fred-md/monthly/current.csv",
                  transform = FALSE,
                  date_start = ym(197201)) %>%
  as_tibble() %>% 
  mutate(LIP = 100*log(INDPRO),
         LCPI = 100*log(CPIAUCSL),
         PI = c(rep(NA, 12), diff(LCPI, 12)),
         DLCPI = c(NA, diff(LCPI)),
         DLIP = c(NA, diff(LIP)),
         `S&P 500` = 100*log(`S&P 500`)) %>% 
  filter(date>ym(197212))
ebp <- read_csv("https://www.federalreserve.gov/econres/notes/feds-notes/ebp_csv.csv")
fred_md$EBP <- ebp$ebp[1:nrow(fred_md)]

fred_md <- list(fred_md, readRDS("local_data/shock_tbl.rds")) %>% 
  reduce(left_join, by = "date")

mp_id <- c("ff4_tc", "mp1_tc", "resid_full", "BRW_monthly",
           "MM_IV1", "MM_IV5", "u1", "Shock", "MPS",
           "MPS_ORTH", "MP1")
mp_type <- c("GK15a", "GK15b", "RR04", "BRW21", "MAR21a", "MAR21b",
             "Jaro22", "AD22", "BS22a", "BS22b", "GSS22")
data_list <- length(mp_id) %>% 
  replicate(fred_md %>% dplyr::select(date, LIP, LCPI, EBP, FEDFUNDS) %>% list)

data_list <- mp_id %>% 
  lapply(function(x) mutate(data_list[[which(mp_id %in% x)]] %>%
                              mutate(fred_md %>% dplyr::select(all_of(x)) %>% rename(MPR = all_of(x))) %>% 
                              mutate(MPR = cumsum(coalesce(MPR, 0)) + MPR*0) %>% 
                              filter(complete.cases(.), date >= ym(199401), date<ym(201401)) %>% 
                              dplyr::select(-date) %>% 
                              mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals)))

sd_list <- tibble(std_dev = lapply(data_list, function(x) apply(x, 2, sd)), mp_type)
data_list <- tibble(data_list = lapply(data_list, function(x) x %>% mutate_all(~.x/sd(.x))),
                    mp_type)
saveRDS(sd_list, "local_data/svarma_sd_list.rds")
saveRDS(data_list, "local_data/svarma_data_list.rds")

