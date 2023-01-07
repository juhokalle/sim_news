# --------------------------------------------- #
# This script downloads the data for estimation #
# --------------------------------------------- #

library(tidyverse)
library(lubridate)
library(fbi)

# download fred data set
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

# download ebp
ebp <- read_csv("https://www.federalreserve.gov/econres/notes/feds-notes/ebp_csv.csv")
fred_md$EBP <- ebp$ebp[1:nrow(fred_md)]

# join shocks to fred table
fred_md <- list(fred_md, readRDS("local_data/shock_tbl.rds")) %>% 
  reduce(left_join, by = "date")

# mp shock identifiers
mp_id <- c("u1", "BRW_monthly", "MPS_ORTH", "MP1", "MP_median")

# names of papers shocks taken from
mp_type <- c("Jaro22", "BRW21", "BS22", "GSS22", "JK20")

# choose baseline variables
dl <- length(mp_id) %>% 
  replicate(fred_md %>% dplyr::select(date, LIP, LCPI, FEDFUNDS) %>% list)

# sample span, and linear detrending
dl <- mp_id %>% 
  lapply(function(x) mutate(dl[[which(mp_id %in% x)]] %>%
                              mutate(fred_md %>% dplyr::select(all_of(x)) %>% rename(MPR = all_of(x))) %>% 
                              # mutate(MPR = cumsum(coalesce(MPR, 0)) + MPR*0) %>% 
                              filter(complete.cases(.), date >= ym(199401), date<=ym(201512)) %>%
                              dplyr::select(-date) %>% 
                              mutate(across(!MPR, ~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals))))

# standardise data and save sd's for later analysis
data_list <- tibble(data_list = lapply(dl, function(x) x %>% mutate_all(~(.x - mean(.x))/sd(.x))),
                    std_dev = lapply(dl, function(x) apply(x, 2, sd)),
                    mp_type)
# save data
saveRDS(data_list, "local_data/svarma_data_list.rds")
