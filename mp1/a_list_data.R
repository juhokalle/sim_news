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
         DLIP = c(NA, diff(LIP))) %>% 
  filter(date>ym(197212))
ebp <- read_csv("https://www.federalreserve.gov/econres/notes/feds-notes/ebp_csv.csv")
fred_md$EBP <- ebp$ebp[1:nrow(fred_md)]

fred_md <- list(fred_md, readRDS("local_data/shock_tbl.rds")) %>% 
  reduce(left_join, by = "date")

data_list <- list()

# GK1
data_list [[1]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, `S&P 500`, FEDFUNDS, ff4_tc) %>%
  mutate(ff4_tc = cumsum(coalesce(ff4_tc, 0)) + ff4_tc*0) %>% 
  filter(complete.cases(.), date >= ym(199401), date<ym(201401)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))

# GK2
data_list[[2]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, `S&P 500`, FEDFUNDS, mp1_tc) %>%
  mutate(mp1_tc = cumsum(coalesce(mp1_tc, 0)) + mp1_tc*0) %>% 
  filter(complete.cases(.), date >= ymd(19940101), date<ymd(20140101)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))

# RR04
data_list[[3]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, `S&P 500`, FEDFUNDS, resid_full) %>% 
  mutate(resid_full = cumsum(coalesce(resid_full, 0)) + resid_full*0) %>%
  filter(complete.cases(.), date >= ymd(19940101), date<ymd(20140101)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))

# BRW21
data_list[[4]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, `S&P 500`, FEDFUNDS, BRW_monthly) %>% 
  mutate(BRW_monthly = cumsum(coalesce(BRW_monthly, 0)) + BRW_monthly*0) %>%
  filter(complete.cases(.), date >= ymd(19940101), date<ymd(20140101)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))

# MAR21a
data_list[[5]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, `S&P 500`, FEDFUNDS, MM_IV1) %>%
  mutate(MM_IV1 = cumsum(coalesce(MM_IV1, 0)) + MM_IV1*0) %>%
  filter(complete.cases(.), date >= ymd(19940101), date<ymd(20140101)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))

# MAR21b
data_list[[6]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, `S&P 500`, FEDFUNDS, MM_IV5) %>% 
  mutate(MM_IV5 = cumsum(coalesce(MM_IV5, 0)) + MM_IV5*0) %>%
  filter(complete.cases(.), date >= ymd(19940101), date<ymd(20140101)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))

# Jaro22
data_list[[7]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, `S&P 500`, FEDFUNDS, u1) %>% 
  mutate(u1 = cumsum(coalesce(u1, 0)) + u1*0) %>% 
  filter(complete.cases(.), date >= ymd(19940101), date<ymd(20140101)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))

# AD22
data_list[[8]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, `S&P 500`, FEDFUNDS, Shock) %>%
  mutate(Shock = cumsum(coalesce(Shock, 0)) + Shock*0) %>%
  filter(complete.cases(.), date >= ymd(19940101), date<ymd(20140101)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))

# BS22a
data_list[[9]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, `S&P 500`, FEDFUNDS, MPS) %>%
  mutate(MPS = cumsum(coalesce(MPS, 0)) + MPS*0) %>%
  filter(complete.cases(.), date >= ymd(19940101), date<ymd(20140101)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))

# BS22b
data_list[[10]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, `S&P 500`, FEDFUNDS, MPS_ORTH) %>%
  mutate(MPS_ORTH = cumsum(coalesce(MPS_ORTH, 0)) + MPS_ORTH*0) %>%
  filter(complete.cases(.), date >= ymd(19940101), date<ymd(20140101)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))

# GSS22
data_list[[11]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, `S&P 500`, FEDFUNDS, MP1) %>%
  mutate(MP1 = cumsum(coalesce(MP1, 0)) + MP1*0) %>%
  filter(complete.cases(.), date >= ymd(19940101), date<ymd(20140101)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))


data_list <- tibble(data_list, type = c("GK1", "GK2", "RR04", "BRW21", "MAR21a", "MAR21b",
                                        "Jaro22", "AD22", "BS22a", "BS22b", "GSS22"))
saveRDS(data_list, "local_data/svarma_data_list.rds")
