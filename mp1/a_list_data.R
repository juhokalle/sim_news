library(tidyverse)
library(lubridate)
library(fbi)
fred_md <- fredmd("http://files.stlouisfed.org/files/htdocs/fred-md/monthly/current.csv",
                  transform = FALSE,
                  date_start = ymd(19720101)) %>%
  as_tibble() %>% 
  mutate(LIP = 100*log(INDPRO),
         LCPI = 100*log(CPIAUCSL),
         PI = c(rep(NA, 12), diff(LCPI, 12)),
         DLCPI = c(NA, diff(LCPI)),
         DLIP = c(NA, diff(LIP))) %>% 
  filter(date>ymd(19721201))
ebp <- read_csv("https://www.federalreserve.gov/econres/notes/feds-notes/ebp_csv.csv")
fred_md$EBP <- ebp$ebp[1:nrow(fred_md)]

fred_md <- list(fred_md, readRDS("local_data/shock_tbl.rds")) %>% 
  reduce(left_join, by = "date")

data_list <- list()
# GK1
data_list [[1]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, FEDFUNDS, ff4_tc) %>%
  mutate(ff4_tc = cumsum(coalesce(ff4_tc, 0)) + ff4_tc*0) %>% 
  filter(complete.cases(.)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))
  # mutate_at(vars(LIP, LCPI, ff4_tc), ~ lm(.x ~ I(1:n())) %>% residuals) %>% 
  # mutate_at(vars(EBP, FEDFUNDS), ~ .x - mean(.x))
# GK2
data_list[[2]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, FEDFUNDS, mp1_tc) %>%
  mutate(mp1_tc = cumsum(coalesce(mp1_tc, 0)) + mp1_tc*0) %>% 
  filter(complete.cases(.)) %>% 
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_all(~ .x/sd(.x))
  # mutate_at(vars(LIP, LCPI, mp1_tc), ~ lm(.x ~ I(1:n())) %>% residuals) %>% 
  # mutate_at(vars(EBP, FEDFUNDS), ~ .x - mean(.x))

data_list <- tibble(data_list, type = c("GK1", "GK2"))
saveRDS(data_list, "local_data/svarma_data_list.rds")
