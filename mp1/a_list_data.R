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

data_list <- list()
# SHORT & QUAD DETREND
data_list [[1]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, FEDFUNDS) %>%
  filter(date >= ymd(19940101), date<ymd(20190101)) %>%
  dplyr::select(-date) %>% 
  mutate_at(vars(LIP, LCPI), ~ lm(.x ~ I(1:n())) %>% residuals) %>% 
  mutate_at(vars(EBP, FEDFUNDS), ~ .x - mean(.x))
# LONG & QUAD DETREND
data_list [[2]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, FEDFUNDS) %>%
  filter(date >= ymd(19940101), date<ymd(20140101)) %>%
  dplyr::select(-date) %>% 
  mutate_at(vars(LIP, LCPI), ~ lm(.x ~ I(1:n())) %>% residuals) %>% 
  mutate_at(vars(EBP, FEDFUNDS), ~ .x - mean(.x))

data_list <- tibble(data_list, length = c("long", "short"))
saveRDS(data_list, "local_data/svarma_data_list.rds")
