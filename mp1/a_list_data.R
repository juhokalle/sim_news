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
# SHORT & STRONG PERSISTENCE
data_list[[1]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, FEDFUNDS) %>%
  filter(date >= ymd(19821001), date<ymd(20080101)) %>%
  mutate_at(vars(LIP, LCPI, FEDFUNDS), ~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_at(vars(EBP), ~ .x - mean(.x)) %>% 
  dplyr::select(-date)
# SHORT & MEDIUM PERSISTENCE
data_list[[2]] <- fred_md %>% dplyr::select(date, LIP, PI, EBP, FEDFUNDS) %>%
  filter(date >= ymd(19821001), date<ymd(20080101)) %>%
  mutate_at(vars(LIP, FEDFUNDS), ~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_at(vars(PI, EBP), ~ .x - mean(.x)) %>% 
  dplyr::select(-date)
# SHORT & LOW PERSISTENCE
data_list[[3]] <- fred_md %>% dplyr::select(date, DLIP, DLCPI, EBP, FEDFUNDS) %>%
  filter(date >= ymd(19821001), date<ymd(20080101)) %>%
  mutate_at(vars(DLIP, DLCPI, EBP), ~.x - mean(.x)) %>% 
  mutate_at(vars(FEDFUNDS), ~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>%
  dplyr::select(-date)
# LONG & STRONG PERSISTENCE
data_list[[4]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, FEDFUNDS) %>%
  filter(date<ymd(20080101)) %>%
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals)
# LONG & MEDIUM PERSISTENCE
data_list[[5]] <- fred_md %>% dplyr::select(date, LIP, PI, EBP, FEDFUNDS) %>%
  filter(date<ymd(20080101)) %>%
  dplyr::select(-date) %>% 
  mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals)
# LONG & LOW PERSISTENCE
data_list[[6]] <- fred_md %>% dplyr::select(date, DLIP, DLCPI, EBP, FEDFUNDS) %>%
  filter(date<ymd(20080101)) %>%
  mutate_at(vars(DLIP), ~.x - mean(.x)) %>%
  mutate_at(vars(DLCPI, EBP, FEDFUNDS), ~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  dplyr::select(-date)
# SHORTEST & STRONG PERSISTENCE
data_list [[7]] <- fred_md %>% dplyr::select(date, LIP, LCPI, EBP, FEDFUNDS) %>%
  filter(date >= ymd(19900101), date<ymd(20080101)) %>%
  dplyr::select(-date) %>% 
  mutate_at(vars(LIP, LCPI, FEDFUNDS),~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_at(vars(EBP), ~ .x - mean(.x))
# SHORTEST & MEDIUM PERSISTENCE
data_list [[8]] <- fred_md %>% dplyr::select(date, LIP, PI, EBP, FEDFUNDS) %>%
  filter(date >= ymd(19900101), date<ymd(20080101)) %>%
  dplyr::select(-date) %>% 
  mutate_at(vars(LIP, PI, FEDFUNDS), ~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_at(vars(EBP), ~ .x - mean(.x))
# SHORTEST & LOW PERSISTENCE
data_list [[9]] <- fred_md %>% dplyr::select(date, DLIP, DLCPI, EBP, FEDFUNDS) %>%
  filter(date >= ymd(19900101), date<ymd(20080101)) %>%
  dplyr::select(-date) %>% 
  mutate_at(vars(DLCPI, FEDFUNDS), ~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) %>% 
  mutate_at(vars(DLIP, EBP), ~ .x - mean(.x))

data_list <- tibble(data_list,
                    expand_grid(length= c("med", "long", "short"),
                                prst = c("str", "med", "low")))
saveRDS(data_list, "local_data/svarma_data_list.rds")