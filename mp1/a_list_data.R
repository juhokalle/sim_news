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
         SP500 = 100*log(`S&P 500`)) %>% 
  filter(date>ym(197212))

# Additional variables ####
# EBP
ebp <- read_csv("https://www.federalreserve.gov/econres/notes/feds-notes/ebp_csv.csv")
# Shadow rates: WX
fred_md$EBP <- ebp$ebp[1:nrow(fred_md)]
WX <- readxl::read_excel("local_data/WuXiaShadowRate.xlsx", 
                              sheet = "Data")[,3]
colnames(WX) <- "WX"
WX$date <- seq(ym(196001), by ="month", length.out = nrow(WX))
# Krippner SSR
SSR <- readxl::read_excel("local_data/SSR_Estimates_202206.xlsx", 
                          sheet = "D. Monthly average SSR series", 
                          skip = 19)[,3]
colnames(SSR) <- "SSR"
SSR$date <- seq(ym(199501), by= "month", length.out = nrow(SSR))

# Shocks ####
fred_md <- list(fred_md, WX, SSR, readRDS("local_data/shock_tbl.rds")) %>% 
  reduce(left_join, by = "date")
fred_md$FEDFUNDS_A = fred_md$FEDFUNDS
zlb <- seq(ym(200812), ym(201512), by = "month")
fred_md$FEDFUNDS_A[fred_md$date %in% zlb] <- SSR$SSR[SSR$date %in% zlb]

# mp shock identifiers
# mp_id <- c("u1", "BRW_monthly", "MPS_ORTH", "MP1", "MP_median")
# mp_id <- c("t_fac", "ffr_fac")
# names of papers shocks taken from
# mp_type <- c("GSS22", "Swanson21")

# choose baseline variables: six variants according to the differencing and policy indic.
data_l <- list(fred_md %>% dplyr::select(date, LIP, LCPI, FEDFUNDS_A),
               fred_md %>% dplyr::select(date, LIP, PI, FEDFUNDS_A),
               fred_md %>% dplyr::select(date, DLIP, DLCPI, FEDFUNDS_A),
               fred_md %>% dplyr::select(date, LIP, LCPI, GS1),
               fred_md %>% dplyr::select(date, LIP, PI, GS1),
               fred_md %>% dplyr::select(date, DLIP, DLCPI, GS1))

# replication thrice: 3 different sample spans
baseline_data <- replicate(n = 3, data_l, simplify = FALSE) %>% unlist(recursive=FALSE)

baseline_data[1:length(data_l)] <- map(baseline_data[1:length(data_l)], 
                                       ~ .x %>% filter(date >= ym(199401),
                                                       date<=ym(201912))
                                       )
baseline_data[(length(data_l)+1):(2*length(data_l))] <- map(baseline_data[(length(data_l)+1):(2*length(data_l))], 
                                                            ~ .x %>% filter(date >= ym(199401), date<=ym(201206))
                                                            )
baseline_data[(2*length(data_l)+1):(3*length(data_l))] <- map(baseline_data[(2*length(data_l)+1):(3*length(data_l))], 
                                                            ~ .x %>% filter(date >= ym(198301), date<=ym(200812))
                                                            )

baseline_data[rep(c(T, T, F), 6)] <- map(baseline_data[rep(c(T, T, F), 6)], 
                                         ~ .x %>% mutate_if(is.numeric, ~ lm(.x ~ I(1:n())) %>% residuals())
                                        )

# dl <- replicate(length(mp_id), baseline_data, simplify = FALSE) %>% 
#   unlist(recursive = FALSE)

# sample span, and linear detrending
# qq <- outer(mp_id, letters[1:4], paste, sep ="_") %>% t %>% c
# dl <- map(qq, ~ mutate(dl[[which(qq %in% .x)]] %>%
#                          inner_join(fred_md %>% 
#                                       dplyr::select(date, all_of(gsub('.{2}$', '', .x))) %>% 
#                                       rename(MPR = all_of(gsub('.{2}$', '', .x))),
#                                     by="date")  %>% 
#                          mutate(MPR =
#                                   if(str_sub(.x,-1)%in%letters[3:4])
#                                     {
#                                     cumsum(coalesce(MPR, 0)) + MPR*0
#                                     } else MPR
#                                 ) %>% 
#                          filter(complete.cases(.)) %>% 
#                          dplyr::select(-date) %>% 
#                          mutate_all(~ lm(.x ~ I(1:n())) %>% residuals())
#                        )
#           )

# standardise data and save sd's for later analysis
dl <- baseline_data %>% map(~ .x %>% dplyr::select(-date))
data_list <- tibble(data_list = map(dl, ~ .x %>% mutate_all(~(.x - mean(.x))/sd(.x))),
                    std_dev = map(dl, ~ apply(.x, 2, sd)))
data_list <- data_list %>% 
  bind_cols(FFR = rep(rep(c(T,F), each = 3), 3),
            GS1 = rep(rep(c(F,T), each = 3), 3),
            long_sample = rep(c(T,F,F), each = 6), 
            short_sample = rep(c(F,T,F), each = 6), 
            fgr_sample = rep(c(F,F,T), each = 6), 
            log_level = rep(c(T,F,F), 6),
            pi = rep(c(F,T,F), 6),
            log_diff = rep(c(F,F,T), 6)
            )
# data_list <- data_list %>% 
#   bind_cols(expand_grid(mp_type, mpr_lvl = c(T,F), log_lvl = c(T,F)))

# save data
saveRDS(data_list, "local_data/svarma_data_list.rds")
