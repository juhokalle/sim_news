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

wx_rate <- readxl::read_excel("local_data/WuXiaShadowRate.xlsx", 
                              sheet = "Data")[,3]
colnames(wx_rate) <- "wx_rate"
wx_rate$date <- seq(lubridate::ym(196001), by ="month", length.out = nrow(wx_rate))
# join shocks to fred table
fred_md <- list(fred_md, wx_rate, readRDS("local_data/shock_tbl.rds")) %>% 
  reduce(left_join, by = "date")
fred_md$FEDFUNDS_A = fred_md$FEDFUNDS
zlb <- seq(ym(200812), ym(201512), by = "month")
fred_md$FEDFUNDS_A[fred_md$date %in% zlb] <- wx_rate$wx_rate[wx_rate$date %in% zlb]

# mp shock identifiers
mp_id <- c("u1", "BRW_monthly", "MPS_ORTH", "MP1", "MP_median", "mp1_tc")

# names of papers shocks taken from
mp_type <- c("Jaro22", "BRW21", "BS22", "GSS22", "JK20", "GK15")

# choose baseline variables
baseline_data <- list(fred_md %>% dplyr::select(date, LIP, LCPI, FEDFUNDS) %>% 
                        filter(date >= ym(199401),
                               date<=ym(201206)),
                      fred_md %>% dplyr::select(date, LIP, PI, FEDFUNDS) %>% 
                        filter(date >= ym(199401),
                               date<=ym(201206)),
                      fred_md %>% dplyr::select(date, DLIP, DLCPI, FEDFUNDS) %>% 
                        filter(date >= ym(199401),
                               date<=ym(201206))
                      )
dl <- replicate(length(mp_id), baseline_data, simplify = FALSE) %>% 
  unlist(recursive = FALSE)

# sample span, and linear detrending
qq <- outer(mp_id, letters[1:3], paste, sep ="_") %>% t %>% c
dl <- map(qq, ~ mutate(dl[[which(qq %in% .x)]] %>%
                         inner_join(fred_md %>% 
                                      dplyr::select(date, all_of(gsub('.{2}$', '', .x))) %>% 
                                      rename(MPR = all_of(gsub('.{2}$', '', .x))),
                                    by="date")  %>% 
                         mutate(MPR = cumsum(coalesce(MPR, 0)) + MPR*0) %>% 
                         filter(complete.cases(.)) %>% 
                         dplyr::select(-date)))
                         # mutate(across(MPR, ~ lm(.x ~ I(1:n())) %>% residuals))))

# standardise data and save sd's for later analysis
data_list <- tibble(data_list = lapply(dl, function(x) x %>% mutate_all(~(.x - mean(.x))/sd(.x))),
                    std_dev = lapply(dl, function(x) apply(x, 2, sd)),
                    mp_type = rep(mp_type, each = 3),
                    log_level = rep(c("all", "pi", "none"), length(mp_id)))
# save data
saveRDS(data_list, "local_data/svarma_data_list.rds")
