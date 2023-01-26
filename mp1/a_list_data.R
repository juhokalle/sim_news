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

# download ebp
ebp <- read_csv("https://www.federalreserve.gov/econres/notes/feds-notes/ebp_csv.csv")
fred_md$EBP <- ebp$ebp[1:nrow(fred_md)]

WX <- readxl::read_excel("local_data/WuXiaShadowRate.xlsx", 
                              sheet = "Data")[,3]
colnames(WX) <- "WX"
WX$date <- seq(lubridate::ym(196001), by ="month", length.out = nrow(WX))
# join shocks to fred table
fred_md <- list(fred_md, WX, readRDS("local_data/shock_tbl.rds")) %>% 
  reduce(left_join, by = "date")
fred_md$FEDFUNDS_A = fred_md$FEDFUNDS
zlb <- seq(ym(200812), ym(201512), by = "month")
fred_md$FEDFUNDS_A[fred_md$date %in% zlb] <- WX$WX[WX$date %in% zlb]

# mp shock identifiers
mp_id <- c("u1", "BRW_monthly", "MPS_ORTH", "MP1", "MP_median")

# names of papers shocks taken from
mp_type <- c("Jaro22", "BRW21", "BS22", "GSS22", "JK20")

# choose baseline variables
baseline_data <- replicate(n = 4,
                           list(fred_md %>% dplyr::select(date, LIP, LCPI, FEDFUNDS_A) %>% 
                                  filter(date >= ym(199401),
                                         date<=ym(201912)),
                                fred_md %>% dplyr::select(date, LIP, PI, FEDFUNDS_A) %>% 
                                  filter(date >= ym(199401),
                                         date<=ym(201912))
                                ),
                           simplify = FALSE) %>% unlist(recursive=FALSE)

baseline_data[c(1,2,5,6)] <- map(baseline_data[c(1,2,5,6)], ~ .x %>% filter(date<=ym(200812)))

dl <- replicate(length(mp_id), baseline_data, simplify = FALSE) %>% 
  unlist(recursive = FALSE)

# sample span, and linear detrending
qq <- outer(mp_id, letters[1:8], paste, sep ="_") %>% t %>% c
dl <- map(qq, ~ mutate(dl[[which(qq %in% .x)]] %>%
                         inner_join(fred_md %>% 
                                      dplyr::select(date, all_of(gsub('.{2}$', '', .x))) %>% 
                                      rename(MPR = all_of(gsub('.{2}$', '', .x))),
                                    by="date")  %>% 
                         mutate(MPR = 
                                  if(str_sub(.x,-1)%in%letters[5:8]){
                                    cumsum(coalesce(MPR, 0)) + MPR*0
                                    } else MPR) %>% 
                         filter(complete.cases(.)) %>% 
                         dplyr::select(-date))
          )

# standardise data and save sd's for later analysis
data_list <- tibble(data_list = lapply(dl, function(x) x %>% mutate_all(~(.x - mean(.x))/sd(.x))),
                    std_dev = lapply(dl, function(x) apply(x, 2, sd)),
                    smpl_s = rep(rep(c(TRUE, FALSE), each=2), 10),
                    mpr_lvl = rep(rep(c(TRUE, FALSE), each = 4), 5),
                    mp_type = rep(mp_type, each = 8),
                    log_level = rep(c(TRUE, FALSE), 20))
# save data
saveRDS(data_list, "local_data/svarma_data_list.rds")
