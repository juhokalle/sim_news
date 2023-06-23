# --------------------------------------------- #
# This script downloads the data for estimation #
# --------------------------------------------- #

library(tidyverse)
library(lubridate)
library(fbi)

# load fred data set
if(!file.exists("local_data/fred_md.rds")){
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
  fred_md$EBP <- ebp$ebp[1:nrow(fred_md)]
  saveRDS(fred_md, "local_data/fred_md.rds")
} else{
  fred_md <- readRDS("local_data/fred_md.rds")
}

# Shadow rates: WX
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

data_list <- map(c("BRW_monthly", "MPS_ORTH", "ffr_fac", "MP1"), 
                 ~ fred_md %>%
                   filter(date>=ym(199401), date<=ym(201912)) %>% 
                   dplyr::select(LIP, LCPI, FEDFUNDS, all_of(.x)) %>% 
                   filter(complete.cases(.))
                 )

# save data
saveRDS(data_list, "local_data/svarma_data_list.rds")
