irf_svar <- readRDS("~/Downloads/irf_svar.rds")
irf_svarma <- readRDS("~/Downloads/irf_svarma(4).rds")
list_svar <- lapply(1:dim(irf_svar)[5], function(x) apply(irf_svar[,,,,x], c(1,2,3), median)[,2,,drop=FALSE])
list_svarma <- lapply(1:dim(irf_svarma)[5], function(x) apply(irf_svarma[,,,,x], c(1,2,3), median)[,2,,drop=FALSE])
