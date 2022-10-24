mod0 <- sim_news(0.5, 0.5, 250, 12)$mod
irf_svar <- readRDS("~/Downloads/irf_svar.rds")
irf_svarma <- readRDS("~/Downloads/irf_svarma(4).rds")
list_svar <- lapply(1:dim(irf_svar)[5], function(x) apply(irf_svar[,,,,x], c(1,2,3), median)[,2,,drop=FALSE])
list_svarma <- lapply(1:dim(irf_svarma)[5], function(x) apply(irf_svarma[,,,,x], c(1,2,3), median)[,2,,drop=FALSE])
th.5 <- unclass(impresp(sim_news(0.5, 0.5, 250, 12)$mod, lag.max = 12, H = mod0$sigma_L)$irf)[,2,,drop=FALSE]
th.9 <- unclass(impresp(sim_news(0.9, 0.5, 250, 12)$mod, lag.max = 12, H = mod0$sigma_L)$irf)[,2,,drop=FALSE]

irf_tbl <- tibble(
  svar11 = list_svar[[1]][1,,],
  svar12 = list_svar[[1]][2,,],
  svar21 = list_svar[[2]][1,,],
  svar22 = list_svar[[2]][2,,],
  svar31 = list_svar[[3]][1,,],
  svar32 = list_svar[[3]][2,,],
  svar41 = list_svar[[4]][1,,],
  svar42 = list_svar[[4]][2,,],
  svarma11 = list_svarma[[1]][1,,],
  svarma12 = list_svarma[[1]][2,,],
  svarma21 = list_svarma[[2]][1,,],
  svarma22 = list_svarma[[2]][2,,],
  svarma31 = list_svarma[[3]][1,,],
  svarma32 = list_svarma[[3]][2,,],
  svarma41 = list_svarma[[4]][1,,],
  svarma42 = list_svarma[[4]][2,,],
  th.51 = th.5[1,,],
  th.52 = th.5[2,,],
  th.91 = th.9[1,,],
  th.92 = th.9[2,,],
  lag = 0:12
)


irf_tbl %>% select(svar42, svarma32, svarma42, th.92, lag) %>% 
  mutate(svar42 = svar42/svar42[3] *irf_tbl$th.92[3],
         svarma32 = svarma32/svarma32[3] *irf_tbl$th.92[3],
         svarma42 = svarma42/svarma42[3] *irf_tbl$th.92[3]) %>% 
  rename("svar" = svar42,
         "svarma1" = svarma32,
         "svarma2" = svarma42,
         "theor" = th.92) %>% pivot_longer(-lag) -> plt_tbl1

plt1 <- plt_tbl1 %>% ggplot(aes(x=lag, y=value, col = name)) + 
  geom_line(aes(linetype=name)) + 
  scale_x_continuous(breaks=plt_tbl1$lag) +
  geom_point(aes(shape=name))
ggsave("~/Documents/Rscripts/sim_news/plt2.pdf", plt1)


