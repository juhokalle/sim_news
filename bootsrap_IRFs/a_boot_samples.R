# Create bootstrap samples for the IRF estimation
# Packages ####
pkgs <- c("svarmawhf", "tidyverse")
void = lapply(pkgs, library, character.only = TRUE)

mb_boot <- function(y, prms, tmpl, b.length=10, nboot=500){
  # Most of this is a direct copy of svars::mb.boot
  u <- get_residuals_once(params_deep = prms, tmpl = tmpl, data_long = y)
  armamod <- armamod_whf(theta = prms, tmpl = tmpl)
  obs <- nrow(y)
  k <- ncol(y)
  errors <- list()
  N <- obs/b.length
  blocks <- array(NA, c(b.length, k, obs - b.length + 1))
  
  for (i in 0:(obs - b.length)) {
    blocks[, , (i + 1)] <- u[(i + 1):(i + b.length), ]
  }
  for (i in 1:nboot) {
    epsilon.star <- matrix(0, b.length * (ceiling(N) + 1), 
                           ncol(u))
    epsilon.star <- list()
    for (kk in 1:(ceiling(N) + 1)) {
      epsilon.star[[kk]] <- blocks[, , floor(runif(1, 1, obs - b.length + 2))]
    }
    epsilon.star <- do.call("rbind", epsilon.star)
    for (s in 1:b.length) {
      b.mean <- colSums(epsilon.star[1:(s + (obs - b.length)), ])/(obs - b.length + 1)
      for (j in 0:floor(N)) {
        epsilon.star[j * b.length + s, ] <- epsilon.star[j * b.length + s, ] - b.mean
      }
    }
    epsilon.star <- epsilon.star[1:obs, ]
    errors[[i]] <- epsilon.star
  }
  ystar <- lapply(errors, function(x) solve_de(lmfd(armamod$polm_ar, armamod$polm_ma), x)$y)
  return(ystar)
}

tbl0 <- readRDS("./local_data/target_model.rds")
ds <- readRDS("./local_data/jobid_20230303/total_data_sim.rds") %>% 
  slice(tbl0$nr) %>%  
  pull(data_list) %>% .[[1]]
bl_vec <- c(5,10,20,50)

arg_list <- map(bl_vec, ~
                  list(y = ds, 
                       prms = tbl0 %>% pull(params_deep_final) %>% .[[1]],
                       tmpl = tbl0 %>% pull(tmpl) %>% .[[1]],
                       b.length = .x,
                       nboot = 1000)
                )

dl <- map(arg_list, ~ do.call(mb_boot, .x)) %>% unlist(recursive = FALSE)

# standardize data for estimation
data_list <- tibble(data_list = lapply(dl, function(x) apply(x, 2, function(xx) (xx-mean(xx))/sd(xx))),
                    std_dev = lapply(dl, function(x) apply(x, 2, sd)),
                    mb_length = rep(bl_vec, each = 1000))

saveRDS(data_list, file = "./local_data/data_list_boot.rds")
