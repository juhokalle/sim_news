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

ds <- readRDS("./local_data/total_data.rds") %>% 
  filter(mp_type == "GSS22") %>% 
  pull(data_list) %>% .[[1]]

tbl0 <- readRDS("./local_data/target_model.rds")
arg_list <- lapply(c("tdist", "sgt"), function(x)
  list(y = ds, 
       prms = tbl0 %>% filter(sd==x) %>% pull(params_deep_final) %>% .[[1]],
       tmpl = tbl0 %>% filter(sd==x) %>% pull(tmpl) %>% .[[1]],
       b.length = 10,
       nboot = 1000
       )
  )

dl <- lapply(arg_list, function(x) do.call(mb_boot, x)) %>% 
  unlist(recursive = FALSE)

# standardize data for estimation
data_list <- tibble(data_list = lapply(dl, function(x) scale(x)),
                    std_dev = lapply(dl, function(x) apply(x, 2, sd)),
                    sd = rep(c("tdist", "sgt"), each = 1000))

saveRDS(data_list, file = "./local_data/data_list_boot.rds")