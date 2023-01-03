# Create bootstrap samples for the IRF estimation
# Packages ####
pkgs <- c("svarmawhf", "tidyverse")
void = lapply(pkgs, library, character.only = TRUE)

X_boot <- function(X, L){
  
  X <- as.matrix(X)
  nobs <- nrow(X)
  nvar <- ncol(X)
  K <- floor(nobs/L)
  blks <- ceiling(stats::runif(K)*K)
  X_boot <- matrix(0, nrow = K*L, ncol = nvar)
  for(i in 1:K){
    X_boot[((i-1)*L + 1):(i*L),] <- X[((blks[i]-1)*L + 1):(blks[i]*L),]
  }
  return(X_boot)
}

total_data <- readRDS("~/Documents/Rscripts/sim_news/local_data/total_data.rds")
ds <- total_data %>% 
  filter(mp_type == "GSS22") %>% 
  pull(data_list) %>% .[[1]]

data_list <- vector("list", 1000)
for(j in 1:1000){
  data_list[[j]] <- X_boot(ds, 37)
}
saveRDS(data_list, file = "./local_data/data_list_boot.rds")