# purpose and the output of the script: ###
# the script 1) simulates data from the model with news shocks
# with different degrees of non-Gaussianity and non-invertibility;
# 2) estimates the SVAR model where the structural shocks are
# identified using non-Gaussianity; 3) saving the dataset and 
# the SVAR impulse response function

# Packages ####
.libPaths(c(.libPaths(), "/proj/juhokois/R/"))
pkgs <- c("lubridate", "xts", "parallel", "svarmawhf", "svars", "fitdistrplus", "sgt", "tidyverse")
void = lapply(pkgs, library, character.only = TRUE)

# Helper functions ####
simu_y = function(model, n.obs, rand.gen = stats::rnorm, n.burnin = 0, ...) {
  
  d = dim(model$sys)
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  
  if ( (m*(p+1)) == 0 ) stop('illegal ARMA system (m=0 or p<0)')
  
  # check 'n.obs'
  n.obs = as.integer(n.obs[1])
  if (n.obs <= 0) stop('"n.obs" must be a positive integer!')
  
  # check 'n.burnin'
  n.burnin = as.integer(n.burnin[1])
  if (n.burnin < 0) stop('"n.burnin" must be a non negative integer!')
  
  # get ARMA parameters
  a = unclass(model$sys$a)
  b = unclass(model$sys$b)
  
  # convert ARMA system
  #    a[0] y[t] + a[1] y[t-1] + ... = b[0] u[t] + b[1] u[t-1] + ...
  # to system of the form
  #    y[t] = a[1] y[t-1] + ... + b[0] u[t] ...
  # and create parameter matrices a = (a[p], ..., a[1]) and b = (b[0], ..., b[q])
  a0 = matrix(a[,,1], nrow = m, ncol = m)
  
  dim(b) = c(m, n*(q+1))
  if ((n*(q+1)) > 0) {
    b = solve(a0, b)
  }
  
  # note for solve_ARMA_cpp we have to reshuffle the AR parameter as:  a = (a[p],...,a[1])
  if (p > 0) {
    a = a[,,(p+1):2, drop = FALSE]
    dim(a) = c(m, m*p)
    a = -solve(a0, a)
  } else {
    a = matrix(0, nrow = m, ncol = 0)
  }
  
  # generate disturbances
  u = matrix(rand.gen((n.obs+n.burnin)*n), nrow = n, ncol = n.obs+n.burnin)
  u = model$sigma_L %*% t(solve(chol(var(t(u[, (n.burnin+1):(n.burnin+n.obs)]))))) %*% u
  # outputs
  y = matrix(0, nrow = m, ncol = n.obs+n.burnin)
  svarmawhf::solve_ARMA_cpp(a, b, u, y, 1)
  
  return(list(y = t(y[, (n.burnin+1):(n.burnin+n.obs), drop = FALSE]),
              u = t(u[, (n.burnin+1):(n.burnin+n.obs), drop = FALSE])))
}
sim_news <- function(beta, rho, nobs, nu){
  theta <- 1/(1-beta*rho)
  ar_pol <- array(c(diag(2),                   # lag 0
                    -rho, -rho*theta, 0, 0),   # lag 1
                  dim = c(2,2,2))
  ma_pol <- array(c(diag(2), # lag 0
                    0, -theta/beta, 0, 1/beta,       # lag 1
                    -1/beta^2, -theta/beta^2, 1/(theta*beta^2), 1/beta^2), # lag 2
                  dim = c(2,2,3))
  bmat <- matrix(c(1, theta, 0, theta*beta^2), 2, 2)
  re_mod <- armamod(sys = lmfd(ar_pol, ma_pol), sigma_L = bmat)
  data_out <- simu_y(model = re_mod,
                     rand.gen =  function(x) stats::rt(x, nu),
                     n.burnin = 2*nobs,
                     n.obs = nobs)
  list(y = data_out, 
       mod = re_mod, 
       prms = c("beta" = beta, "rho" = rho, "nobs" = nobs, "nu" = nu))
}

# Simulation params ####
n_ahead = 12
mc_n <- 101
sim_prm <- expand.grid(beta=c(0.5,0.9), rho = 0.5, nobs = 250, nu = c(12,48))
data_list <- vector("list", mc_n*nrow(sim_prm))
irf_svar <- array(NA, c(2,2,n_ahead+1,mc_n,nrow(sim_prm)))

# Simulation and data save ####
for(prm_ix in 1:nrow(sim_prm)){
  
  mc_ix <- 1
  
  while(mc_ix <= mc_n){
    
    DATASET <- do.call(sim_news, sim_prm[prm_ix,])$y$y
    
    # SVAR ####
    var_obj <- vars::VAR(DATASET, type = "none", lag.max = 12, ic = "SC") %>% 
      suppressWarnings()
    rest_mat <- diag(4)
    rest_mat[3,3] <- 0 
    svar_obj <- svars::id.ngml(var_obj, 
                               stage3 = TRUE,
                               restriction_matrix = rest_mat) %>% try(silent = TRUE)
    
    if(!inherits(svar_obj, 'try-error')){
    
      irf_svar[,,,mc_ix,prm_ix] <- svar_obj %>% 
        vars::irf(n.ahead = n_ahead+1) %>%
        .$irf %>%
        unlist %>%
        .[-c(1:(n_ahead+1))] %>% 
        array(c(n_ahead+1, 2, 2)) %>%
        aperm(c(3,2,1))
      data_list[[(prm_ix-1)*mc_n+mc_ix]] <- DATASET
      mc_ix <- mc_ix + 1
      if(mc_ix%%20==0) cat(paste("MC run:", mc_ix, "w/ prm_set:", prm_ix, "\n"))
    }
  }
}

data_list <- tibble(data_list, mc_ix = rep(1:mc_n, 4), prm_ix = rep(1:4, each = mc_n))
saveRDS(irf_svar, file = "./local_data/irf_svar.rds")
saveRDS(data_list, file = "./local_data/data_list.rds")
