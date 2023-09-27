# SIMULATIONS
simu_y = function(model, n.obs, rand.gen = stats::rnorm, n.burnin = 0, ...)
{
  
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
  u = u - tcrossprod(rowMeans(u[,(n.burnin+1):(n.burnin+n.obs)]), rep(1, n.obs+n.burnin))
  # outputs
  y = matrix(0, nrow = m, ncol = n.obs+n.burnin)
  svarmawhf::solve_ARMA_cpp(a, b, u, y, 1)
  
  return(list(y = t(y[, (n.burnin+1):(n.burnin+n.obs), drop = FALSE]),
              u = t(u[, (n.burnin+1):(n.burnin+n.obs), drop = FALSE])))
}

sim_news <- function(beta, rho, 
                     no_sim = FALSE, nobs = 250, nu = 3, 
                     rg_fun = function(x) stats::rt(x, nu))
{
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
  if(no_sim) return(re_mod)
  data_out <- simu_y(model = re_mod,
                     rand.gen =  rg_fun,
                     n.burnin = 1000,
                     n.obs = nobs)
  return(list(y = data_out, 
              mod = re_mod, 
              prms = c("beta" = beta, "rho" = rho, "nobs" = nobs, "nu" = nu))
         )
}

mixed_marg_dists <- function(dim_out, nu)
{
  function(x) c(replicate(x/dim_out, c(rnorm(2), sqrt(1-2/nu)*stats::rt(dim_out-2, nu))))
}

get_simu_model <- function(mod_tmpl, imp_mat)
{
  root_flag <- TRUE
  scl_prm <- 1
  n_unst <- k_kappa2nunst(q = mod_tmpl$input_orders$MAorder,
                          dim_out = mod_tmpl$input_orders$dim_out,
                          k = mod_tmpl$input_orders$k,
                          kappa = mod_tmpl$input_orders$kappa)
  
  while(root_flag && scl_prm > 0){
    mod_prms <- rnorm(mod_tmpl$n_par, mean = 0, sd = scl_prm)
    arma_whf <- fill_tmpl_whf_rev(mod_prms, mod_tmpl)
    arma_std <- armamod_whf(mod_prms, mod_tmpl)
    Bmat <- diag(diag(imp_mat))
    ma_polm <- arma_std$polm_ma%r%(solve(unclass(arma_std$polm_ma)[,,1])%*%imp_mat%*%solve(Bmat))
    # Check root constraints
    ar_flag <- any(abs(eigen(companion_matrix(arma_whf$polm_ar))$val) > 1)
    ma_bwd_flag <- any(abs(eigen(companion_matrix(arma_whf$polm_ma_bwd))$val) > 1)
    ma_fwd_flag <- any(abs(eigen(companion_matrix(arma_whf$polm_ma_fwd))$val) > 1)
    ma_flag <- !(sum(abs(eigen(companion_matrix(ma_polm))$val)>1)==n_unst)
    root_flag <- ar_flag || ma_bwd_flag || ma_fwd_flag || ma_flag
    scl_prm <- scl_prm-.05
  }
  armamod(lmfd(arma_std$polm_ar, ma_polm), Bmat)
}

# BOOTSTRAP
mb_boot <- function(y, prms, tmpl, b.length=10, nboot=500)
{
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

# SHOCK ID & LABELLING
get_llf <- function(p, q, kappa, k, dtype, sd)
{
  total_data <- readRDS("local_data/total_data.rds")
  data_i <- total_data %>% 
    filter(mp_type==dtype) %>% 
    pull(data_list) %>% .[[1]]
  dim_out <- dim(data_i)[2]
  tmpl_i <- tmpl_whf_rev(dim_out = dim_out,
                         ARorder = p,
                         MAorder = q,
                         kappa = kappa,
                         k = k,
                         shock_distr = "sgt")
  ll_whf_factory(data_wide = t(data_i), tmpl = tmpl_i, shock_distr = sd)
  
}

id_news_shox <- function(irf_arr, policy_var)
{
  dim_in <- dim(irf_arr)[2]
  news_ix <- which.min(abs(unclass(irf_arr)[policy_var,,1]))
  rot_mat <- diag(dim_in)[, c((1:dim_in)[-news_ix], news_ix)]
  irf_tmp <- irf_arr%r%rot_mat
  # Impose positive responses to the shocks: 1) conventional shock(s)
  diag0 <- unclass(irf_tmp)[1:(dim_in-1), 1:(dim_in-1), 1]
  if(dim_in > 2) diag0 <- diag(diag0)
  rot_mat <- rot_mat%*%diag(c(sign(diag0), 1))
  
  # 2) news shock
  peak_ix <- which.max(abs(unclass(irf_tmp)[policy_var, dim_in, -1]))
  if(unclass(irf_tmp)[policy_var, dim_in, peak_ix+1] < 0){
    rot_mat <- rot_mat%*%diag(c(rep(1, dim_in-1), -1))
  }
  rot_mat
}

rotmat <- function(x, n)
{
  rot_ix <- combn(n,2)
  final_mat <- diag(n)
  
  for(jj in 1:ncol(rot_ix)){
    rot_mat <- matrix(c(cos(x[jj]), -sin(x[jj]), sin(x[jj]), cos(x[jj])), 2, 2)
    temp_mat <- diag(n)
    temp_mat[rot_ix[,jj], rot_ix[,jj]] <- rot_mat
    final_mat <- final_mat %*% temp_mat
  }
  final_mat
}

ff <- function(x, zero_ix, input_mat)
{
  row_ix <- if(is.null(dim(zero_ix))) zero_ix[1] else zero_ix[,1]
  col_ix <- if(is.null(dim(zero_ix))) zero_ix[2] else zero_ix[,2]
  sum(sqrt(diag(input_mat[row_ix,] %*% rotmat(x, ncol(input_mat))[, col_ix])^2))
}

id_policy_shox <- function(irf_arr, policy_var, n_shox)
{
  
  dim_out <- dim(irf_arr)[1]
  n_ahead <- dim(irf_arr)[3]
  fevd_obj <- get_fevd(irf_arr, int_var = policy_var)[[1]]
  combs <- combn(dim_out, n_shox)
  res_mt <- matrix(NA, n_ahead, ncol(combs))
  
  for(i in 1:n_ahead){
    for(j in 1:ncol(combs)){
      res_mt[i, j] <- sum(fevd_obj[i, combs[,j]])
    }
  }
  combs[,which.max(colMeans(res_mt))]
}

optim_zr <- function(input_mat, zr_ix, opt_it = TRUE)
{
  n_var <- dim(input_mat)[1]
  
  if(opt_it){
    # chooses the zero restriction based on minimum distance (measured by F-rnom) 
    # from unconstrained estimator
    opt_ix <- rep(0, n_var)
    rotmat_ix <- array(NA, dim = c(n_var, n_var, n_var))
    n_rest <- length(zr_ix)
    for(jj in 1:n_var){
      # UNCOSTRAINED OPTIM: BFGS
      zero_ix = if(n_rest==1) c(zr_ix, jj) else matrix(c(zr_ix, rep(jj, n_rest)), n_rest, n_rest)
      optim(par = rep(0, (n_var^2-n_var)/2),
            fn = ff, 
            zero_ix = zero_ix,
            input_mat = input_mat,
            method = "BFGS") -> opt_obj
      
      rotmat_ix[,,jj] <- rotmat(opt_obj$par, n_var)
      opt_ix[jj] <- norm(input_mat - input_mat%*%rotmat_ix[,,jj], "F")
      return(rotmat_ix[,,which.min(opt_ix)])
    }
  } else {
    # ex ante fixed the zero restriction
    optim(par = rep(0, (n_var^2-n_var)/2),
          fn = ff, 
          zero_ix = zr_ix,
          input_mat = input_mat,
          method = "BFGS") -> opt_obj
    
    return(rotmat(opt_obj$par, n_var))
  }
}

zero_rest <- function(target_mat, zero_ix){
  
  rot_angle <- ifelse(zero_ix[2]==1, 
                      -target_mat[zero_ix[1], zero_ix[2]+1],
                      target_mat[zero_ix[1], zero_ix[2]-1])
  rot_angle <- atan(target_mat[zero_ix[1], zero_ix[2]]/rot_angle)
  sub_mat <- c(cos(rot_angle), sin(rot_angle), -sin(rot_angle), cos(rot_angle))
  rot_mat <- diag(1, dim(target_mat)[2])
  tmp_ix <- ifelse(zero_ix[2]==1, zero_ix[2], zero_ix[2]-1)
  rot_mat[tmp_ix:(tmp_ix+1), tmp_ix:(tmp_ix+1)] <- sub_mat
  rot_mat
}


get_fevd <- function (irf_arr, int_var = NULL, by_arg = NULL)
{
  irf_arr <- unclass(irf_arr)
  nvar <- dim(irf_arr)[1]
  n.ahead <- dim(irf_arr)[3]
  fe <- list()
  if(is.null(int_var)){
    for (i in 1:nvar) {
      fe[[i]] <- as.data.frame(t(irf_arr[i, , ]))
    }
  } else{
    fe[[1]] <- as.data.frame(t(irf_arr[int_var, , ]))
  }
  fe2 <- fe
  seq_grid <- if(is.null(by_arg)) 1:n.ahead else seq(1,n.ahead, by = by_arg)
  for (i in 1:length(fe)) {
    for (j in seq_grid) {
      fe2[[i]][j, ] <- (colSums(fe[[i]][j:1, ]^2)/sum(fe[[i]][j:1, ]^2)) * 100
    }
    fe2[[i]] <- fe2[[i]][seq_grid,]
  }
  return(fe2)
}

norm_irf <- function(irf_arr, 
                     norm_int, 
                     norm_pos=dim(irf_arr)[1], 
                     norm_scm = c("surp", "news"))
{
  
  if(norm_scm=="surp"){
    irf_arr <- irf_arr[,1,,drop=FALSE]/irf_arr[norm_pos,1,1]*norm_int
  } else if(norm_scm=="news"){
    irf_arr <- irf_arr[,1,,drop=FALSE]/max(irf_arr[norm_pos,1,])*norm_int
  }
  irf_arr
}

get_rest_irf <- function(tbl_slice, ...)
{
  nobs <- tbl_slice$nobs
  sd_mat <- readRDS("local_data/svarma_data_list.rds") %>% 
    filter(mp_type%in%tbl_slice$mp_type) %>%
    pull(std_dev) %>% .[[1]] %>%  diag
  
  irf_out <- sd_mat %r% tbl_slice$irf[[1]]
  
  opt_obj <- optim_zr(irf_out[,,1], ...)
  llf0 <- do.call(get_llf, tbl_slice %>% 
                    select(p, q, kappa, k, mp_type, sd) %>% 
                    rename(dtype = mp_type)
  )
  params0 <- tbl_slice$params_deep_final[[1]]
  
  params0[tbl_slice$params_deep_final[[1]] %in% tbl_slice$B_mat[[1]]] <- c(tbl0$B_mat[[1]]%*%opt_obj$rmat)
  pval <- 1-pchisq(2*nobs*(llf0(params0)-tbl0$value_final), opt_obj$n_rest)
  
  irf_out <- irf_whf(params0, tbl_slice$tmpl[[1]], 48)
  
  return(list(irf = irf_out, pval = pval, rmat = opt_obj$rmat))
}

choose_perm_sign <- function(target_mat, cand_mat, type = c("frob", "dg_abs", "min_rmse"))
{
  nvar <- ncol(cand_mat)
  sign_ix <- nvar %>% replicate(list(c(-1,1))) %>% expand.grid
  perm_ix <- nvar %>% 
    replicate(list(1:nvar)) %>% 
    expand.grid %>% 
    filter(apply(., 1, n_distinct)==nvar)
  cr0 <- 1e25
  for(j in 1:nrow(sign_ix)){
    for(jj in 1:nrow(perm_ix)){
      rt_mat <- diag(sign_ix[j, ]) %*% diag(nvar)[, unlist(perm_ix[jj, ])]
      x1 <- if(type=="min_rmse") cand_mat %r% rt_mat else cand_mat %*% rt_mat
      if(type%in%c("frob", "min_rmse")){
        cr1 <- sqrt(mean(c(target_mat - x1)^2))
      } else if(type=="dg_abs"){
        cr1 <- if(all(diag(x1)>0)) -abs(prod(diag(x1))) else 1e25
      }
      if(cr1 < cr0){
        x_opt <- x1
        rt_opt <- rt_mat
        cr0 <- cr1
      }
    }
  }
  rt_opt
}

id_sign <- function(chol_irf, ndraws=1e3, sign_rest, max_draws=1e5){
  irf_dim <- if(length(dim(chol_irf))==2) c(dim(chol_irf), 1) else dim(chol_irf)
  dim(chol_irf) <- irf_dim
  nvar <- irf_dim[1]
  perm_ix <- nvar %>% 
    replicate(list(1:nvar)) %>% 
    expand.grid %>% 
    filter(apply(., 1, n_distinct)==nvar) %>% 
    as.matrix()
  imp_arr <- array(data = NA, dim = c(irf_dim, ndraws))
  i <- 1
  loop_counter <- 1
  while(i <= ndraws && loop_counter <= max_draws){
    qr_obj <- qr(matrix(rnorm(nvar^2), nvar, nvar))
    qr_q <- qr.Q(qr_obj)
    qr_r <- qr.R(qr_obj)
    rot_comp <- FALSE
    j <- 1
    while(!rot_comp){
      ort_mat <- qr_q%*%diag(sign(diag(qr_r)))%*%diag(nvar)[, perm_ix[j,]]
      imp_arr[,,,i] <- apply(chol_irf, 3, function(x) x%*%ort_mat)
      sign_ix <- which(is.finite(sign_rest))
      rot_ok <- min(diag(sign_rest[sign_ix])%*%imp_arr[,,,i][sign_ix])>0
      rot_comp <- rot_ok || j == nrow(perm_ix)
      j <- j + 1
    }
    if(rot_ok) i <- i + 1
    loop_counter <- loop_counter + 1
  }
  if(ndraws>i){
    stop("Failed to find the requested number of impulse responses matching the sign restrictions.")
  } else{
    drop(imp_arr)
  }
}

id_mixed_old <- function(chol_irf,
                         rest_mat,
                         news_rest = NULL,
                         ndraws = 1e3,
                         max_draws = ndraws*100,
                         irf_cb = c(0.14, 0.86),
                         replace_md = TRUE,
                         verbose = FALSE){
  # Extract IRF dimension
  irf_dim <- dim(chol_irf)
  # Number of variables
  nvar <- irf_dim[1]
  # Initialize output to arrays: IRFs and resulting orthogonal matrices
  imp_arr <- array(NA, c(irf_dim, ndraws))
  omat_arr <- array(NA, c(nvar, nvar, ndraws))
  # Initialization
  search_complete <- FALSE
  i <- 1
  max_i <- 1
  # Print progress bar
  if(verbose) cat("Search in progress")  
  while(!search_complete){
    # Generate random orth. matrix using QR decomp
    qr_obj <- qr(matrix(rnorm(nvar^2), nvar, nvar))
    qr_q <- qr.Q(qr_obj)
    qr_r <- qr.R(qr_obj)
    ort_mat <- qr_q%*%diag(sign(diag(qr_r)))
    # Calculate structural IRFs
    imp_arr[, , , i] <- apply(chol_irf, 3, function(x) x%*%ort_mat)
    # Obtain zero restrictions
    zr_ix <- which(rest_mat==0, arr.ind = TRUE)[,1:2]
    if(length(zr_ix)!=0){
      # Impose zero restrictions via Givens rotation
      ort_mat <- ort_mat%*%zero_rest(imp_arr[, , 1, i], zr_ix)
      # Update IRF output
      imp_arr[, , , i] <- apply(chol_irf, 3, function(x) x%*%ort_mat)
      # Update orth. matrix output
      omat_arr[, , i] <- ort_mat
    }
    # Index of elements with only sign restrictions, i.e. length of the tested prms
    sign_ix <- which(rest_mat!=0)
    # Index of allowed column permutations 
    no_zr_ix <- apply(rest_mat, 2, function(x) !(0 %in% x))
    perm_sign_list <- get_perm_sign_mat(nvar, which(no_zr_ix))
    # Initialize signed column permutation search loop: 
    perm_sign_complete <- FALSE 
    j <- 1
    while(!perm_sign_complete){
      # Sign and column reversal
      ort_mat <- ort_mat%*%perm_sign_list[[j]]
      # Update IRF output
      imp_arr[, , , i] <- apply(chol_irf, 3, function(x) x%*%ort_mat)
      # Update orth. matrix output
      omat_arr[, , i] <- ort_mat
      # check constraints
      rot_ok <- min(diag(rest_mat[sign_ix])%*%imp_arr[, , , i][sign_ix])>0
      if(rot_ok && !is.null(news_rest)){
        # checks whether the sign restriction wrt news shock is satisfied 
        sign_news <- news_rest[[1]] # save sign response to a news shock
        # Extract the impulse responses of target variable to news shock
        news_resp <- imp_arr[news_rest[[2]][1], news_rest[[2]][2], ,i]
        # Check response to news shock 
        rot_ok <- rot_ok && sign(news_resp[which.max(abs(news_resp))]) == sign_news
      }
      perm_sign_complete <- rot_ok || j == length(perm_sign_list)
      j <- j + 1
    }

    if(rot_ok){
      # Increase loop index if restrictions satisfied
      i <- i + 1
      if(verbose){
        # progress bar
        grid_check <- seq((ndraws/20), ndraws, by=ndraws/20)
        if(i %in% grid_check) cat(".")
        if(i %in% grid_check[seq(5, 20, by=5)]) cat(paste0(100*i/ndraws, "%"))
        if(i == grid_check[length(grid_check)]) cat("\n")
      }
    }
    # Update maximum draws counter
    max_i <- max_i + 1
    # Boolean for outer loop
    search_complete <- i > ndraws || max_i > max_draws
  }
  if(ndraws >= i) warning("Failed to find the requested number of 
                          impulse responses matching the sign restrictions.")
  if(i==1){
    list_out <- "No rotations found."
  } else{
    # Collect output, note that 'i - 1' returns the true number of draws 
    # matching the restrictions
    list_out <- list(irf = drop(imp_arr[, , , 1:(i - 1)]), 
                     bmat = omat_arr[, , 1:(i - 1)],
                     draws_ok = i - 1,
                     draws_total = max_i - 1) 
    # Append the output list with confidence bands if applicable
    if(!is.null(irf_cb)){
      list_out$irf_qt <- apply(X = imp_arr, 
                               MARGIN = c(1, 2, 3), 
                               FUN = quantile, 
                               probs = c(irf_cb[1], 0.5, irf_cb[2]),
                               na.rm = TRUE,
                               ) %>% 
        aperm(c(2, 3, 4, 1))
      if(replace_md) list_out$irf_qt[,,,2] <- fry_pagan_mt(list_out$irf)
    }
  }
  list_out
}

id_mixed_new <- function(chol_irf,
                         rest_mat,
                         news_rest = NULL,
                         ndraws = 1e3,
                         max_draws = ndraws*100,
                         irf_cb = c(0.14, 0.86),
                         replace_md = TRUE,
                         verbose = FALSE){
  # Extract IRF dimension
  irf_dim <- dim(chol_irf)
  # Number of variables
  nvar <- irf_dim[1]
  # Initialize output to arrays: IRFs and resulting orthogonal matrices
  imp_arr <- array(NA, c(irf_dim, ndraws))
  omat_arr <- array(NA, c(nvar, nvar, ndraws))
  
  # Subsetting index of elements with only sign restrictions, i.e. length of the tested prms
  sign_ix <- which(rest_mat!=0)
  
  # Index of allowed column permutations 
  no_zr_ix <- apply(rest_mat, 2, function(x) !(0 %in% x))
  perm_sign_arr <- get_perm_sign_mat(nvar, which(no_zr_ix)) %>% 
    abind::abind(along = 3)

  # Obtain zero restrictions
  zr_ix <- which(rest_mat==0, arr.ind = TRUE)[,1:2]

  # Initialization
  search_complete <- FALSE
  i <- 1
  max_i <- 1
  
  # Print progress bar
  if(verbose) cat("Search in progress")  
  while(!search_complete){
    # Generate random orth. matrix using QR decomp
    qr_obj <- qr(matrix(rnorm(nvar^2), nvar, nvar))
    qr_q <- qr.Q(qr_obj)
    qr_r <- qr.R(qr_obj)
    # Ensure uniqueness
    ort_mat <- qr_q%*%diag(sign(diag(qr_r)))
    
    # Multiply by Givens matrix if zero restricitons imposed
    if(length(zr_ix)!=0) ort_mat <- ort_mat%*%zero_rest(chol_irf[, , 1]%*%ort_mat, zr_ix)
    
    # Collect arguments for cpp function
    arg_list <- list(chol_irf = chol_irf, 
                     ort_mat = ort_mat, 
                     rest_arr = rest_mat, 
                     sign_ix = sign_ix, 
                     perm_sign_arr = perm_sign_arr,
                     news_info = if(is.null(news_rest)) rep(0, 3) else news_rest)
    search_rotations <- do.call(search_rotations_cpp, arg_list)
    
    if(search_rotations$success){
      # Save results 
      omat_arr[, , i] <- search_rotations$omat
      imp_arr[, , , i] <- search_rotations$irf_tmp
      # Increase loop index if restrictions satisfied
      i <- i + 1
      if(verbose){
        # progress bar
        grid_check <- seq((ndraws/20), ndraws, by=ndraws/20)
        if(i %in% grid_check) cat(".")
        if(i %in% grid_check[seq(5, 20, by=5)]) cat(paste0(100*i/ndraws, "%"))
        if(i == grid_check[length(grid_check)]) cat("\n")
      }
    }
    # Update maximum draws counter
    max_i <- max_i + 1
    # Boolean for outer loop
    search_complete <- i > ndraws || max_i > max_draws
  }
  if(ndraws >= i) warning("Failed to find the requested number of 
                          impulse responses matching the sign restrictions.")
  if(i==1){
    list_out <- "No rotations found."
  } else{
    # Collect output, note that 'i - 1' returns the true number of draws 
    # matching the restrictions
    list_out <- list(irf = drop(imp_arr[, , , 1:(i - 1)]), 
                     bmat = omat_arr[, , 1:(i - 1)],
                     draws_ok = i - 1,
                     draws_total = max_i - 1) 
    # Append the output list with confidence bands if applicable
    if(!is.null(irf_cb)){
      list_out$irf_qt <- apply(X = imp_arr, 
                               MARGIN = c(1, 2, 3), 
                               FUN = quantile, 
                               probs = c(irf_cb[1], 0.5, irf_cb[2]),
                               na.rm = TRUE
                               ) %>% 
        aperm(c(2, 3, 4, 1))
      if(replace_md) list_out$irf_qt[,,,2] <- fry_pagan_mt(list_out$irf)
    }
  }
  list_out
}

Rcpp::cppFunction(depends = "RcppArmadillo",
  'List search_rotations_cpp(
  const arma::cube& chol_irf,
  const arma::mat& ort_mat,
  const arma::cube& rest_arr, 
  const arma::uvec& sign_ix,
  const arma::cube& perm_sign_arr, 
  const arma::uvec& news_info) {
  
  // preallocation
  arma::cube irf_temp(size(chol_irf));
  arma::mat ort_temp(size(ort_mat)), dg_rest(sign_ix.n_elem, sign_ix.n_elem);
  arma::colvec irf_vec(chol_irf.n_elem), irf_vec_sub(sign_ix.n_elem);
  arma::colvec news_resp(chol_irf.n_slices), abs_news_resp(chol_irf.n_slices);
  double news_rest = news_info[0], row_ix = news_info[1]-1, col_ix = news_info[2]-1;
  
  // subset restrictions and place them on main diagonal of diagonal matrix
  arma::colvec rest_vec = arma::vectorise(rest_arr);
  dg_rest.diag(0) = rest_vec.elem(sign_ix-1);

  // initialization
  bool rot_ok = false, perm_sign_complete = false;
  int i = 0;
  
  while(!perm_sign_complete){
  	// calculate rotated IRFs
  	ort_temp = ort_mat*perm_sign_arr.slice(i);
  	for(int j = 0; j < chol_irf.n_slices; j++){
  		irf_temp.slice(j) = chol_irf.slice(j)*ort_temp;
  	}
  	// vectorize restrictions and irf array
  	irf_vec = arma::vectorise(irf_temp);
  	// subset restricted elements into column vector
  	irf_vec_sub = irf_vec.elem(sign_ix-1);
  	// check restrictions
  	rot_ok = min(dg_rest*irf_vec_sub) > 0;
  	
  	// news restriction
  	if(sum(news_info)!=0 && rot_ok){
      news_resp = irf_temp.subcube(arma::span(row_ix), arma::span(col_ix), arma::span());
  		abs_news_resp = abs(news_resp);
  		double max_ix = abs_news_resp.index_max();
  		rot_ok = as_scalar(news_resp.row(max_ix)*news_rest) > 0;
  	 }
  	
  	// check if the search for signed rotation is complete or continues
  	perm_sign_complete = rot_ok || i == perm_sign_arr.n_slices-1;
  	i++;
  }

  return List::create(Named("success") = rot_ok, Rcpp::Named("omat") = ort_temp, Rcpp::Named("irf_tmp") = irf_temp);
  }'
)

# RESULTS: FIGS & TBLS
plot_irf <- function(irf_arr, var_name = NULL, shock_name = NULL, date_break = 12)
{
  n_var <- dim(irf_arr)[1]
  if(is.null(var_name)) var_name <- letters[1:n_var]
  n_shock <- dim(irf_arr)[2]
  if(is.null(shock_name)) shock_name <- paste0("e_", 1:n_shock)
  n_ahead <- dim(irf_arr)[3]-1
  if(length(dim(irf_arr))==4){
    n_qt <- dim(irf_arr)[4]
    irf_tbl <- map(1:n_qt, ~ c(irf_arr[,,,.x]))
    names(irf_tbl) <- paste0("irf", 1:n_qt)
    irf_tbl <- as_tibble(irf_tbl)
  } else{
    n_qt <- 1
    irf_tbl <- tibble(irf = c(irf_arr))
  }
  
  irf_tbl %>% 
    bind_cols(variable = var_name %>% 
                factor(levels = var_name) %>% 
                rep((n_ahead+1)*n_shock),
              shock = shock_name %>% 
                factor(levels=shock_name) %>% 
                rep(each=n_var) %>% 
                rep(n_ahead+1),
              months = rep(0:n_ahead, each = n_var*n_shock)) %>% 
    pivot_longer(starts_with("irf")) %>% 
    ggplot(aes(x=months, y=value)) +
    geom_line(aes(linetype = name)) + 
    scale_linetype_manual(values = abs(-(n_qt%/%2):(n_qt%/%2)) + 1) + #aes(linetype = name)) +
    geom_hline(yintercept = 0, size = 0.15) +
    facet_grid(variable ~ shock, scales = "free_y") +
    facet_rep_grid(variable ~ shock, scales = "free_y", repeat.tick.labels = 'left') +
    scale_x_continuous(breaks = seq(date_break, n_ahead, by = date_break), expand = c(0,0)) +
    scale_y_continuous(n.breaks = 5) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 6, angle = 0),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=6),
          panel.grid.major = element_blank(),
          panel.spacing = unit(.5, "lines"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          axis.line = element_line(colour = "black"))
}

get_stats_tbl <- function(x)
{
  nvar <- ncol(x)
  res <- list()
  for(j in 1:nvar){
    sw_test <- shapiro.test(x[,j])
    jb_test <- tsoutliers::JarqueBera.test(x[,j])
    lb_lvl <- Box.test(x[,j], lag = 24, type = "Ljung-Box")
    lb_abs <- Box.test(abs(x[,j]), lag = 24, type = "Ljung-Box")
    lb_sq <- Box.test(x[,j]^2, lag = 24, type = "Ljung-Box")
    res[[j]] <- c(sw_test$statistic,
                  sw_test$p.value,
                  jb_test[[1]]$statistic,
                  jb_test[[1]]$p.value, 
                  lb_lvl$statistic,
                  lb_lvl$p.value, 
                  lb_abs$statistic, 
                  lb_abs$p.value,
                  lb_sq$statistic,
                  lb_sq$p.value) %>% round(3)
    res[[j]][!1:length(res[[j]])%%2] <- paste0("(", res[[j]][!1:length(res[[j]])%%2], ")")
  }
  res <- reduce(res, cbind)
  rownames(res) <- c("SW-test", "", "JB-test", "",
                     "LB 1", "", 
                     "LB 2", "",
                     "LB 3", "")
  res_ltx <- xtable::xtable(res)
  res_ltx
}

get_cor_tbl <- function(x)
{
  comb_i <- combn(ncol(x),2)
  l_corr <- cor(x, method = "pearson") %>% round(3)
  r_corr <- cor(x, method = "spearman") %>% round(3)
  for(j in 1:ncol(comb_i)){
    
    l_corr[comb_i[2, j], comb_i[1, j]] <- paste0("(", 
                                                 cor.test(x[, comb_i[1,j]], 
                                                          x[, comb_i[2,j]], 
                                                          method = "pearson")$p.value %>% round(3),
                                                 ")")
    r_corr[comb_i[2, j], comb_i[1, j]] <- paste0("(",
                                                 cor.test(x[, comb_i[1,j]], 
                                                          x[, comb_i[2,j]], 
                                                          method = "spearman")$p.value %>% round(3),
                                                 ")")
    
  }
  list(noquote(cbind(l_corr, r_corr)), xtable::xtable(cbind(l_corr, r_corr)))
}

get_qqplots <- function(x)
{
  plot_list <- list()
  for(j in 1:ncol(x)){
    df <- data.frame(y = x[,j])
    p <- ggplot(df, aes(sample = y))
    p <- p + 
      stat_qq() + 
      stat_qq_line() +
      ylab("") +
      xlab("") +
      if(j==1){
        ggtitle(bquote(u[1]))
      } else if (j==2){
        ggtitle(bquote(u[2]))
      } else if (j==3){
        ggtitle(bquote(u[3]))
      } else if (j==4){
        ggtitle(bquote(u[4]))
      }
    plot_list[[j]] <- p + theme(plot.title = element_text(hjust = 0.5))
  }
  plot_list
}

sim_comparison <- function(irf_arr, shock_ix, qntl, prms)
{
  n_ahead <- dim(irf_arr)[3]-1
  res_final <- expand_grid(prms, model = c("VAR", "VARMA"))
  res_final$irf <- t(expand_grid(1:nrow(prms), 1:2)) %>% 
    data.frame() %>% 
    as.list %>% 
    map(~apply(X = irf_arr[shock_ix[1], shock_ix[2],, .x[1],, .x[2]],
               MARGIN = 1,
               FUN = function(x) quantile(x, qntl)))
  res_final <- res_final %>%
    mutate(irf_lb = map(irf, ~.x[1,]),
           irf_ub = map(irf, ~.x[2,])) %>% 
    dplyr::select(-irf) %>%
    unnest_longer(starts_with("irf")) %>% 
    mutate(lag = rep(0:n_ahead, n()/(n_ahead+1)))
  
  res_final <- tibble(prms) %>% 
    mutate(irf_md = map(.x = beta, 
                        ~with(sim_news(beta = .x, rho = 0.5, no_sim = TRUE), 
                              pseries(sys, n_ahead)%r%sigma_L) %>%
                          unclass %>% 
                          .[shock_ix[1],shock_ix[2],]
    )
    ) %>% 
    unnest_longer("irf_md") %>% 
    mutate(lag = rep(0:n_ahead, n()/(n_ahead+1)),
           model = "MOD_0") %>% 
    bind_rows(res_final)
  
  plt_list <- list()
  for(jj in 1:nrow(prms))
  {
    dt_plt <- res_final %>% filter(beta==prms$beta[jj], nu==prms$nu[jj])
    dt_plt <- tibble(irf0 = dt_plt %>% filter(model=="MOD_0") %>% pull(irf_md),
                     varma_lb = dt_plt %>% filter(model=="VARMA") %>% pull(irf_lb),
                     varma_ub = dt_plt %>% filter(model=="VARMA") %>% pull(irf_ub),
                     var_lb = dt_plt %>% filter(model=="VAR") %>% pull(irf_lb),
                     var_ub = dt_plt %>% filter(model=="VAR") %>% pull(irf_ub),
                     lag = dt_plt %>% filter(model=="MOD_0") %>% pull(lag))
    plt_list[[jj]] <- dt_plt %>% ggplot(aes(x=lag)) + 
      geom_line(aes(y=irf0)) + 
      geom_line(aes(y=varma_lb), linetype="dashed", size = 0.25) +
      geom_line(aes(y=varma_ub), linetype="dashed", size = 0.25) +
      geom_line(aes(y=var_lb), linetype="dotted", size = 0.25) +
      geom_line(aes(y=var_ub), linetype="dotted", size = 0.25) +
      geom_ribbon(aes(ymin=varma_lb, ymax = varma_ub), alpha = 0.1, fill="blue") + 
      geom_ribbon(aes(ymin=var_lb, ymax = var_ub), alpha = 0.1, fill = "red") + 
      scale_x_continuous(expand = c(0,0)) +
      ggtitle(bquote(.(paste("(", letters[c(1,3,2,4)][jj], "): ", "")) ~
                       beta == .(paste(prms$beta[jj], ",", sep = "")) ~
                       m == .(prms$nu[jj])
      )
      ) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")) +
      ylab("") + 
      xlab("")
    
  }
  plt_list
}

# OTHER
pmap_tmpl_whf_rev = function(dim_out = DIM_OUT, p, q, kappa, k, shock_distr = "gaussian", ...)
{
  tmpl_whf_rev(dim_out = DIM_OUT, ARorder = p, MAorder = q, kappa = kappa, k = k, shock_distr = shock_distr)
}

pmap_get_residuals_once = function(params_deep_final, tmpl, data_list, ...)
{
  get_residuals_once(params_deep = params_deep_final, tmpl = tmpl, data_long = data_list)
}

hlp_parallel = function(list_input)
{
  return(create_results_list(theta_init = list_input[[1]], 
                             tmpl       = list_input[[2]],
                             params     = params, 
                             DATASET    = list_input[[3]],
                             shock_distr = list_input[[4]]))
}

k_kappa2nunst <- function(q, dim_out, k, kappa){
  
  kk <- c(k, kappa)
  hlp_tbl <- tibble(n_unst = map(q, ~0:(dim_out*.x))) %>% 
    unnest(n_unst) %>% 
    mutate(n_st = dim_out * q - n_unst) %>% 
    mutate(kappa = n_unst %/% dim_out,
           k = n_unst %% dim_out)
  hlp_tbl %>% 
    filter(k == kk[1], kappa == kk[2]) %>% 
    dplyr::select(n_unst) %>% 
    as.double
}

# Tools related to structrual models
get_struc_mat <- function(model_type = c("dynamic", "static", "static_small", "r_smooth"), 
                          param_list = NULL){

  if(model_type == "dynamic"){
    
    if(is.null(param_list)) param_list <- list(beta = .99, kappa = .05, gamma = .5,
                                               delta_x = .1, alfa = .5, 
                                               tau_pi = 1.8, tau_x = .5, tau_r = .6, 
                                               rho_pi = .5, rho_x = .5, rho_r = .5)
    param_vec <- c("beta", "alfa", "kappa", "gamma", "delta_x", 
                   "tau_x", "tau_pi", "tau_r", 
                   "rho_x", "rho_pi", "rho_r")
    check_vec <- param_vec%in%names(param_list)
    if(!all(check_vec)) stop("The following parameters are missing: ", 
                             paste(param_vec[!check_vec], collapse = ", "))
    
    Kmat <- with(param_list, matrix(c(1, -kappa, -tau_x*(1-tau_r),
                                      0, 1, -tau_pi*(1-tau_r),
                                      delta_x, 0, 1), 3, 3))
    
    Amat <- with(param_list, matrix(c(1-gamma, 0, 0,
                                      0, alfa/(1+alfa*beta), 0,
                                      0, 0, tau_r), 3, 3))
    Bmat <- with(param_list, matrix(c(gamma, 0, 0,
                                      delta_x, beta/(1+alfa*beta), 0,
                                      0, 0, 0), 3, 3))
    Hmat <- diag(3)
    
    Dmat <- with(param_list, diag(c(rho_x, rho_pi, rho_r)))
    
  } else if(model_type == "static"){
    
    if(is.null(param_list)) param_list <- list(beta = 0.995, kappa = 0.005,
                                                phi_pi = 1.5, phi_y = 0.1,
                                                sigma_d = 1.6, sigma_s = 0.95, sigma_m = 0.23)
    param_vec <- c("beta", "kappa",
                   "phi_pi", "phi_y",
                   "sigma_d", "sigma_s", "sigma_m")
    check_vec <- param_vec%in%names(param_list)
    if(!all(check_vec)) stop("The following parameters are missing: ", 
                             paste(param_vec[!check_vec], collapse = ", "))
    
    Kmat <- with(param_list, matrix(c(1, -kappa, -phi_y,
                                      0, 1, -phi_pi,
                                      1, 0, 1), 3, 3))
    
    Amat <- matrix(0, 3, 3)
    
    Bmat <- with(param_list, matrix(c(1, 0, 0,
                                      1, beta, 0,
                                      0, 0, 0), 3, 3))
    
    Hmat <- with(param_list, diag(c(sigma_d, sigma_s, sigma_m)))
    
    Dmat <- matrix(0, 3, 3)
  
  } else if(model_type == "static_small"){
    
    if(is.null(param_list)) param_list <- list(sigma = .4, gamma = .75, psi = 2, beta = 0.995,
                                               sigma_d = 1.6, sigma_s = 0.95, sigma_m = 0.23)
    
    param_vec <- c("sigma", "gamma",
                   "psi", "beta",
                   "sigma_d", "sigma_s", "sigma_m")
    
    check_vec <- param_vec%in%names(param_list)
    
    if(!all(check_vec)) stop("The following parameters are missing: ", 
                             paste(param_vec[!check_vec], collapse = ", "))
    
    
    Kmat <- with(param_list, matrix(c(1, -gamma, 0,
                                      0, 1, -psi,
                                      sigma, 0, 1), 3, 3))
    Amat <- matrix(0, 3, 3)
    
    Bmat <- with(param_list, matrix(c(1, 0, 0,
                                      sigma, beta, 0,
                                      0, 0, 0),
                                    3, 3)
                 )
    
    Hmat <- with(param_list, diag(c(sigma_d, sigma_s, sigma_m)))
    
    Dmat <- matrix(0,3,3)
    
  } else if(model_type == "r_smooth"){
    
    if(is.null(param_list)) param_list <- list(beta = .99, kappa = .75, tau = 2.08^-1,
                                               psi_pi = 2.19, psi_x = 0.3, rho_R = 0.84,
                                               sigma_g = .21, sigma_z = 1.16, sigma_R = 0.24)
    
    param_vec <- c("tau", "kappa", "rho_R", "beta", 
                   "psi_x", "psi_pi",
                   "sigma_g", "sigma_z", "sigma_R")
    check_vec <- param_vec%in%names(param_list)
    if(!all(check_vec)) stop("The following parameters are missing: ", 
                             paste(param_vec[!check_vec], collapse = ", "))
    
    Kmat <- with(param_list, matrix(c(1, -kappa, -(1-rho_R)*psi_x,
                                      0, 1, 0,
                                      tau, 0, 1), 3, 3))
    
    Amat <- with(param_list, matrix(c(rep(0,8), rho_R), 3, 3))
    
    Bmat <- with(param_list, matrix(c(1, 0, 0,
                                      tau, beta, (1-rho_R)*psi_pi, 
                                      0, 0, 0), 3, 3))
    
    Hmat <- with(param_list, matrix(c(sigma_g, 0, 0, 
                                      0, -kappa*sigma_z, -(1-rho_R)*psi_x*sigma_z,
                                      0, 0, sigma_R), 3, 3))
    Dmat <- matrix(0, 3, 3)
    
  } else{
    stop("Supply a valid model type.")
  }
  
  list(Kmat=Kmat, Amat=Amat, Bmat=Bmat, Hmat=Hmat, Dmat=Dmat)
}

# The system is of the form
# K * x[t] = A * x[t-1] + B * E(x[t+1]|I[t])
#                + H * w[t]
# w[t] = D * w[t-1] + v[t]

solve_re_mod_bp <- function(Kmat, Amat, Bmat, Hmat, Dmat, eps_val = 1e-9){
  
  # Transform System to Canonical Form:
  # x[t] = Ahat * x[t-1] + Bhat * E(x[t+1]|I[t]) 
  #        + Hhat * w[t]
  Ahat = solve(Kmat, Amat)
  Bhat = solve(Kmat, Bmat)
  Hhat = solve(Kmat, Hmat)
  
  dim1 = dim(Amat)[1]
  dim2 = dim(Amat)[2]
  
  # Compute Matrix C Using Brute-Force Iterative Procedure
  Cmat = diag(dim1)       # Initial Condition
  Fmat = diag(dim1)       # Initial Condition
  iter <- 1
  while(iter == 1 || crit1 >= eps_val || crit2 >= eps_val){
    Fi = solve(diag(dim1)-Bhat%*%Cmat, Bhat)
    Ci = solve(diag(dim1)-Bhat%*%Cmat, Ahat)
    crit1 = max(abs(Fi-Fmat))
    crit2 = max(abs(Ci-Cmat))
    Cmat = Ci 
    Fmat = Fi
    iter = iter+1
    if(iter > 500){ 
      stop("The brute-force iterative procedure did not converge after 500 iterations.")
    }
  }
  
  Gmat <- solve(diag(dim1)-Bhat%*%Cmat, Hhat)
  Pmat <- matrix(solve(diag(dim1^2) - t(Dmat)%x%Fmat, c(Gmat)), dim1, dim1)
  var_polm <- abind::abind(diag(dim1),
                           -Pmat%*%Dmat%*%solve(Pmat) - Cmat, 
                           Pmat%*%Dmat%*%solve(Pmat)%*%Cmat,
                           along = 3)
  armamod(lmfd(a = polm(var_polm), b = diag(dim1)), sigma_L = Pmat)
}

plot_irfs <- function(){
  for(prm_ix in 1:nrow(prms)){
    for(var_ix in 1:2){
      for(shock_ix in 1:2){
        irf_arr <- tt_full %>%
          filter(beta == prms[prm_ix,1],
                 nu == prms[prm_ix,2],
                 q==2) %>% 
          # mutate(value_aic = map2_dbl(.x = value_final, .y = tmpl, ~  .x + .y$n_par * 2/250)) %>%
          # group_by(q) %>%
          # slice_min(value_aic) %>%
          # ungroup() %>% 
          # filter(mb_length == prms[prm_ix]) %>%
          mutate(irf = map(.x=irf,~.x %>% unclass)) %>%
          pull(irf) %>%
          abind::abind(along=4)
        irf_qt <- apply(X = irf_arr[var_ix,shock_ix,,],
                        1, 
                        quantile, 
                        probs= c(.05, .14, .5, .86, .95))
        n_quantiles <- 1:dim(irf_qt)[1]
        n_ahead <- dim(irf_qt)[2]-1
        plot(0:n_ahead,
             irf_qt[median(n_quantiles),], 
             type = "l", 
             ylim = c(min(irf_qt),
                      max(irf_qt)),
             # main = paste0("mb_len = ", prms[prm_ix])
             main = paste("beta: ", prms[prm_ix, 1],
                          " nu: ", prms[prm_ix, 2])
        )
        lines(0:n_ahead, irf_qt[1,], col = "red")
        lines(0:n_ahead, irf_qt[max(n_quantiles),], col = "red")
        if(length(n_quantiles)==5){
          lines(0:n_ahead, irf_qt[2,], col = "blue")
          lines(0:n_ahead, irf_qt[4,], col = "blue")
        }
        abline(h=0, lty = "dashed")
      }
    }
  }
}

get_res_app1 <- function(res_tbl)
{
  nobs <- sort(unique(res_tbl$nobs))
  row_1 <- paste0("\\multicolumn{2}{c}{$T=", nobs, "$}", collapse = " & ")
  col_ix <- 2 + seq(1, 2*length(nobs), by = 2)
  row_2 <- paste0("\\cmidrule(lr){", col_ix, "-", col_ix+1, "}")
  row_3 <- paste0(rep(c(" SVARMA ", " SVAR "), length(nobs)), collapse ="&")
  tmp_tbl <- res_tbl %>% 
    pivot_wider(names_from = c(nobs, p),
                values_from = mad_md) %>% 
    mutate(beta = factor(beta),
           nu = factor(nu)) %>% 
    mutate_if(is.numeric, round, digits=3) %>% 
    as.matrix
  rest_str <- sapply(1:length(tmp_tbl), 
                     function(x) paste(t(tmp_tbl)[x], 
                                       ifelse(test = x%%ncol(tmp_tbl), 
                                              yes = " & ", 
                                              no = "\\\\ \n \\addlinespace \n"
                                              )
                                       )
                     ) %>% 
    paste(collapse = "")
  
  cat("\\begin{table} \n",
      "\\centering \n",
      paste0("\\begin{tabular}{l*{", 2+2*length(nobs), "}{c}} \n"),
      "\\toprule \n",
      "& & ", row_1, "\\\\ \n",
      row_2, "\n",
      "$\\beta$ & $m$ &", row_3, "\\\\ \n",
      "\\midrule \n",
      rest_str,
      "\\bottomrule \n",
      "\\end{tabular} \n",
      "\\end{table}"
  )
}

get_perm_sign_mat <- function(nvar, free_cols = 1:nvar){
  
  # Create a list of size 2^nvar containing all distinct matrices with 
  # 1 and -1 on the diagonal 
  sign_mat_list <- nvar %>% 
    replicate(list(c(1,-1))) %>% 
    expand.grid %>% 
    apply(1, function(x) diag(x), simplify = FALSE)
  
  # Helper matrix
  perm_col_mat <- matrix(1:nvar, nvar, factorial(length(free_cols)))
  
  # Create all permutations of the free columns, where free column
  # refers to absence of zero restrictions
  perm_free_col <- length(free_cols) %>% 
    replicate(list(free_cols)) %>% 
    expand.grid %>% 
    filter(apply(., 1, n_distinct)==length(free_cols)) %>% 
    t

  # Replace the non-fixed columns with permutations 
  # into columns of the helper matrix 
  perm_col_mat[(1:nvar)%in%free_cols, ] <- perm_free_col
  
  # Create all allowed permutation matrices 
  perm_mat_list <- apply(X = perm_col_mat, 
                         MARGIN = 2, 
                         FUN = function(x) diag(nvar)[,x], 
                         simplify = FALSE)
  
  # Create all allowed signed permutation matrices into a list of size (2^n*free_cols!)
  lapply(sign_mat_list, function(x) lapply(perm_mat_list, function(y) x%*%y)) %>% 
    unlist(recursive=FALSE)
}


fry_pagan_mt <- function(irf_arr, rest_mat = NULL, ct_hor = dim(irf_arr)[3], qt = 0.5){

  # Check if the median target is calculated only wrt to the restricted shocks
  if(is.null(rest_mat)){
    tgt_shocks <- 1:dim(irf_arr)[2] 
  } else {
    tgt_shocks <- which(apply(rest_mat, 2, function(x) any(is.finite(x))))
  }
  # Calculate the median scaled by respective standard deviation across draws
  scaled_irf_qt <- apply(X = irf_arr[, tgt_shocks, 1:ct_hor, , drop=FALSE], 
                         MARGIN = c(1, 2, 3), 
                         FUN = function(x) quantile(x, qt, na.rm=TRUE)/sd(x, na.rm = TRUE))
  # Subtract the scaled median from every draw, square this gap, and calculate the 
  # squared sum of these gaps for every draw
  sq_sum_error <- apply(X = irf_arr[, tgt_shocks, 1:ct_hor, , drop = FALSE], 
                        MARGIN = 4, 
                        FUN = function(x) sum((x - scaled_irf_qt)^2, na.rm=TRUE))
  # Return the draw which minimizes the quadratic loss function
  irf_arr[, , , which.min(sq_sum_error)]
}

