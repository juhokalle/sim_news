pkgs = c("tidyverse", "svarmawhf")
select <- dplyr::select
void = lapply(pkgs, library, character.only = TRUE)
params <- list(PATH = "local_data/jobid_",
               JOBID = "20221118")
get_llf <- function(p, q, kappa, k, dtype){
  total_data <- readRDS("local_data/total_data.rds")
  data_i <- total_data %>% 
    filter(type==dtype) %>% 
    pull(data_list) %>% .[[1]]
  dim_out <- dim(data_i)[2]
  tmpl_i <- tmpl_whf_rev(dim_out = dim_out,
                         ARorder = p,
                         MAorder = q,
                         kappa = kappa,
                         k = k,
                         shock_distr = "sgt")
  ll_whf_factory(data_wide = t(data_i), tmpl = tmpl_i, shock_distr = "sgt")
  
}


sftp::sftp_connect(server = "turso.cs.helsinki.fi",
                   folder = "/proj/juhokois/sim_news/local_data/",
                   username = "juhokois",
                   password = "***") -> scnx
file_ix <- 1
file_dl <- NULL
while(!inherits(file_dl, 'try-error')){

  zz <- if(file_ix<10) "00" else if(file_ix < 100) "0"
  sftp::sftp_download(paste0("jobid_20221118/arrayjob_", zz, file_ix, ".rds"),
                      tofolder = "/local_data/",
                      sftp_connection = scnx) %>%
    try() %>% suppressWarnings() -> file_dl
  file_ix <- file_ix + 1
}
sftp::sftp_download(file = "total_data.rds",
                    tofolder = "/local_data/",
                    sftp_connection = scnx)
vec_files = list.files(paste0(params$PATH, params$JOBID))
vec_files = vec_files[grepl("arrayjob", vec_files)]
SCRIPT_PARAMS = readRDS(paste0(params$PATH, params$JOBID, "/", vec_files[1]))[[1]]$results_list$script_params
DIM_OUT = SCRIPT_PARAMS$DIM_OUT
  
pmap_tmpl_whf_rev = function(dim_out = DIM_OUT, p, q, kappa, k, shock_distr = "sgt", ...){
  tmpl_whf_rev(dim_out = DIM_OUT, ARorder = p, MAorder = q, kappa = kappa, k = k, shock_distr = shock_distr)
}

pmap_get_residuals_once = function(params_deep_final, tmpl, data_list, ...){
  get_residuals_once(params_deep = params_deep_final, tmpl = tmpl, data_long = data_list)
}

tt_full <- tibble()
for (ix_file in seq_along(vec_files)){
  file_this = readRDS(paste0(params$PATH, params$JOBID, "/",
                             vec_files[ix_file]))
  
  SCRIPT_PARAMS_this = file_this[[1]]$results_list$script_params
  
  IX_ARRAY_JOB_this = SCRIPT_PARAMS_this$IX_ARRAY_JOB
  N_MODS_this = with(SCRIPT_PARAMS_this, N_MODS_PER_CORE * N_CORES)
  
  tt_full =  
    enframe(file_this) %>% 
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>% 
    unnest_wider(results_list) %>% 
    select(nr, params_deep_final, value_final, input_integerparams) %>% 
    mutate(n_params = map_int(params_deep_final, length)) %>% 
    unnest_wider(input_integerparams) %>% 
    mutate(readRDS("local_data/total_data.rds") %>% slice(nr)) %>% 
    mutate(punish_aic = map2_dbl(.x = data_list, .y = n_params, ~ .y * 2/nrow(.x))) %>% 
    mutate(punish_bic = map2_dbl(.x = data_list, .y = n_params, ~ .y * log(nrow(.x))/nrow(.x))) %>% 
    mutate(value_aic = value_final + punish_aic) %>% 
    mutate(value_bic = value_final + punish_bic) %>% 
    mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
    mutate(res = pmap(., pmap_get_residuals_once)) %>% 
    mutate(B_mat = map2(params_deep_final, tmpl, 
                        ~fill_tmpl_whf_rev(theta = .x, 
                                           tmpl = .y)$B)) %>% 
    mutate(shocks = map2(res, B_mat, ~ solve(.y, t(.x)) %>% t())) %>%
    select(nr, value_aic, value_bic, value_final, p, q, B_mat, shocks, n_st, n_unst, type, kappa, k, params_deep_final) %>% 
    bind_rows(tt_full)
    #mutate(cov_shocks = map(shocks, function(x){y = abs(cov(x) - diag(DIM_OUT)); names(y) = paste0("cov_el_", letters[1:(DIM_OUT^2)]); y})) %>% 
    #unnest_wider(cov_shocks) %>% 
    #mutate(cov_el_sum = rowSums(across(contains("cov_el")))) # %>% select(-tmpl, -starts_with("punish"), -res, -B_mat)
}

#tt_full = reduce(tibble_list, bind_rows)
tt = tt_full %>% 
  mutate(rk_aic = rank(value_aic),
         rk_bic = rank(value_bic),
         rk_mle = rank(value_final)) %>% 
  arrange(value_aic) %>% 
  mutate(p_plus_q = p+q)

THRESHOLD_SW = 0.05 

# filter good models by flag == 0
# H_0: Normality -> good models have small p-values
tt = tt %>% 
  mutate(sw = map(shocks, ~apply(.x, 2, FUN = function(x){shapiro.test(x)$p.value}))) %>% 
  mutate(sw_flag = map_int(sw, ~sum(.x > THRESHOLD_SW))) %>% # one component may be Gaussian
  mutate(sw_pval_sum = map_dbl(sw, sum)) %>% 
  unnest_wider(sw, names_sep = "_pval") %>% 
  arrange(desc(sw_pval_sum))

tt %>% 
  pull(sw_flag) %>% table()

THRESHOLD_JB = 0.05 

# filter good moodels by flag == 0
# H_0: Normality -> good models have small p-values
tt = tt %>% 
  mutate(jb = map(shocks, ~apply(.x, 2, FUN = function(x){tsoutliers::JarqueBera.test(x)[[1]]$p.value}))) %>% 
  mutate(jb_flag = map_int(jb, ~sum(.x > THRESHOLD_JB))) %>% # one component may be Gaussian
  mutate(jb_pval_sum = map_dbl(jb, sum)) %>% 
  unnest_wider(jb, names_sep = "_pval") %>% 
  arrange(desc(jb_pval_sum))

tt %>% 
  pull(jb_flag) %>% table()

tt = tt %>% mutate(normality_flag = sw_flag + jb_flag)
tt %>% pull(normality_flag) %>% table()

THRESHOLD_LB = 0.05

# filter good moodels by flag == 0
# H_0: No autocorrelation of (transformation of) residuals
# -> good models have high p-values
tt = tt %>% 
  mutate(lb = map2(.x = shocks, .y = p_plus_q, ~ apply(.x, 2, FUN = function(x){ Box.test(x, lag = 24, type = "Ljung-Box")$p.value }))) %>% 
  mutate(lb_flag = map_lgl(lb, ~ any(.x < THRESHOLD_LB))) %>% 
  mutate(lb_pval_sum = map_dbl(lb, sum)) %>%
  unnest_wider(lb, names_sep = "_pval") %>% 
  
  mutate(lb_abs = map2(.x = shocks, .y = p_plus_q, ~ apply(.x, 2, FUN = function(x){ Box.test(abs(x), lag = 24, type = "Ljung-Box")$p.value }))) %>% 
  mutate(lb_abs_flag = map_lgl(lb_abs, ~any(.x < THRESHOLD_LB))) %>% 
  mutate(lb_abs_pval_sum = map_dbl(lb_abs, sum)) %>%
  unnest_wider(lb_abs, names_sep = "_pval") %>% 
  
  mutate(lb_sq = map2(.x = shocks, .y = p_plus_q, ~ apply(.x, 2, FUN = function(x){ Box.test(x^2, lag = 24, type = "Ljung-Box")$p.value }))) %>% 
  mutate(lb_sq_flag = map_lgl(lb_sq, ~any(.x < THRESHOLD_LB))) %>% 
  mutate(lb_sq_pval_sum = map_dbl(lb_sq, sum)) %>%
  unnest_wider(lb_sq, names_sep = "_pval") %>% 
  
  mutate(lb_all_pval_sum = lb_pval_sum + lb_abs_pval_sum + lb_sq_pval_sum) %>% 
  arrange(lb_all_pval_sum)

tt = tt %>% mutate(indep_flag = lb_flag + lb_abs_flag + lb_sq_flag)
tt %>% pull(indep_flag) %>% table()

tt <- tt %>% mutate(norm_indep_flag = indep_flag+normality_flag)

tt %>%
  #mutate(n_params = map_int(params_deep_final, length)) %>% 
  filter(norm_indep_flag==0) %>% 
  group_by(type, p_plus_q) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = type, values_from = n)

tt %>% filter(type=="MAR21a") %>% filter(norm_indep_flag==0) %>% arrange(value_bic)

rotmat <- function(x, n){
  
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

ff <- function(x, zero_ix, input_mat){
  
  row_ix <- if(is.null(dim(zero_ix))) zero_ix[1] else zero_ix[,1]
  col_ix <- if(is.null(dim(zero_ix))) zero_ix[2] else zero_ix[,2]
  sum(sqrt(diag(input_mat[row_ix,] %*% rotmat(x, ncol(input_mat))[, col_ix])^2))
}

optim_zr <- function(input_mat, row_ix){
  
  n_var <- dim(input_mat)[1]
  opt_ix <- rep(0, n_var)
  rotmat_ix <- array(NA, dim = c(n_var, n_var, n_var))
  n_rest <- length(row_ix)
  
  for(jj in 1:n_var){
    # UNCOSTRAINED OPTIM: BFGS
    zero_ix = if(n_rest==1) c(row_ix, jj) else matrix(c(row_ix, rep(jj, n_rest)), n_rest, n_rest)
    optim(par = rep(0, (n_var^2-n_var)/2),
          fn = ff, 
          zero_ix = zero_ix,
          input_mat = input_mat,
          method = "BFGS") -> opt_obj
    
    rotmat_ix[,,jj] <- rotmat(opt_obj$par, n_var)
    opt_ix[jj] <- norm(input_mat - input_mat%*%rotmat_ix[,,jj], "F")
  }
  rotmat_ix[,,which.min(opt_ix)]
}

get_rest_irf <- function(tbl_slice, rest_ix){
  
  fred_md <- fbi::fredmd("http://files.stlouisfed.org/files/htdocs/fred-md/monthly/current.csv",
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
  ebp <- read_csv("https://www.federalreserve.gov/econres/notes/feds-notes/ebp_csv.csv")
  fred_md$EBP <- ebp$ebp[1:nrow(fred_md)]
  ds <- fred_md %>% 
    dplyr::select(date, LIP, LCPI, EBP, FEDFUNDS) %>% 
    filter(date >= ym(199401), date<ym(201401))
  
  MPR <- if(tbl_slice$type=="GK1"){
    readRDS("local_data/shock_tbl.rds") %>% filter(date %in% ds$date) %>% dplyr::select(mp1_tc)  %>% .[[1]]
  } else if(tbl_slice$type=="GK2"){
    readRDS("local_data/shock_tbl.rds") %>% filter(date %in% ds$date) %>% dplyr::select(ff4_tc) %>% .[[1]]
  } else if(tbl_slice$type=="MAR21a"){
    readRDS("local_data/shock_tbl.rds") %>% filter(date %in% ds$date) %>% dplyr::select(MM_IV1) %>% .[[1]]
  } else if(tbl_slice$type=="MAR21b"){
    readRDS("local_data/shock_tbl.rds") %>% filter(date %in% ds$date) %>% dplyr::select(MM_IV5) %>% .[[1]]
  } else if(tbl_slice$type=="Jaro22"){
    readRDS("local_data/shock_tbl.rds") %>% filter(date %in% ds$date) %>% dplyr::select(u1) %>% .[[1]]
  } else if(tbl_slice$type=="AD22"){
    readRDS("local_data/shock_tbl.rds") %>% filter(date %in% ds$date) %>% dplyr::select(Shock) %>% .[[1]]
  } else if(tbl_slice$type=="BS22a"){
    readRDS("local_data/shock_tbl.rds") %>% filter(date %in% ds$date) %>% dplyr::select(MPS) %>% .[[1]]
  } else if(tbl_slice$type=="BS22b"){
    readRDS("local_data/shock_tbl.rds") %>% filter(date %in% ds$date) %>% dplyr::select(MPS_ORTH) %>% .[[1]]
  } else if(tbl_slice$type=="BRW21"){
    readRDS("local_data/shock_tbl.rds") %>% filter(date %in% ds$date) %>% dplyr::select(BRW_monthly) %>% .[[1]]
  } else if(tbl_slice$type=="GSS22"){
    readRDS("local_data/shock_tbl.rds") %>% filter(date %in% ds$date) %>% dplyr::select(MP1) %>% .[[1]]
  } else if(tbl_slice$type=="RR04"){
    readRDS("local_data/shock_tbl.rds") %>% filter(date %in% ds$date) %>% dplyr::select(resid_full) %>% .[[1]]
  }
  
  ds$MPR <- MPR
  ds %>% mutate(MPR = cumsum(coalesce(MPR, 0)) + MPR*0) %>% 
    filter(complete.cases(.)) %>% 
    dplyr::select(-date) %>% 
    mutate_all(~ lm(.x ~ I(1:n()) + I((1:n())^2)) %>% residuals) -> ds
  nobs <- nrow(ds)
  sd_mat <- apply(ds, 2, sd) %>% diag()
  
  irf_out <- sd_mat %r% tbl_slice$irf[[1]]
  
  rmat <- optim_zr(irf_out[,,1], rest_ix)
  llf0 <- do.call(get_llf, tbl_slice %>% select(p, q, kappa, k, type) %>% rename(dtype = type))
  params0 <- tbl_slice$params_deep_final[[1]]
  
  params0[tbl_slice$params_deep_final[[1]] %in% tbl_slice$B_mat[[1]]] <- c(tbl0$B_mat[[1]]%*%rmat)
  pval <- 1-pchisq(2*nobs*(llf0(params0)-tbl0$value_final), length(rest_ix))

  irf_out <- irf_whf(params0, tbl_slice$tmpl[[1]], 48)
  
  return(list(irf = irf_out, pval = pval, rmat = rmat))
}

tbl0 <- tt %>% filter(nr==1105) %>% 
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~ irf_whf(.x, .y, n_lags = 48)))

irf0 <- get_rest_irf(tbl0, rest_ix = 6)
irf0$irf %>% plot
