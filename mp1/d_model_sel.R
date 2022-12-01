# ----------------------------------------------------------------------- #
# This script creates model output: diagnostics, model selection, IRFs... #
# ----------------------------------------------------------------------- #

pkgs = c("tidyverse", "svarmawhf")
void = lapply(pkgs, library, character.only = TRUE)
select <- dplyr::select
params <- list(PATH = "local_data/jobid_",
               JOBID = "20221118")

# functions for the analysis
norm_irf <- function(irf_arr, norm_int, norm_pos=dim(irf_arr)[1], norm_scm = c("surp", "news"))
{
  
  if(norm_scm=="surp"){
    irf_arr <- irf_arr[,1,,drop=FALSE]/irf_arr[norm_pos,1,1]*norm_int
  } else if(norm_scm=="news"){
    irf_arr <- irf_arr[,1,,drop=FALSE]/max(irf_arr[norm_pos,1,])*norm_int
  }
  irf_arr
}

plot_irf <- function(irf_arr, var_name, shock_name){
  n_var <- dim(irf_arr)[1]
  n_shock <- dim(irf_arr)[2]
  n_ahead <- dim(irf_arr)[3]-1
  tibble(irf = c(irf_arr),
         variable = var_name %>% 
           factor(levels = var_name) %>% 
           rep((n_ahead+1)*n_shock),
         shock = shock_name %>% 
           factor(levels=shock_name) %>% 
           rep(each=n_var) %>% 
           rep(n_ahead+1),
         months = rep(0:n_ahead, each = n_var*n_shock)) %>% 
    ggplot(aes(x=months, y = irf)) +
    geom_line() +
    geom_hline(yintercept = 0, size = 0.15) +
    facet_grid(variable ~ shock, scales = "free_y") +
    facet_rep_grid(variable ~ shock, scales = "free_y", repeat.tick.labels = 'left') +
    scale_x_continuous(breaks = seq(12, n_ahead, by = 12), expand = c(0,0)) +
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

get_llf <- function(p, q, kappa, k, dtype)
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
  ll_whf_factory(data_wide = t(data_i), tmpl = tmpl_i, shock_distr = "sgt")
  
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

optim_zr <- function(input_mat, row_ix)
{
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

get_rest_irf <- function(tbl_slice, rest_ix)
{
  nobs <- tbl_slice$nobs
  sd_mat <- readRDS("local_data/svarma_data_list.rds") %>% 
    filter(mp_type%in%tbl_slice$mp_type) %>%
    pull(std_dev) %>% .[[1]] %>%  diag
  
  irf_out <- sd_mat %r% tbl_slice$irf[[1]]
  
  rmat <- optim_zr(irf_out[,,1], rest_ix)
  llf0 <- do.call(get_llf, tbl_slice %>% select(p, q, kappa, k, mp_type) %>% rename(dtype = mp_type))
  params0 <- tbl_slice$params_deep_final[[1]]
  
  params0[tbl_slice$params_deep_final[[1]] %in% tbl_slice$B_mat[[1]]] <- c(tbl0$B_mat[[1]]%*%rmat)
  pval <- 1-pchisq(2*nobs*(llf0(params0)-tbl0$value_final), length(rest_ix))
  
  irf_out <- irf_whf(params0, tbl_slice$tmpl[[1]], 48)
  
  return(list(irf = irf_out, pval = pval, rmat = rmat))
}

pmap_tmpl_whf_rev = function(dim_out = DIM_OUT, p, q, kappa, k, shock_distr = "sgt", ...)
{
  tmpl_whf_rev(dim_out = DIM_OUT, ARorder = p, MAorder = q, kappa = kappa, k = k, shock_distr = shock_distr)
}

pmap_get_residuals_once = function(params_deep_final, tmpl, data_list, ...)
{
  get_residuals_once(params_deep = params_deep_final, tmpl = tmpl, data_long = data_list)
}

get_fevd <- function (irf_arr)
{
  nvar <- dim(irf_arr)[1]
  n.ahead <- dim(irf_arr)[3]
  fe <- list()
  for (i in 1:nvar) {
    fe[[i]] <- as.data.frame(t(irf_arr[i, , ]))
  }
  fe2 <- fe
  for (i in 1:length(fe)) {
    for (j in 1:n.ahead) {
      fe2[[i]][j, ] <- (colSums(fe[[i]][j:1, ]^2)/sum(fe[[i]][j:1, ]^2)) * 100
    }
  }
  return(fe2)
}

# sftp::sftp_connect(server = "turso.cs.helsinki.fi",
#                    folder = "/proj/juhokois/sim_news/local_data/",
#                    username = "juhokois",
#                    password = "securitasLOKKIv4lvoo") -> scnx
# sftp::sftp_download(file = "jobid_20221118.zip",
#                     tofolder = "/local_data/",
#                     sftp_connection = scnx)
# sftp::sftp_download(file = "total_data.rds",
#                     tofolder = "/local_data/",
#                     sftp_connection = scnx)
vec_files = list.files(paste0(params$PATH, params$JOBID))
vec_files = vec_files[grepl("arrayjob", vec_files)]
SCRIPT_PARAMS = readRDS(paste0(params$PATH, params$JOBID, "/", vec_files[1]))[[1]]$results_list$script_params
DIM_OUT = SCRIPT_PARAMS$DIM_OUT
  
tt_full <- tibble()
for (ix_file in seq_along(vec_files))
{
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
    mutate(nobs = map_dbl(data_list, ~nrow(.x))) %>% 
    mutate(punish_aic = n_params * 2/nobs) %>% 
    mutate(punish_bic = n_params * log(nobs)/nobs) %>% 
    mutate(value_aic = value_final + punish_aic) %>% 
    mutate(value_bic = value_final + punish_bic) %>% 
    mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
    mutate(res = pmap(., pmap_get_residuals_once)) %>% 
    mutate(B_mat = map2(params_deep_final, tmpl, 
                        ~fill_tmpl_whf_rev(theta = .x, 
                                           tmpl = .y)$B)) %>% 
    mutate(shocks = map2(res, B_mat, ~ solve(.y, t(.x)) %>% t())) %>%
    select(nr, p, q, kappa, k, n_st, n_unst, value_final, value_aic, value_bic, nobs, mp_type, B_mat, shocks, params_deep_final) %>% 
    bind_rows(tt_full)
    #mutate(cov_shocks = map(shocks, function(x){y = abs(cov(x) - diag(DIM_OUT)); names(y) = paste0("cov_el_", letters[1:(DIM_OUT^2)]); y})) %>% 
    #unnest_wider(cov_shocks) %>% 
    #mutate(cov_el_sum = rowSums(across(contains("cov_el")))) # %>% select(-tmpl, -starts_with("punish"), -res, -B_mat)
}
#saveRDS(tt_full, "local_data/results_total_20221122.rds")
#tt_full <- readRDS("local_data/results_199401_201312.rds")
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

tt %>% pull(norm_indep_flag) %>% table

tt %>%
  #mutate(n_params = map_int(params_deep_final, length)) %>% 
  filter(norm_indep_flag==0) %>%
  filter(nobs==240) %>% 
  group_by(mp_type, p_plus_q) %>%
  summarise(n=n()) %>% 
  pivot_wider(names_from = mp_type, values_from = n) %>% 
  arrange(p_plus_q)

tt %>% filter(mp_type=="BS22a") %>% 
  filter(norm_indep_flag==0) %>% 
  arrange(value_aic)

tbl0 <- tt %>% filter(nr==1120) %>% 
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~ irf_whf(.x, .y, n_lags = 48)))

sd_mat <- readRDS("local_data/svarma_data_list.rds") %>% 
  filter("BS22a"%in%tbl0$mp_type) %>%
  pull(std_dev) %>% .[[1]] %>%  diag

irf_out <- sd_mat %r% tbl0$irf[[1]]

irf_bs <- map2(.x = list(irf_out[,1,,drop=FALSE],
                         irf_out[,5,,drop=FALSE]),
                .y = c("news", "surp"), 
                ~norm_irf(irf_arr = .x, norm_scm = .y, norm_pos = 4, norm_int = 0.25)) %>%
  reduce(abind::abind, along=2)

irf_bs[c(1:2,4),,] %>%
  plot_irf(var_name = c("LogIP", "LogCPI", "FFR"),
           shock_name = c("News", "Surprise")) -> p1

ggsave(filename = "paper_output/IRF1.pdf", plot = p1)

tt %>% filter(mp_type=="BRW21") %>% 
  filter(norm_indep_flag==0) %>% 
  arrange(value_aic)

tbl0 <- tt %>% filter(nr==2248) %>% 
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~ irf_whf(.x, .y, n_lags = 48)))

sd_mat <- readRDS("local_data/svarma_data_list.rds") %>% 
  filter("BRW21"%in%tbl0$mp_type) %>%
  pull(std_dev) %>% .[[1]] %>%  diag

irf_out <- sd_mat %r% tbl0$irf[[1]]

irf_brw <- map2(.x = list(irf_out[,1,,drop=FALSE],
                         irf_out[,3,,drop=FALSE]),
               .y = c("news", "surp"), 
               ~norm_irf(irf_arr = .x, norm_scm = .y, norm_pos = 4, norm_int = 0.25)) %>%
  reduce(abind::abind, along=2)

irf_bs[c(1:2,4),,] %>%
  plot_irf(var_name = c("LogIP", "LogCPI", "FFR"),
           shock_name = c("News", "Surprise")) -> p2

ggsave(filename = "paper_output/IRF2.pdf", plot = p2)
