# preliminaries ####
main_path <- "local_data/jobid_20220912/"
vec_files <- list.files(main_path)

# modify file names such that they are ordered correctly 1,2,...
for(ix_file in seq_along(vec_files)){
  file_sub <- sub(".*_", "", sub("\\..*", "", vec_files[ix_file]))
  if((nchar(file_sub)==1 && length(vec_files) < 100) || (nchar(file_sub)==2 && length(vec_files) >= 100)){
    file_sub <- paste0("0", file_sub)
  } else if(nchar(file_sub)==1 && length(vec_files) > 100){
    file_sub <- paste0("00", file_sub)
  }
  file.rename(paste0(main_path, vec_files[ix_file]),
              paste0(main_path, "arrayjob_", file_sub, ".rds"))
}

tibble_list = vector("list", length(vec_files))

vec_files <- list.files(main_path)
SCRIPT_PARAMS = readRDS(paste0(main_path, vec_files[1]))[[1]]$results_list$script_params
DIM_OUT = SCRIPT_PARAMS$DIM_OUT

# helper functions ####
permute_chgsign = function(irf_array, 
                           perm = rep(1, dim(irf_array)[1]), 
                           sign = rep(1, dim(irf_array)[1])){
  
  dim_out = dim(irf_array)[1]
  
  perm_mat = diag(dim_out)
  perm_mat = perm_mat[,perm]
  
  sign = diag(sign)
  
  ll = map(1:dim(irf_array)[3], ~ irf_array[,,.x] %*% perm_mat %*% sign)
  
  irf_array = array(0, dim = c(dim_out, dim_out, dim(irf_array)[3]))
  for (ix in 1:dim(irf_array)[3]){
    irf_array[,,ix] = ll[[ix]] 
  }
  
  return(irf_array)
}
rotmat <- function(x, n){
  
  rot_mat <- matrix(c(cos(x), -sin(x), sin(x), cos(x)), 2,2)
  rot_ix <- combn(n,2)
  final_mat <- diag(n)
  
  for(jj in 1:ncol(rot_ix)){
    temp_mat <- diag(n)
    temp_mat[rot_ix[,jj], rot_ix[,jj]] <- rot_mat
    final_mat <- final_mat %*% temp_mat
  }
  
  final_mat
}


ff <- function(x){
  abs(input_mat[nrow(input_mat),] %*% rotmat(x, ncol(input_mat))[, ncol(input_mat)])
}

for (ix_file in seq_along(vec_files)){
  file_this = readRDS(paste0(main_path, vec_files[ix_file]))
  
  SCRIPT_PARAMS_this = file_this[[1]]$results_list$script_params
  N_NOISE_PARAMS_SGT = DIM_OUT^2 + DIM_OUT*3
  
  IX_ARRAY_JOB_this = SCRIPT_PARAMS_this$IX_ARRAY_JOB
  N_MODS_this = with(SCRIPT_PARAMS_this, N_MODS_PER_CORE * N_CORES)
  
  tibble_list[[ix_file]] =  
    enframe(file_this) %>% 
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>% 
    unnest_wider(results_list) %>% 
    select(nr, params_deep_final, value_final, input_integerparams) %>% 
    mutate(n_params = map_int(params_deep_final, length)) %>% 
    mutate(n_params_sys = n_params - N_NOISE_PARAMS_SGT) %>% # SGT has (too) many noise parameters, so I want to check what happens when one uses only the number of system parameters as punishment
    unnest_wider(input_integerparams) %>% 
    mutate(punish_aic = n_params * 2/SCRIPT_PARAMS_this$N_OBS) %>% 
    mutate(punish_bic = n_params * log(SCRIPT_PARAMS_this$N_OBS)/SCRIPT_PARAMS_this$N_OBS) %>% 
    mutate(value_aic = value_final + punish_aic) %>% 
    mutate(value_bic = value_final + punish_bic) %>% 
    select(-value_final, -starts_with("punish"))
}

tt_full = reduce(tibble_list, bind_rows)

tt_full <- bind_cols(tt_full, tt %>% select(setdiff(names(tt), names(tt_full))))
tt_full <- tt_full %>% 
  group_by(prm_ix, mc_ix) %>% 
  slice_min(order_by = value_bic) %>% 
  ungroup %>%
  mutate(irf = map2(params_deep_final, tmpl, ~irf_whf(.x, .y, n_lags = 12)))

n_ahead <- 12
irf_all <- array(0, c(2, 2, n_ahead+1, max(tt_full$mc_ix), max(tt_full$prm_ix)))
# 
for(para_ix in 1:4){
  
  for(mc_i in 1:101){
    irf0 <- tt_full %>% 
      filter(prm_ix==para_ix) %>% 
      select(irf) %>%
      slice(mc_i) %>% 
      deframe %>% .[[1]] %>% 
      unclass
    perm_ix <- order(apply(abs(sapply(1:dim(irf0)[3], function(x) diag(irf0[,,x]))/diag(irf0[,,1])), 1, max))
    sign_ix <- sapply(1:dim(irf0)[2], function(x) ifelse(abs(min(irf0[,x,]))>abs(max(irf0[,x,])), -1, 1))
    irf0 <- permute_chgsign(irf0, perm = perm_ix, sign = sign_ix[perm_ix])
    irf_all[,,,mc_i,para_ix] <- (pseries(irf0, n_ahead)%r%rotmat(atan(-irf0[1,2,1]/irf0[1,1,1]),2)) %>% unclass
  }
}

apply(irf_all[,,,,1], c(1,2,3), median) %>% pseries(12) %>% plot
