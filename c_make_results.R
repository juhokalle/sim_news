# SVARMA ####
tt = tt %>%
  # template
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  # generate initial values and likelihood functions (we can use the same template for initial values and likelihood fct bc both have no parameters for density)
  mutate(theta_init = map2(tmpl, data_list, ~get_init_armamod_whf_random(.y, .x))) 

# Parallel setup ####
tt_optim_parallel = tt %>% 
  select(theta_init, tmpl, data_list)

params_parallel = lapply(1:nrow(tt_optim_parallel),
                         function(i) t(tt_optim_parallel)[,i])

cl = makeCluster(params$N_CORES, type = "FORK")
mods_parallel_list <- clusterApply(cl, params_parallel, fun = hlp_parallel)
stopCluster(cl)

top_mod = 
  enframe(mods_parallel_list) %>% 
  unnest_wider(value) %>% 
  unnest_wider(results_list) %>% 
  mutate(n_params = map_int(params_deep_final, length)) %>%
  unnest_wider(input_integerparams) %>%
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>%
  mutate(value_aic = value_final + n_params * 2/params$N_OBS) %>% 
  mutate(value_bic = value_final + n_params * log(params$N_OBS)/params$N_OBS) %>% 
  arrange(value_aic) %>% 
  slice(1)

irf0 <- irf_whf(theta = top_mod$params_deep_final[[1]],
                tmpl = top_mod$tmpl[[1]],
                n_lags = n_ahead) %>% unclass
perm_ix <- order(apply(abs(sapply(1:dim(irf0)[3], function(x) diag(irf0[,,x]))/diag(irf0[,,1])), 1, max))
sign_ix <- sapply(1:dim(irf0)[2], function(x) ifelse(abs(min(irf0[,x,]))>abs(max(irf0[,x,])), -1, 1))
irf0 <- permute_chgsign(irf0, perm = perm_ix, sign = sign_ix[perm_ix])
irf_svarma[,,,mc_ix,jj] <- (pseries(irf0, n_ahead)%r%rotmat(atan(-irf0[1,2,1]/irf0[1,1,1]))) %>% unclass

mc_ix <- mc_ix + 1
#if(mc_ix%%10==0){cat(paste("MC round", mc_ix), "\n"); Sys.time()-stt}
print(Sys.time()-stt)