pkgs = c("tidyverse")
select <- dplyr::select
void = lapply(pkgs, library, character.only = TRUE)
params <- list(PATH = "local_data/jobid_",
               JOBID = "20230126")


vec_files = list.files(paste0(params$PATH, params$JOBID))
vec_files = vec_files[grepl("arrayjob", vec_files)]
vec_files %>% tail

tibble_list = vector("list", length(vec_files))

for (ix_file in seq_along(vec_files)){
  
  file_this = readRDS(paste0(params$PATH, params$JOBID, "/",
                             vec_files[ix_file]))
  
  # file_this = readRDS(paste0("local_data/p_whf/ukko_bq_jobid_", params$JOBID, "/",
  #                     vec_files[ix_file]))
  
  SCRIPT_PARAMS_this = file_this[[1]]$results_list$script_params
  
  IX_ARRAY_JOB_this = SCRIPT_PARAMS_this$IX_ARRAY_JOB
  N_MODS_this = with(SCRIPT_PARAMS_this, N_MODS_PER_CORE * N_CORES)
  
  hlp_init = 
    enframe(file_this) %>% 
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>%
    unnest_wider(results_list) %>%
    select(-params_deep_final, -value_final, -script_params) %>% 
    pivot_longer(c(initial, gaussian, ic, laplace, sgt), names_to = "est", values_to = "v") %>% 
    filter(est == "initial") %>% 
    unnest_wider(v, names_sep = "_") %>% 
    mutate(ix = 1) %>% 
    rename(theta = v_theta,
           value_laplace = v_value_laplace)
  
  hlp_gaussian_BFGS = 
    enframe(file_this) %>% 
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>%
    unnest_wider(results_list) %>%
    select(-params_deep_final, -value_final, -script_params) %>% 
    pivot_longer(c(initial, gaussian, ic, laplace, sgt), names_to = "est", values_to = "v") %>% 
    filter(est == "gaussian") %>% 
    mutate(est = "gaussian_BFGS") %>% 
    unnest_longer(v,indices_to = "ix") %>% 
    mutate(ix = ix*2) %>% 
    unnest_wider(v, names_sep = "_") %>% 
    select(-v_NM) %>% 
    unnest_wider(v_BFGS, names_sep = "_") %>% 
    select(-contains("v_BFGS_msg")) %>% 
    rename(convergence = v_BFGS_convergence,
           theta = v_BFGS_theta,
           value_laplace = v_BFGS_value_laplace,
           value_gaussian = v_BFGS_value,
           duration = v_BFGS_duration)
  
  hlp_gaussian_NM = 
    enframe(file_this) %>% 
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>%
    unnest_wider(results_list) %>%
    select(-params_deep_final, -value_final, -script_params) %>% 
    pivot_longer(c(initial, gaussian, ic, laplace, sgt), names_to = "est", values_to = "v") %>% 
    filter(est == "gaussian") %>% 
    mutate(est = "gaussian_NM") %>% 
    unnest_longer(v,indices_to = "ix") %>% 
    mutate(ix = ix*2+1) %>% 
    unnest_wider(v, names_sep = "_") %>% 
    select(-v_BFGS) %>% 
    unnest_wider(v_NM, names_sep = "_") %>% 
    rename(convergence = v_NM_convergence,
           theta = v_NM_theta,
           value_laplace = v_NM_value_laplace,
           value_gaussian = v_NM_value,
           duration = v_NM_duration)
  
  hlp_ic = 
    enframe(file_this) %>% 
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>%
    unnest_wider(results_list) %>%
    select(-params_deep_final, -value_final, -script_params) %>% 
    pivot_longer(c(initial, gaussian, ic, laplace, sgt), names_to = "est", values_to = "v") %>% 
    filter(est == "ic") %>% 
    unnest_wider(v, names_sep = "_") %>% 
    mutate(ix = 2 + 2*length(file_this[[1]]$results_list$gaussian)) %>% 
    rename(theta = v_theta,
           value_laplace = v_value_laplace)
  
  hlp_laplace_BFGS = 
    enframe(file_this) %>% 
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>%
    unnest_wider(results_list) %>%
    select(-params_deep_final, -value_final, -script_params) %>% 
    pivot_longer(c(initial, gaussian, ic, laplace, sgt), names_to = "est", values_to = "v") %>% 
    filter(est == "laplace") %>% 
    mutate(est = "laplace_BFGS") %>% 
    unnest_longer(v,indices_to = "ix") %>% 
    mutate(ix = 1 + ix*2+2*length(file_this[[1]]$results_list$gaussian)) %>% 
    unnest_wider(v, names_sep = "_") %>% 
    select(-v_NM) %>% 
    unnest_wider(v_BFGS, names_sep = "_") %>% 
    select(-contains("v_BFGS_msg")) %>% 
    rename(convergence = v_BFGS_convergence,
           theta = v_BFGS_theta,
           value_laplace = v_BFGS_value_laplace,
           duration = v_BFGS_duration) %>% 
    select(-v_BFGS_value)
  
  hlp_laplace_NM = 
    enframe(file_this) %>% 
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>%
    unnest_wider(results_list) %>%
    select(-params_deep_final, -value_final, -script_params) %>% 
    pivot_longer(c(initial, gaussian, ic, laplace, sgt), names_to = "est", values_to = "v") %>% 
    filter(est == "laplace") %>% 
    mutate(est = "laplace_NM") %>% 
    unnest_longer(v,indices_to = "ix") %>% 
    mutate(ix = ix*2 + 2 + 2*length(file_this[[1]]$results_list$gaussian)) %>% 
    unnest_wider(v, names_sep = "_") %>% 
    select(-v_BFGS) %>% 
    unnest_wider(v_NM, names_sep = "_") %>% 
    rename(convergence = v_NM_convergence,
           theta = v_NM_theta,
           value_laplace = v_NM_value_laplace,
           duration = v_NM_duration) %>% 
    select(-v_NM_value)
  
  hlp_sgt_BFGS =
    enframe(file_this) %>% 
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>%
    unnest_wider(results_list) %>%
    select(-params_deep_final, -value_final, -script_params) %>% 
    pivot_longer(c(initial, gaussian, ic, laplace, sgt), names_to = "est", values_to = "v") %>% 
    filter(est == "sgt") %>% 
    mutate(est = "sgt_BFGS") %>% 
    unnest_longer(v,indices_to = "ix") %>% 
    mutate(ix = 1+ix*2+2*(length(file_this[[1]]$results_list$gaussian) + length(file_this[[1]]$results_list$laplace))) %>% 
    unnest_wider(v, names_sep = "_") %>% 
    select(-v_NM) %>% 
    mutate(v_BFGS = map(v_BFGS, ~ .x[!names(.x)%in%"msg"])) %>% 
    unnest_wider(v_BFGS, names_sep = "_") %>% 
    rename(convergence = v_BFGS_convergence,
           theta = v_BFGS_theta,
           value_laplace = v_BFGS_value_laplace,
           duration = v_BFGS_duration,
           value_sgt = v_BFGS_value)
  
  hlp_sgt_NM =
    enframe(file_this) %>%
    rename(nr = name) %>% 
    mutate(nr = nr + (IX_ARRAY_JOB_this-1)*N_MODS_this) %>% 
    unnest_wider(value) %>%
    unnest_wider(results_list) %>%
    select(-params_deep_final, -value_final, -script_params) %>% 
    pivot_longer(c(initial, gaussian, ic, laplace, sgt), names_to = "est", values_to = "v") %>% 
    filter(est == "sgt") %>% 
    mutate(est = "sgt_NM") %>% 
    unnest_longer(v,indices_to = "ix") %>% 
    mutate(ix = 2+ix*2+2*(length(file_this[[1]]$results_list$gaussian) + length(file_this[[1]]$results_list$laplace))) %>% 
    unnest_wider(v, names_sep = "_") %>% 
    select(-v_BFGS) %>% 
    mutate(v_NM = map(v_NM, ~ .x[!names(.x)%in%"msg"])) %>% 
    unnest_wider(v_NM, names_sep = "_") %>% 
    rename(convergence = v_NM_convergence,
           theta = v_NM_theta,
           value_laplace = v_NM_value_laplace,
           duration = v_NM_duration,
           value_sgt = v_NM_value) 
  
  tibble_list[[ix_file]] = 
    bind_rows(hlp_init,
              hlp_gaussian_BFGS,
              hlp_gaussian_NM,
              hlp_ic,
              hlp_laplace_BFGS,
              hlp_laplace_NM,
              hlp_sgt_BFGS,
              hlp_sgt_NM) %>% 
    group_by(nr) %>% 
    arrange(ix, .by_group = TRUE)
  
}

tt_full = reduce(tibble_list, bind_rows)

tt_full %>% 
  pull(convergence) %>% table()

tt_full %>% 
  filter(convergence == 10) %>% 
  ungroup() %>% 
  unnest_wider(input_integerparams)

tt_full %>% 
  filter(convergence == 10) %>% 
  slice(1, .preserve = 4) %>% 
  ungroup() %>% 
  unnest_wider(input_integerparams) %>% 
  arrange(!desc(p), desc(q))

tt_full %>% 
  group_by(nr, input_integerparams) %>% 
  filter(convergence == 1) %>% 
  mutate(count_NM = grepl("_NM", est)) %>% 
  ungroup() %>% 
  summarise(mean(count_NM))

tt_full %>% 
  group_by(nr, input_integerparams) %>% 
  filter(convergence == 1) %>% 
  filter(grepl("_BFGS", est)) %>% 
  unnest_wider(input_integerparams) %>% 
  mutate(p_plus_q = p + q) %>%
  arrange(desc(p_plus_q)) %>% 
  pull(ix) %>% table()

tt_full %>% 
  ungroup() %>% 
  filter(convergence == 2) %>% 
  mutate(count_gaussian = grepl("gaussian", est),
         count_laplace = grepl("laplace", est),
         count_sgt = grepl("sgt", est)) %>% 
  summarise(across(contains("count"), ~sum(.x)))
