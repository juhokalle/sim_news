# ----------------------------------------------------------------------- #
# This script creates model output: diagnostics, model selection, IRFs... #
# ----------------------------------------------------------------------- #

pkgs = c("tidyverse", "svarmawhf")
void = lapply(pkgs, library, character.only = TRUE)
select <- dplyr::select
params <- list(PATH = "local_data/jobid_",
               JOBID = "20230303")
source("list_of_functions.R")

# sftp::sftp_connect(server = "turso.cs.helsinki.fi",
#                    folder = "/proj/juhokois/sim_news/local_data/",
#                    username = "juhokois",
#                    password = "***") -> scnx
# sftp::sftp_download(file = "jobid_20221215.zip",
#                     tofolder = "/local_data/",
#                     sftp_connection = scnx)

n_ahead <- 8
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
    mutate(readRDS(paste0(params$PATH, params$JOBID, "/total_data_sim.rds")) %>% slice(nr)) %>% 
    mutate(punish_aic = n_params * 2/nobs) %>% 
    mutate(punish_bic = n_params * log(nobs)/nobs) %>% 
    mutate(value_aic = value_final + punish_aic) %>% 
    mutate(value_bic = value_final + punish_bic) %>% 
    mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% # CHECK THIS, should have shock_distr="tdist"
    mutate(res = pmap(., pmap_get_residuals_once)) %>% 
    mutate(B_mat = map2(params_deep_final, tmpl, 
                        ~fill_tmpl_whf_rev(theta = .x, 
                                           tmpl = .y)$B)) %>% 
    mutate(shocks = map2(res, B_mat, ~ solve(.y, t(.x)) %>% t())) %>%
    select(nr, p, q, kappa, k, n_st, n_unst, value_final, value_aic, value_bic, nobs, beta, rho, nu, mc_ix, B_mat, shocks, params_deep_final) %>% 
    bind_rows(tt_full)
  #mutate(cov_shocks = map(shocks, function(x){y = abs(cov(x) - diag(DIM_OUT)); names(y) = paste0("cov_el_", letters[1:(DIM_OUT^2)]); y})) %>% 
  #unnest_wider(cov_shocks) %>% 
  #mutate(cov_el_sum = rowSums(across(contains("cov_el")))) # %>% select(-tmpl, -starts_with("punish"), -res, -B_mat)
}

tt_opt <- tt_full %>% 
  filter(value_final!=1e25) %>% 
  group_by(mc_ix, n_unst, beta, nu) %>% 
  slice_min(value_aic) %>% 
  ungroup() %>% 
  mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>% 
  #mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~irf_whf(.x, .y, n_ahead)))
  #mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~irf_unique(.x, .y, lag.max=n_ahead)))
  mutate(b0 = map2(.x = params_deep_final, .y = tmpl,
                   ~ armamod_whf(.x, .y) %>%
                     .$polm_ma %>%
                     unclass %>%
                     .[,,1]
                   )
         ) %>%
  mutate(b0 = map2(.x = b0, .y = B_mat, ~ .x%*%.y)) %>%
  mutate(rmat = map(.x = b0, ~ choose_perm_sign(cand_mat = .x, type = "dg_abs")[[2]])) %>%
  mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~ irf_whf(theta = .x, tmpl = .y, n_lags = n_ahead))) %>% 
  mutate(irf = map2(.x = irf, .y = rmat, ~ .x%r%.y))
  # mutate(irf = map(.x = irf, ~ id_news_shox(unclass(.x))))

mc_n <- unique(tt_opt$mc_ix)
prms <- expand.grid(beta = unique(tt_opt$beta), nu = unique(tt_opt$nu))
# dim_out x dim_out x n_ahead x {(beta_1,nu_1), ...,(beta_2,nu_2)} x {mc_1, mc_2, ..., mc_n} x {VAR, VARMA}
irf_arr <- array(NA, c(DIM_OUT, DIM_OUT, n_ahead+1, nrow(prms), length(mc_n), 2))

margs_plus_font <- theme(plot.title = element_text(size = 8, hjust = 0.5),
                         plot.margin=grid::unit(c(1, 1, 0, 0), "mm"),
                         axis.text.x = element_text(size = 6),
                         axis.text.y = element_text(size = 6))

plot1 <- map2(.x = sim_comparison(irf_arr, c(2, 2), c(.14, .86), prms),
              .y = list(c(-0.2, 1.5), c(-0.25, 2.2),
                        c(-0.2, 1.5), c(-0.25, 2.2)),
                 ~ .x + margs_plus_font + ylim(.y)
                )

plot2 <- map2(.x = sim_comparison(irf_arr, c(1, 2), c(.14,.86), prms),
              .y = list(c(-0.3,1.2), c(-0.25,1.2),
                        c(-0.3,1.2), c(-0.25,1.2)),
                 ~ .x + margs_plus_font + ylim(.y)
              )

ggsave(filename = "./paper_output/sim_plt.pdf",
       gridExtra::marrangeGrob(plot1, nrow = 2, ncol = 2, top = NULL),
       width = 7.5,
       height = 5)
ggsave(filename = "./paper_output/sim_plt_a.pdf",
       plot1[[1]],
       width = 7.5,
       height = 5)
ggsave(filename = "./paper_output/sim_plt_b.pdf",
       plot1[[3]],
       width = 7.5,
       height = 5)
ggsave(filename = "./paper_output/sim_plt_d.pdf",
       plot1[[4]],
       width = 7.5,
       height = 5)
ggsave(filename = "./paper_output/appendix1.pdf",
       gridExtra::marrangeGrob(plot2, nrow = 2, ncol = 2, top = NULL),
       width = 7.5,
       height = 5)
