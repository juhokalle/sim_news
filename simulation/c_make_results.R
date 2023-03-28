# ----------------------------------------------------------------------- #
# This script creates model output: diagnostics, model selection, IRFs... #
# ----------------------------------------------------------------------- #

source("list_of_functions.R")
pkgs = c("tidyverse", "svarmawhf")
void = lapply(pkgs, library, character.only = TRUE)
select <- dplyr::select
params <- list(PATH = "local_data/jobid_",
               JOBID = "202303270")

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

tibble_list <- list()
for (ix_file in seq_along(vec_files))
{
  file_this = readRDS(paste0(params$PATH, params$JOBID, "/",
                             vec_files[ix_file]))
  
  SCRIPT_PARAMS_this = file_this[[1]]$results_list$script_params
  
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
    unnest_wider(input_integerparams) %>% 
    mutate(readRDS(paste0(params$PATH, params$JOBID, "/total_data_sim.rds")) %>% slice(nr)) %>% 
    # mutate(punish_aic = n_params * 2/nobs) %>% 
    # mutate(punish_bic = n_params * log(nobs)/nobs) %>% 
    # mutate(value_aic = value_final + punish_aic) %>% 
    # mutate(value_bic = value_final + punish_bic) %>% 
    mutate(shock_distr = "tdist") %>% 
    mutate(tmpl = pmap(., pmap_tmpl_whf_rev)) %>%
    mutate(res = pmap(., pmap_get_residuals_once)) %>% 
    select(nr, p, q, kappa, k, n_st, n_unst, beta, rho, nu, mc_ix, value_final,
           res,
           params_deep_final,
           tmpl)
}

tt_opt <- reduce(tibble_list, bind_rows) %>% 
  filter(value_final!=1e25) %>%
  mutate(irf = map2(.x = params_deep_final, .y = tmpl, ~ irf_whf(.x, .y, n_lags = 8))) %>% 
  mutate(B_mat = map2(.x = params_deep_final, .y = tmpl,
                      ~fill_tmpl_whf_rev(theta = .x,
                                         tmpl = .y)$B)) %>%
  mutate(aux_ix = map2(.x = params_deep_final, .y = B_mat, ~ .x%in%c(.y))) %>%
  mutate(rmat = map(.x = irf, ~ id_news_shox(cand_mat = .x, policy_var = 1))) %>%
  mutate(rmat = map2(.x = irf, .y = rmat, ~ .y%*%optim_zr(unclass(.x)[,,1]%*%.y, c(1,2), opt_it = FALSE))) %>%
  mutate(irf = map2(.x = irf, .y = rmat, ~ .x%r%.y)) %>% 
  mutate(B_mat = map2(.x = B_mat, .y = rmat, ~ .x%*%.y)) %>%
  # Shocks, rotated now appropriately
  mutate(shocks = map2(res, B_mat, ~ solve(.y, t(.x)) %>% t())) %>%
  # # Shock covariance
  mutate(Sigma = map(.x = shocks, ~ cov(.x))) %>%
  # Normalize to one std dev shocks
  mutate(irf = map2(.x = irf, .y = Sigma, ~ .x%r%diag(diag(.y)^-.5)))
  # Finally multiply the irf by the resp. std devs to get units in the original scale
  # mutate(irf = map2(.x = std_dev, .y = irf, ~  diag(.x)%r%.y))

mc_n <- unique(tt_opt$mc_ix)
prms <- expand.grid(beta = unique(tt_opt$beta), nu = unique(tt_opt$nu))
# dim_out x dim_out x n_ahead x {(beta_1,nu_1), ...,(beta_2,nu_2)} x {mc_1, mc_2, ..., mc_n} x {VAR, VARMA}
irf_arr <- array(NA, c(DIM_OUT, DIM_OUT, n_ahead+1, nrow(prms), length(mc_n), 2))

margs_plus_font <- theme(plot.title = element_text(size = 8, hjust = 0.5),
                         plot.margin = grid::unit(c(1, 1, 0, 0), "mm"),
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
