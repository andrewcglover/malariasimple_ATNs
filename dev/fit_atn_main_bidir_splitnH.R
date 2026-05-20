# fit_atn_main_bidir_splitnH.R
# Bidirectional ATN blocking-probability driver, SPLIT-nH variant.
#
# Identical in structure to fit_atn_main_bidir.R, but uses the split-nH
# Stan files so nH_pre and nH_post are estimated separately on each side.
# Intended to be compared against the shared-nH bidir fits via LOO
# (see compare_atn_main_bidir_loo.R).
#
# Outputs:
#   atn_fit_tba_main_bidir_splitnH.rds
#   atn_fit_tra_main_bidir_splitnH.rds
#   atn_main_bidir_splitnH_pairs_<tag>.png
#   atn_main_bidir_splitnH_trace_<tag>.png
#   ppc_atn_main_bidir_splitnH_<tag>_ctl.png
#   ppc_atn_main_bidir_splitnH_<tag>_trt.png
# where <tag> in {tba, tra}.

suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(bayesplot)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
})

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)
set.seed(2025)

# ---- Toggles ----
RUN_TBA <- TRUE
RUN_TRA <- TRUE

# ---- Paths ----
here_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  fn <- sub("--file=", "", args[grep("--file=", args)])
  if (length(fn) && nzchar(fn)) {
    dirname(normalizePath(fn))
  } else if (requireNamespace("rstudioapi", quietly = TRUE) &&
             rstudioapi::isAvailable()) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    getwd()
  }
}, error = function(e) getwd())
setwd(here_dir)

data_path <- "./exp_data/mos_level_int_summary_pre_post_full.csv"

# ---- Priors (single source of truth: priors_atn_main.R) ----
source("priors_atn_main.R")
prior_data <- prior_stan_data(priors)

# ---- Data ----
mos <- read.csv(data_path)

stopifnot(
  all(c("exp_id", "pre_time", "post_time", "int", "group") %in% names(mos)),
  all(mos$group %in% c("control", "200mg")),
  all(mos$int >= 0),
  all(abs(mos$pre_time + mos$post_time) < 1e-6)
)

# Bidir fit: keep all rows.  side = 0 if exposure before infection
# (pre_time >= 0), side = 1 if exposure after infection (pre_time < 0).
# s = abs(pre_time)/24 in days, always non-negative.
mos_fit <- mos

rows <- mos_fit %>%
  distinct(exp_id, pre_time) %>%
  arrange(exp_id, pre_time) %>%
  mutate(row_idx = row_number(),
         s       = abs(pre_time) / 24,
         side    = as.integer(pre_time < 0))

cat("Row-level breakdown by side:\n")
print(rows %>% group_by(side) %>%
        summarise(n_rows = n(),
                  unique_s = paste(sort(unique(s)), collapse = ", "),
                  .groups = "drop"))

mos_fit <- mos_fit %>%
  left_join(rows %>% select(exp_id, pre_time, row_idx),
            by = c("exp_id", "pre_time")) %>%
  mutate(arm_int = ifelse(group == "control", 0L, 1L))

N <- nrow(rows)
M <- nrow(mos_fit)
K <- N

# Per-row prevalence summary (TBA fit).
prev_summary <- mos_fit %>%
  group_by(row_idx, arm_int) %>%
  summarise(n = n(), pos = sum(int > 0), .groups = "drop") %>%
  arrange(row_idx, arm_int)

n_ctl   <- as.integer(prev_summary$n[prev_summary$arm_int == 0])
pos_ctl <- as.integer(prev_summary$pos[prev_summary$arm_int == 0])
n_trt   <- as.integer(prev_summary$n[prev_summary$arm_int == 1])
pos_trt <- as.integer(prev_summary$pos[prev_summary$arm_int == 1])
stopifnot(length(n_ctl) == N, length(n_trt) == N)

cat(sprintf("Bidir split-nH fit on %d mosquitoes across %d rows (%d pre, %d post).\n",
            M, N, sum(rows$side == 0), sum(rows$side == 1)))

# Stan data lists (data + prior hyperparameters from priors_atn_main.R).
stan_data_tba <- c(list(
  N       = N,
  K       = K,
  exp_id  = rows$row_idx,
  s       = rows$s,
  side    = rows$side,
  n_ctl   = n_ctl,
  pos_ctl = pos_ctl,
  n_trt   = n_trt,
  pos_trt = pos_trt
), prior_data)

stan_data_tra <- c(list(
  N       = N,
  K       = K,
  M       = M,
  exp_id  = rows$row_idx,
  s       = rows$s,
  side    = rows$side,
  row_id  = mos_fit$row_idx,
  arm     = mos_fit$arm_int,
  y       = mos_fit$int
), prior_data)

# ---- Configs ----
configs <- list(
  tba = list(
    run        = RUN_TBA,
    stan_file  = "int_fit_linear_hill_bidir_splitnH.stan",
    label      = "TBA bidir split-nH (linear x Hill)",
    stan_data  = stan_data_tba,
    key_vars   = c("mu_p", "r_min", "b_max",
                   "s_half_pre", "log_s_half_pre",
                   "s_half_post", "log_s_half_post",
                   "nH_pre", "log_nH_pre",
                   "nH_post", "log_nH_post",
                   "sigma_exp"),
    pair_vars  = c("mu_p", "r_min",
                   "log_s_half_pre", "log_s_half_post",
                   "log_nH_pre", "log_nH_post", "sigma_exp"),
    likelihood = "linear"
  ),
  tra = list(
    run        = RUN_TRA,
    stan_file  = "int_fit_nb_rate_hill_bidir_splitnH.stan",
    label      = "TRA bidir split-nH (NB rate-thinning x Hill)",
    stan_data  = stan_data_tra,
    key_vars   = c("mu_L", "r_min", "b_max",
                   "s_half_pre", "log_s_half_pre",
                   "s_half_post", "log_s_half_post",
                   "nH_pre", "log_nH_pre",
                   "nH_post", "log_nH_post",
                   "phi", "reciprocal_phi", "sigma_exp"),
    pair_vars  = c("mu_L", "r_min",
                   "log_s_half_pre", "log_s_half_post",
                   "log_nH_pre", "log_nH_post",
                   "reciprocal_phi", "sigma_exp"),
    likelihood = "nb_rate"
  )
)

# ---- Fit one variant ----
fit_variant <- function(cfg, tag) {
  cat(sprintf("\n\n==== Fitting: %s ====\n", cfg$label))
  mod <- stan_model(file = cfg$stan_file)
  fit <- sampling(
    mod,
    data    = cfg$stan_data,
    seed    = 2025,
    chains  = 4,
    iter    = 3000,
    warmup  = 1000,
    control = list(adapt_delta = 0.999, max_treedepth = 14),
    refresh = 500
  )

  cat("\n-- HMC diagnostics --\n")
  check_hmc_diagnostics(fit)

  cat("\n-- Key parameter summary --\n")
  print(summary(fit, pars = cfg$key_vars)$summary)

  cat("\n-- OLRE offsets (beta_exp) --\n")
  print(summary(fit, pars = "beta_exp")$summary)

  post_arr <- as.array(fit, pars = cfg$pair_vars)
  p_pairs <- mcmc_pairs(post_arr,
                        off_diag_args = list(size = 0.4, alpha = 0.25))
  ggsave(sprintf("atn_main_bidir_splitnH_pairs_%s.png", tag), p_pairs,
         width = 10, height = 10, dpi = 200)

  p_trace <- mcmc_trace(post_arr)
  ggsave(sprintf("atn_main_bidir_splitnH_trace_%s.png", tag), p_trace,
         width = 11, height = 7, dpi = 200)

  if (cfg$likelihood == "linear") {
    pp_ctl <- rstan::extract(fit, pars = "pos_ctl_rep")$pos_ctl_rep
    pp_trt <- rstan::extract(fit, pars = "pos_trt_rep")$pos_trt_rep

    ppc_ctl <- ppc_intervals(y = pos_ctl, yrep = pp_ctl,
                             x = seq_len(N)) +
      labs(title    = sprintf("Control arm PPC -- %s", cfg$label),
           subtitle = "x = row index", x = "row", y = "oocyst-positive count")
    ggsave(sprintf("ppc_atn_main_bidir_splitnH_%s_ctl.png", tag), ppc_ctl,
           width = 9, height = 5, dpi = 200)

    ppc_trt <- ppc_intervals(y = pos_trt, yrep = pp_trt,
                             x = seq_len(N)) +
      labs(title    = sprintf("Treatment arm PPC -- %s", cfg$label),
           subtitle = "x = row index", x = "row", y = "oocyst-positive count")
    ggsave(sprintf("ppc_atn_main_bidir_splitnH_%s_trt.png", tag), ppc_trt,
           width = 9, height = 5, dpi = 200)
  } else {
    yrep_mat <- rstan::extract(fit, pars = "y_rep")$y_rep
    agg_idx  <- split(seq_len(M),
                      interaction(mos_fit$row_idx, mos_fit$arm_int,
                                  drop = TRUE))
    obs_tot <- do.call(rbind, lapply(agg_idx, function(ix) {
      data.frame(
        row_idx = mos_fit$row_idx[ix[1]],
        arm_int = mos_fit$arm_int[ix[1]],
        total   = sum(mos_fit$int[ix])
      )
    }))
    rep_tot <- sapply(agg_idx,
                      function(ix) rowSums(yrep_mat[, ix, drop = FALSE]))
    colnames(rep_tot) <- NULL

    idx_c <- which(obs_tot$arm_int == 0)
    idx_t <- which(obs_tot$arm_int == 1)

    ppc_ctl <- ppc_intervals(y = obs_tot$total[idx_c],
                             yrep = rep_tot[, idx_c, drop = FALSE],
                             x = obs_tot$row_idx[idx_c]) +
      labs(title    = sprintf("Control arm PPC (row totals) -- %s",
                              cfg$label),
           subtitle = "per-row sum of mosquito-level oocyst counts",
           x = "row", y = "row total (oocysts)")
    ggsave(sprintf("ppc_atn_main_bidir_splitnH_%s_ctl.png", tag), ppc_ctl,
           width = 9, height = 5, dpi = 200)

    ppc_trt <- ppc_intervals(y = obs_tot$total[idx_t],
                             yrep = rep_tot[, idx_t, drop = FALSE],
                             x = obs_tot$row_idx[idx_t]) +
      labs(title    = sprintf("Treatment arm PPC (row totals) -- %s",
                              cfg$label),
           subtitle = "per-row sum of mosquito-level oocyst counts",
           x = "row", y = "row total (oocysts)")
    ggsave(sprintf("ppc_atn_main_bidir_splitnH_%s_trt.png", tag), ppc_trt,
           width = 9, height = 5, dpi = 200)
  }

  # r(s) summary at key s values, on each side, using side-specific nH.
  cat("\n-- r(s) at key s values (smaller r = stronger blocking) --\n")
  post <- rstan::extract(fit, pars = c("r_min",
                                       "s_half_pre", "s_half_post",
                                       "nH_pre", "nH_post"))
  for (label in c("pre", "post")) {
    cat(sprintf("[side = %s]\n", label))
    s_half_d <- if (label == "pre") post$s_half_pre else post$s_half_post
    nH_d     <- if (label == "pre") post$nH_pre     else post$nH_post
    for (s_val in c(0, 0.5, 1, 2, 3)) {
      shn  <- s_half_d ^ nH_d
      kern <- shn / (shn + s_val ^ nH_d)
      r_s  <- 1 - (1 - post$r_min) * kern
      cat(sprintf("  s = %.2f d  2.5%% / 50%% / 97.5%%:  %.3f / %.3f / %.3f\n",
                  s_val,
                  quantile(r_s, 0.025, names = FALSE),
                  quantile(r_s, 0.5,   names = FALSE),
                  quantile(r_s, 0.975, names = FALSE)))
    }
  }

  saveRDS(fit, file = sprintf("atn_fit_%s_main_bidir_splitnH.rds", tag))
  write_priors_json(sprintf("atn_fit_%s_main_bidir_splitnH_priors.json", tag),
                    priors)
  cat(sprintf(paste0("\nSaved: atn_fit_%s_main_bidir_splitnH.rds, ",
                     "atn_fit_%s_main_bidir_splitnH_priors.json + ",
                     "diagnostic PNGs.\n"),
              tag, tag))

  invisible(fit)
}

# ---- Run enabled variants ----
for (nm in names(configs)) {
  cfg <- configs[[nm]]
  if (cfg$run) {
    fit_variant(cfg, nm)
  } else {
    cat(sprintf("\nSkipping: %s\n", cfg$label))
  }
}

cat("\n\nAll done.\n")
