# plot_atn_main_bidir_splitnH.R
# Inhibition-probability posterior plots for the BIDIRECTIONAL SPLIT-nH
# main-analysis ATN blocking-probability fits produced by
# fit_atn_main_bidir_splitnH.R.
#
# Identical in structure to plot_atn_main_bidir.R but reads (nH_pre,
# nH_post) from the fit and uses the side-specific value in the kernel.
#
# Each enabled variant generates 4 PNGs:
#   atn_main_bidir_splitnH_posterior_draws_<tag>.png
#   atn_main_bidir_splitnH_posterior_draws_with_data_<tag>.png
#   atn_main_bidir_splitnH_posterior_cri_<tag>.png
#   atn_main_bidir_splitnH_posterior_cri_with_data_<tag>.png
# where <tag> in {tba, tra}.
#
# The *_cri_with_data_<tag>.png plots additionally overlay the prior 95%
# band on b(s) as dashed lines (same colour as the posterior, alpha = 0.3).
# The split-nH prior gives an independent log_nH draw on each side, so the
# resulting prior band on the post side is wider than in the shared-nH
# variant -- a useful direct visual of how much extra "shape freedom" the
# split-nH model is offering.

suppressPackageStartupMessages({
  library(rstan)
  library(ggplot2)
  library(dplyr)
})

set.seed(2025)

# ---- Toggles ----
PLOT_TBA <- TRUE
PLOT_TRA <- TRUE

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

# ---- Priors (single source of truth: priors_atn_main.R) ----
source("priors_atn_main.R")

data_path <- "./exp_data/mos_level_int_summary_pre_post.csv"

# ---- Curves grid (single uninterrupted grid; continuous at 0) ----
x_full <- seq(-10, 10, by = 0.02)

# ---- Model curves (Hill kernel, side-dependent s_half AND nH) ----
kernel_hill <- function(s_half, nH, s) s_half^nH / (s_half^nH + s^nH)

compute_b_all <- function(fit) {
  post <- rstan::extract(fit, pars = c("r_min",
                                       "s_half_pre", "s_half_post",
                                       "nH_pre", "nH_post"))
  n_draws <- length(post$r_min)
  abs_x   <- abs(x_full)
  is_post <- x_full > 0   # x = 0 routed via pre branch (kernel = 1 either way)
  b_all   <- matrix(NA_real_, nrow = n_draws, ncol = length(x_full))
  for (k in seq_len(n_draws)) {
    r_min_k <- post$r_min[k]
    s_pre_k <- post$s_half_pre[k]
    s_pst_k <- post$s_half_post[k]
    nH_pre_k  <- post$nH_pre[k]
    nH_post_k <- post$nH_post[k]
    kern    <- ifelse(is_post,
                      kernel_hill(s_pst_k, nH_post_k, abs_x),
                      kernel_hill(s_pre_k, nH_pre_k,  abs_x))
    b_all[k, ] <- (1 - r_min_k) * kern
  }
  b_all
}

# ---- Prior 95% band on b(s) (bidir split-nH kernel) ------------------
# Draws (r_min, s_half_pre, s_half_post, nH_pre, nH_post) directly from
# the priors in priors_atn_main.R via sample_kernel_prior(variant =
# "bidir_split").  log_nH_pre and log_nH_post are independent draws here,
# so the resulting prior band on the post side is wider than in the
# shared-nH variant.
prior_b_band <- function(n_prior = 10000, seed = 2025) {
  pp      <- sample_kernel_prior(n_prior, variant = "bidir_split",
                                 seed = seed)
  abs_x   <- abs(x_full)
  is_post <- x_full > 0
  b_all   <- matrix(NA_real_, nrow = n_prior, ncol = length(x_full))
  for (k in seq_len(n_prior)) {
    kern <- ifelse(is_post,
                   kernel_hill(pp$s_half_post[k], pp$nH_post[k], abs_x),
                   kernel_hill(pp$s_half_pre[k],  pp$nH_pre[k],  abs_x))
    b_all[k, ] <- pp$b_max[k] * kern
  }
  data.frame(
    x     = x_full,
    lower = 100 * apply(b_all, 2, quantile, probs = 0.025),
    upper = 100 * apply(b_all, 2, quantile, probs = 0.975)
  )
}

df_prior <- prior_b_band()

# ---- Data overlays ----
mos <- read.csv(data_path)
stopifnot(all(mos$group %in% c("control", "200mg")),
          all(mos$int >= 0))
mos$x_days <- mos$post_time / 24

B   <- 4000
eps <- 1e-8

bootstrap_tra_row <- function(y_ctl, y_trt) {
  b_draws <- numeric(B)
  for (b in seq_len(B)) {
    yc <- sample(y_ctl, length(y_ctl), replace = TRUE)
    yt <- sample(y_trt, length(y_trt), replace = TRUE)
    mc <- mean(yc)
    mt <- mean(yt)
    b_draws[b] <- 1 - mt / max(mc, eps)
  }
  b_draws
}

jeffreys_tba_row <- function(y_ctl, y_trt) {
  pos_c <- sum(y_ctl > 0); n_c <- length(y_ctl)
  pos_t <- sum(y_trt > 0); n_t <- length(y_trt)
  p_c   <- rbeta(B, pos_c + 0.5, n_c - pos_c + 0.5)
  p_t   <- rbeta(B, pos_t + 0.5, n_t - pos_t + 0.5)
  1 - p_t / pmax(p_c, eps)
}

rows <- mos %>%
  distinct(exp_id, pre_time, post_time, x_days) %>%
  arrange(exp_id)

build_overlay <- function(draw_fn) {
  bind_rows(lapply(seq_len(nrow(rows)), function(i) {
    r_exp  <- rows$exp_id[i]
    r_pre  <- rows$pre_time[i]
    r_xd   <- rows$x_days[i]
    sub    <- mos[mos$exp_id == r_exp & mos$pre_time == r_pre, , drop = FALSE]
    y_ctl  <- sub$int[sub$group == "control"]
    y_trt  <- sub$int[sub$group == "200mg"]

    b_draws <- draw_fn(y_ctl, y_trt)
    data.frame(
      x     = r_xd,
      b_est = max(median(b_draws), 0),
      b_lo  = max(quantile(b_draws, 0.025, names = FALSE), 0),
      b_hi  = min(quantile(b_draws, 0.975, names = FALSE), 1)
    )
  }))
}

overlay_tba <- build_overlay(jeffreys_tba_row)
overlay_tra <- build_overlay(bootstrap_tra_row)

jitter_w <- 0.15
add_jitter <- function(df, seed = 2025) {
  set.seed(seed)
  df$x <- df$x + runif(nrow(df), -jitter_w / 2, jitter_w / 2)
  df
}
overlay_tba <- add_jitter(overlay_tba)
overlay_tra <- add_jitter(overlay_tra)

# ---- Shared plot components ----
y_floor_pct <- -12.0

make_ann_layers <- function() {
  list(
    annotate("segment", x = 0.2, xend = 10, y = -5.0, yend = -5.0,
             arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
             colour = "black"),
    annotate("text",    x = 4.6, y = -10.0,
             label = "Exposure after infection", vjust = 0),
    annotate("segment", x = -0.2, xend = -10, y = -5.0, yend = -5.0,
             arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
             colour = "black"),
    annotate("text",    x = -4.6, y = -10.0,
             label = "Exposure before infection", vjust = 0)
  )
}

point_colour <- "grey20"
prior_colour <- "dodgerblue"
prior_alpha  <- 0.3
prior_lw     <- 0.5

data_layers <- function(overlay_df) {
  list(
    geom_errorbar(
      data = overlay_df,
      aes(x = x, ymin = 100 * b_lo, ymax = 100 * b_hi),
      colour = point_colour, width = 0.1, alpha = 0.5, linewidth = 0.4,
      inherit.aes = FALSE
    ),
    geom_point(
      data = overlay_df,
      aes(x = x, y = 100 * b_est),
      colour = point_colour, size = 2.2, alpha = 0.85, stroke = 0.8,
      inherit.aes = FALSE
    )
  )
}

prior_layers <- function() {
  list(
    geom_line(
      data = df_prior, aes(x = x, y = lower),
      colour = prior_colour, linetype = "dashed",
      linewidth = prior_lw, alpha = prior_alpha,
      inherit.aes = FALSE
    ),
    geom_line(
      data = df_prior, aes(x = x, y = upper),
      colour = prior_colour, linetype = "dashed",
      linewidth = prior_lw, alpha = prior_alpha,
      inherit.aes = FALSE
    )
  )
}

common_x <- labs(x = "Time of exposure relative to infection (days)")
common_y <- coord_cartesian(ylim = c(y_floor_pct, 100))

theme_transparent <- theme(
  plot.background  = element_rect(fill = "transparent", colour = NA),
  panel.background = element_rect(fill = "transparent", colour = NA)
)

# ---- Per-variant plotting ----
plot_variant <- function(tag, y_lab, overlay_df) {
  rds_path <- sprintf("atn_fit_%s_main_bidir_splitnH.rds", tag)
  if (!file.exists(rds_path)) {
    warning(sprintf("Skipping '%s': %s not found.", tag, rds_path))
    return(invisible(NULL))
  }
  cat(sprintf("\n-- %s --\n", tag))
  fit       <- readRDS(rds_path)
  b_all     <- compute_b_all(fit)
  n_draws   <- nrow(b_all)
  b_all_pct <- 100 * b_all

  idx_sub <- sample(n_draws, min(100, n_draws))
  df_draws <- do.call(rbind, lapply(idx_sub, function(k) {
    data.frame(x = x_full, b = b_all_pct[k, ], draw = k)
  }))
  df_cri <- data.frame(
    x      = x_full,
    median = apply(b_all_pct, 2, median),
    lower  = apply(b_all_pct, 2, quantile, probs = 0.025),
    upper  = apply(b_all_pct, 2, quantile, probs = 0.975)
  )

  y_lab_layer <- labs(y = y_lab)

  p_draws_base <- ggplot(df_draws, aes(x = x, y = b, group = draw)) +
    geom_line(alpha = 0.18, colour = "dodgerblue") +
    make_ann_layers() +
    common_x + y_lab_layer + common_y +
    theme_minimal() + theme_transparent

  p_draws_data <- ggplot(df_draws, aes(x = x, y = b, group = draw)) +
    geom_line(alpha = 0.18, colour = "dodgerblue") +
    data_layers(overlay_df) +
    make_ann_layers() +
    common_x + y_lab_layer + common_y +
    theme_minimal() + theme_transparent

  p_cri_base <- ggplot(df_cri, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.3, fill = "dodgerblue") +
    geom_line(aes(y = median), colour = "dodgerblue", linewidth = 0.8) +
    make_ann_layers() +
    common_x + y_lab_layer + common_y +
    theme_minimal() + theme_transparent

  p_cri_data <- ggplot(df_cri, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.3, fill = "dodgerblue") +
    geom_line(aes(y = median), colour = "dodgerblue", linewidth = 0.8) +
    prior_layers() +
    data_layers(overlay_df) +
    make_ann_layers() +
    common_x + y_lab_layer + common_y +
    theme_minimal() + theme_transparent

  ggsave(sprintf("atn_main_bidir_splitnH_posterior_draws_%s.png", tag),
         p_draws_base, width = 8, height = 5, dpi = 450, bg = "transparent")
  ggsave(sprintf("atn_main_bidir_splitnH_posterior_draws_with_data_%s.png", tag),
         p_draws_data, width = 8, height = 5, dpi = 450, bg = "transparent")
  ggsave(sprintf("atn_main_bidir_splitnH_posterior_cri_%s.png", tag),
         p_cri_base,   width = 8, height = 5, dpi = 450, bg = "transparent")
  ggsave(sprintf("atn_main_bidir_splitnH_posterior_cri_with_data_%s.png", tag),
         p_cri_data,   width = 8, height = 5, dpi = 450, bg = "transparent")
  cat(sprintf("Saved 4 plots for '%s'.\n", tag))
}

# ---- Dispatch ----
variants <- list(
  list(run = PLOT_TBA, tag = "tba",
       y_lab = "TBA (%)", overlay = overlay_tba),
  list(run = PLOT_TRA, tag = "tra",
       y_lab = "TRA (%)", overlay = overlay_tra)
)

for (v in variants) {
  if (v$run) {
    plot_variant(v$tag, v$y_lab, v$overlay)
  } else {
    cat(sprintf("Skipping: %s\n", v$tag))
  }
}

cat("\nAll done.\n")
