# convert_tra_to_field_tba_bidir_splitnH.R
# BIDIRECTIONAL SPLIT-nH counterpart to convert_tra_to_field_tba_bidir.R:
# convert the bidir split-nH lab-fitted TRA posterior (from
# fit_atn_main_bidir_splitnH.R, nb_rate_hill bidir split-nH model) into
# field TBA via Bompard's closed-form relationship, Challenger et al. 2023
# JID 228:212-23 (SI p.4):
#
#   TBA(TRA; m, k) = [1 / (1 - (k/(k+m))^k)] *
#                    [ (k/(k+m*(1-TRA)))^k - (k/(k+m))^k ]
#
# Identical pipeline to the shared-nH bidir conversion script, but reads
# (nH_pre, nH_post) from the fit and uses the side-specific value in the
# kernel.
#
# Default Bompard parameters (Challenger 2023 SI, Burkina Faso wild
# mosquitoes, Bompard 2020):
#   m_central = 0.000157,  k = 0.00000495
#   sensitivity sweep on m: m/4, m, 2*m  (k held fixed).
#
# Inputs:
#   atn_fit_tra_main_bidir_splitnH.rds  (bidir split-nH TRA fit)
#   atn_fit_tba_main_bidir_splitnH.rds  (bidir split-nH lab TBA fit;
#                                        optional, three-curve plot only)
#
# Outputs (base, no overlay):
#   atn_main_bidir_splitnH_field_tba_cri.png
#   atn_main_bidir_splitnH_field_tba_sensitivity.png
#   atn_main_bidir_splitnH_three_curve_comparison.png
#   atn_main_bidir_splitnH_field_tba_summary.csv
#
# Outputs (with empirical overlays + prior 95% band):
#   atn_main_bidir_splitnH_field_tba_cri_with_data.png
#   atn_main_bidir_splitnH_field_tba_sensitivity_with_data.png   (overlay only;
#                                                                 no prior)
#   atn_main_bidir_splitnH_three_curve_comparison_with_data.png  (3-panel facet)
#
# Prior overlay convention: dashed lines at the 2.5%/97.5% prior quantiles
# of b(s), shown on *_cri_with_data and *_three_curve_comparison_with_data
# (omitted from the m-sensitivity plot to avoid clutter when comparing
# m-curves).  Same colour as the corresponding posterior, alpha = 0.3.

suppressPackageStartupMessages({
  library(rstan)
  library(ggplot2)
  library(dplyr)
  library(ggokabeito)
  library(viridisLite)
  library(ggnewscale)
})

okabe_ito <- c(
  black = "#000000",
  orange = "#E69F00",
  sky_blue = "#56B4E9",
  bluish_green = "#009E73",
  yellow = "#F0E442",
  blue = "#0072B2",
  vermillion = "#D55E00",
  purple = "#CC79A7"
)
#okabe_ito["sky_blue"]

okabe_ito_dark <- c(
  black = "#000000",
  orange = "#A66F00",
  sky_blue = "#2F7FAE",
  bluish_green = "#00664D",
  yellow = "#B3A300",
  blue = "#004C7A",
  vermillion = "#A04300",
  purple = "#8F4F73"
)

set.seed(2025)

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

# ---- Bompard parameters ----
m_central <- 0.000157
k_bompard <- 0.00000495
m_low     <- m_central / 4
m_high    <- m_central * 2

# ---- Grid (single, continuous through origin) ----
x_full <- seq(-10, 10, by = 0.02)

# ---- Hill kernel ----
kernel_hill <- function(s_half, nH, s) s_half^nH / (s_half^nH + s^nH)

# ---- Numerically stable Bompard TRA -> TBA (same as shared-nH version) ----
tba_from_tra <- function(tra, m, k) {
  tra_safe <- pmin(pmax(tra, 0), 1 - 1e-12)
  log_a    <- k * (log(k) - log(k + m))
  log_b    <- k * (log(k) - log(k + m * (1 - tra_safe)))
  num      <- exp(log_a) * expm1(log_b - log_a)
  denom    <- -expm1(log_a)
  out      <- num / denom
  out[tra_safe <= 0] <- 0
  pmin(pmax(out, 0), 1)
}

# ---- Helper: load split-nH bidir Hill fit and compute b(s) draws ----
# Reads s_half_pre, s_half_post, nH_pre, nH_post and selects kernel by
# sign of x.  Continuous at x = 0 since both kernels equal 1 there.
compute_b_all <- function(fit) {
  post    <- rstan::extract(fit, pars = c("r_min",
                                          "s_half_pre", "s_half_post",
                                          "nH_pre", "nH_post"))
  n_draws <- length(post$r_min)
  abs_x   <- abs(x_full)
  is_post <- x_full > 0
  b_all   <- matrix(NA_real_, nrow = n_draws, ncol = length(x_full))
  for (k_d in seq_len(n_draws)) {
    r_min_k   <- post$r_min[k_d]
    s_pre_k   <- post$s_half_pre[k_d]
    s_pst_k   <- post$s_half_post[k_d]
    nH_pre_k  <- post$nH_pre[k_d]
    nH_post_k <- post$nH_post[k_d]
    kern    <- ifelse(is_post,
                      kernel_hill(s_pst_k, nH_post_k, abs_x),
                      kernel_hill(s_pre_k, nH_pre_k,  abs_x))
    b_all[k_d, ] <- (1 - r_min_k) * kern
  }
  b_all
}

# ---- Prior 95% band on b(s) (bidir split-nH kernel) ----------------------
# Sample directly from the Stan-file priors (split-nH variant):
#   r_min            ~ Beta(1, 1)
#   log_s_half_pre   ~ N(log(1), 0.75)
#   log_s_half_post  ~ N(log(1), 0.75)
#   log_nH_pre       ~ N(log(2), 0.5)
#   log_nH_post      ~ N(log(2), 0.5)
# OLRE / mu priors do not affect b(s) and are excluded.  The same lab-scale
# prior on b(s) applies to TRA and TBA fits.

prior_b_matrix <- function(n_prior = 10000, seed = 2025) {
  set.seed(seed)
  r_min       <- rbeta(n_prior, 0.5, 0.5)
  b_max       <- 1 - r_min
  # Tighter s_half (centered at 1 day, ~95% in [0.46, 2.17])
  s_half_pre  <- exp(rnorm(n_prior, log(2), 1))
  s_half_post <- exp(rnorm(n_prior, log(2), 1))
  # Higher, modestly wider nH (median 5, ~95% in [2.3, 10.9])
  nH_pre      <- exp(rnorm(n_prior, log(5.0), 1))
  nH_post     <- exp(rnorm(n_prior, log(5.0), 1))
  abs_x       <- abs(x_full)
  is_post     <- x_full > 0
  b_all       <- matrix(NA_real_, nrow = n_prior, ncol = length(x_full))
  for (k in seq_len(n_prior)) {
    kern <- ifelse(is_post,
                   kernel_hill(s_half_post[k], nH_post[k], abs_x),
                   kernel_hill(s_half_pre[k],  nH_pre[k],  abs_x))
    b_all[k, ] <- b_max[k] * kern
  }
  b_all
}
# prior_b_matrix <- function(n_prior = 10000, seed = 2025) {
#   set.seed(seed)
#   r_min       <- rbeta(n_prior, 1, 1)
#   b_max       <- 1 - r_min
#   s_half_pre  <- exp(rnorm(n_prior, log(1.0), 0.75))
#   s_half_post <- exp(rnorm(n_prior, log(1.0), 0.75))
#   nH_pre      <- exp(rnorm(n_prior, log(2.0), 0.5))
#   nH_post     <- exp(rnorm(n_prior, log(2.0), 0.5))
#   abs_x       <- abs(x_full)
#   is_post     <- x_full > 0
#   b_all       <- matrix(NA_real_, nrow = n_prior, ncol = length(x_full))
#   for (k in seq_len(n_prior)) {
#     kern <- ifelse(is_post,
#                    kernel_hill(s_half_post[k], nH_post[k], abs_x),
#                    kernel_hill(s_half_pre[k],  nH_pre[k],  abs_x))
#     b_all[k, ] <- b_max[k] * kern
#   }
#   b_all
# }

quant_band <- function(b_mat) {
  data.frame(
    x     = x_full,
    lower = 100 * apply(b_mat, 2, quantile, probs = 0.025),
    upper = 100 * apply(b_mat, 2, quantile, probs = 0.975)
  )
}

prior_b_lab    <- prior_b_matrix()
prior_b_field  <- matrix(tba_from_tra(as.vector(prior_b_lab),
                                      m_central, k_bompard),
                         nrow = nrow(prior_b_lab), ncol = ncol(prior_b_lab))
df_prior_lab   <- quant_band(prior_b_lab)
df_prior_field <- quant_band(prior_b_field)

# ---- Load TRA fit ----
tra_path <- "atn_fit_tra_main_bidir_splitnH.rds"
stopifnot(file.exists(tra_path))
fit_tra   <- readRDS(tra_path)
b_tra_all <- compute_b_all(fit_tra)
n_draws   <- nrow(b_tra_all)
cat(sprintf("Bidir split-nH TRA fit loaded: %d posterior draws.\n", n_draws))

# ---- Apply Bompard conversion at three m values ----
convert_to_field_tba <- function(tra_mat, m, k) {
  matrix(tba_from_tra(as.vector(tra_mat), m, k),
         nrow = nrow(tra_mat), ncol = ncol(tra_mat))
}

b_tba_field      <- convert_to_field_tba(b_tra_all, m_central, k_bompard)
b_tba_field_low  <- convert_to_field_tba(b_tra_all, m_low,     k_bompard)
b_tba_field_high <- convert_to_field_tba(b_tra_all, m_high,    k_bompard)

# ---- Posterior summary helper (output in percent, 0-100) ----
summarise_b <- function(b_mat) {
  data.frame(
    x      = x_full,
    median = 100 * apply(b_mat, 2, median),
    lower  = 100 * apply(b_mat, 2, quantile, probs = 0.025),
    upper  = 100 * apply(b_mat, 2, quantile, probs = 0.975)
  )
}

df_tra       <- summarise_b(b_tra_all)
df_tba_field <- summarise_b(b_tba_field)
df_tba_low   <- summarise_b(b_tba_field_low)
df_tba_high  <- summarise_b(b_tba_field_high)

# ---- Lab TBA fit (optional) ----
tba_path   <- "atn_fit_tba_main_bidir_splitnH.rds"
df_tba_lab <- NULL
if (file.exists(tba_path)) {
  fit_tba_lab <- readRDS(tba_path)
  b_tba_lab   <- compute_b_all(fit_tba_lab)
  df_tba_lab  <- summarise_b(b_tba_lab)
  cat(sprintf("Bidir split-nH lab TBA fit loaded: %d posterior draws.\n",
              nrow(b_tba_lab)))
} else {
  warning(sprintf("%s not found; three-curve plot will omit lab TBA.",
                  tba_path))
}

# ---- Empirical overlays ---------------------------------------------------
data_path <- "./exp_data/mos_level_int_summary_pre_post_full.csv"
mos       <- read.csv(data_path)
stopifnot(all(mos$group %in% c("control", "200mg")), all(mos$int >= 0))
mos$x_days <- mos$post_time / 24

B_boot <- 4000
eps    <- 1e-8

bootstrap_tra_row <- function(y_ctl, y_trt) {
  b_draws <- numeric(B_boot)
  for (b in seq_len(B_boot)) {
    yc <- sample(y_ctl, length(y_ctl), replace = TRUE)
    yt <- sample(y_trt, length(y_trt), replace = TRUE)
    b_draws[b] <- 1 - mean(yt) / max(mean(yc), eps)
  }
  b_draws
}

jeffreys_tba_row <- function(y_ctl, y_trt) {
  pos_c <- sum(y_ctl > 0); n_c <- length(y_ctl)
  pos_t <- sum(y_trt > 0); n_t <- length(y_trt)
  p_c   <- rbeta(B_boot, pos_c + 0.5, n_c - pos_c + 0.5)
  p_t   <- rbeta(B_boot, pos_t + 0.5, n_t - pos_t + 0.5)
  1 - p_t / pmax(p_c, eps)
}

field_tba_draws_row <- function(y_ctl, y_trt) {
  tba_from_tra(bootstrap_tra_row(y_ctl, y_trt), m_central, k_bompard)
}

rows <- mos %>%
  distinct(exp_id, pre_time, post_time, x_days) %>%
  arrange(exp_id)

build_overlay <- function(draw_fn) {
  bind_rows(lapply(seq_len(nrow(rows)), function(i) {
    r_exp <- rows$exp_id[i]
    r_pre <- rows$pre_time[i]
    r_xd  <- rows$x_days[i]
    sub   <- mos[mos$exp_id == r_exp & mos$pre_time == r_pre, , drop = FALSE]
    yc    <- sub$int[sub$group == "control"]
    yt    <- sub$int[sub$group == "200mg"]
    d     <- draw_fn(yc, yt)
    data.frame(
      x     = r_xd,
      b_est = max(median(d), 0),
      b_lo  = max(quantile(d, 0.025, names = FALSE), 0),
      b_hi  = min(quantile(d, 0.975, names = FALSE), 1)
    )
  }))
}

overlay_tra       <- build_overlay(bootstrap_tra_row)
overlay_field_tba <- build_overlay(field_tba_draws_row)
overlay_tba_lab   <- build_overlay(jeffreys_tba_row)

# jitter_w <- 0.15
jitter_w <- 0.3
#add_jitter <- function(df, seed = 2025) {
add_jitter <- function(df, seed = 123) {
  set.seed(seed)
  df$x <- df$x + runif(nrow(df), -jitter_w / 2, jitter_w / 2)
  df
}
overlay_tra       <- add_jitter(overlay_tra)
overlay_field_tba <- add_jitter(overlay_field_tba)
overlay_tba_lab   <- add_jitter(overlay_tba_lab)

# ---- Shared plot components ----
y_floor_pct <- -12.0

make_ann_layers <- function() {
  list(
    annotate("segment", x = 0.2, xend = 4, y = -5.0, yend = -5.0,
             arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
             colour = "black"),
    annotate("text",    x = 2.1, y = -10.0,
             label = "Exposure after infection", vjust = 0),
    annotate("segment", x = -0.2, xend = -4, y = -5.0, yend = -5.0,
             arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
             colour = "black"),
    annotate("text",    x = -2.1, y = -10.0,
             label = "Exposure before infection", vjust = 0)
  )
}

common_x <- labs(x = "Time of exposure relative to infection (days)")
common_y <- coord_cartesian(ylim = c(y_floor_pct, 100))

theme_transparent <- theme(
  plot.background  = element_rect(fill = "transparent", colour = NA),
  panel.background = element_rect(fill = "transparent", colour = NA)
)

point_colour <- "grey20"
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

# Single-panel field-TBA prior band (darkorange, matching p_field).
prior_layers_field <- function() {
  list(
    geom_line(
      data = df_prior_field, aes(x = x, y = lower),
      colour = "darkorange", linetype = "dashed",
      linewidth = prior_lw, alpha = prior_alpha,
      inherit.aes = FALSE
    ),
    geom_line(
      data = df_prior_field, aes(x = x, y = upper),
      colour = "darkorange", linetype = "dashed",
      linewidth = prior_lw, alpha = prior_alpha,
      inherit.aes = FALSE
    )
  )
}

# ---- Plot 1: field TBA posterior median + 95% CrI ----
field_TBA_col <- ggokabeito::palette_okabe_ito()[2]
p_field <- ggplot(df_tba_field, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.3, fill = field_TBA_col) +
  geom_line(aes(y = median), colour = field_TBA_col, linewidth = 0.8) +
  make_ann_layers() +
  common_x + labs(y = "Field TBA (%)") + common_y +
  theme_minimal() + theme_transparent

ggsave("atn_main_bidir_splitnH_field_tba_cri.png", p_field,
       width = 8, height = 5, dpi = 450, bg = "transparent")

# ---- Plot 2: sensitivity to m (m/4, m, 2m) ----
mako_disc <- mako(6)  # pick number of categories
mako_disc <- mako_disc[3:5]

sens_disc <- c("#395D9CFF", "#56B4E9", "#60CEACFF")
sens_disc <- c("slateblue2", "#56B4E9", "seagreen3")

df_sens <- bind_rows(
  df_tba_low   %>% mutate(label = sprintf("m = %.4g", m_low)),
  df_tba_field %>% mutate(label = sprintf("m = %.4g",       m_central)),
  df_tba_high  %>% mutate(label = sprintf("m = %.4g",   m_high))
)
df_sens$label <- factor(df_sens$label, levels = unique(df_sens$label))

p_sens <- ggplot(df_sens, aes(x = x, colour = label, fill = label)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  geom_line(aes(y = median), linewidth = 0.8) +
  scale_colour_manual(values = sens_disc) +
  scale_fill_manual(values = sens_disc) +
  make_ann_layers() +
  common_x + labs(y = "Field TBA (%)", colour = NULL, fill = NULL) +
  common_y +
  theme_minimal() + theme_transparent +
  theme(legend.position = "bottom")

ggsave("atn_main_bidir_splitnH_field_tba_sensitivity.png", p_sens,
       width = 8, height = 5.5, dpi = 450, bg = "transparent")

# ---- Plot 3: three-curve comparison (legend mode) ----
df_compare <- bind_rows(
  df_tra       %>% mutate(curve = "Lab TRA"),
  df_tba_field %>% mutate(curve = "Field TBA"),
  if (!is.null(df_tba_lab)) df_tba_lab %>% mutate(curve = "Lab TBA")
)
df_compare$curve <- factor(df_compare$curve, levels = unique(df_compare$curve))

p_compare <- ggplot(df_compare, aes(x = x, colour = curve, fill = curve)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  geom_line(aes(y = median), linewidth = 0.8) +
  make_ann_layers() +
  common_x + labs(y = "Blocking probability (%)",
                  colour = NULL, fill = NULL) +
  common_y +
  theme_minimal() + theme_transparent +
  # scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  # scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_colour_manual(values = c("#C89B2B", "#EE7733", "#9A4D5E")) +
  scale_fill_manual(values = c("#C89B2B", "#EE7733", "#9A4D5E")) +
  theme(legend.position = "bottom")

ggsave("atn_main_bidir_splitnH_three_curve_comparison.png", p_compare,
       width = 8, height = 5.5, dpi = 450, bg = "white")

# ---- Plot 1b: field TBA CrI with overlay + prior band ----
p_field_data <- p_field + prior_layers_field() + data_layers(overlay_field_tba)
ggsave("atn_main_bidir_splitnH_field_tba_cri_with_data.png", p_field_data,
       width = 8, height = 5, dpi = 450, bg = "transparent")

# ---- Plot 2b: sensitivity with central-m overlay (no prior band) ----
p_sens_data <- p_sens + data_layers(overlay_field_tba)
ggsave("atn_main_bidir_splitnH_field_tba_sensitivity_with_data.png", p_sens_data,
       width = 8, height = 5.5, dpi = 450, bg = "transparent")

# ---- Plot 3b: three-curve comparison faceted with overlays + prior bands ----
# Each panel shows one model curve, its model-matched empirical overlay,
# and the panel's prior 95% band as dashed lines (lab-scale prior for
# Lab TRA / Lab TBA, field-scale prior at central m for Field TBA).
# Panel y-axis labels (e.g. "Lab TRA (%)") are rendered via strip.position
# = "left" + a labeller; no shared y-axis title.
panel_levels <- c("Lab TRA", "Field TBA", "Lab TBA")
panel_labels <- c("Lab TRA"   = "Lab Transmission Reduction Activity (%)",
                  "Field TBA" = "Field Transmission Blocking Activity (%)",
                  "Lab TBA"   = "Lab Transmission Blocking Activity (%)")

curve_blocks <- list(
  df_tra       %>% mutate(panel = panel_levels[1]),
  df_tba_field %>% mutate(panel = panel_levels[2])
)
overlay_blocks <- list(
  overlay_tra       %>% mutate(panel = panel_levels[1]),
  overlay_field_tba %>% mutate(panel = panel_levels[2])
)
prior_blocks <- list(
  df_prior_lab   %>% mutate(panel = panel_levels[1]),
  df_prior_field %>% mutate(panel = panel_levels[2])
)
if (!is.null(df_tba_lab)) {
  curve_blocks[[3]]   <- df_tba_lab     %>% mutate(panel = panel_levels[3])
  overlay_blocks[[3]] <- overlay_tba_lab %>% mutate(panel = panel_levels[3])
  prior_blocks[[3]]   <- df_prior_lab   %>% mutate(panel = panel_levels[3])
}
df_facet_curve   <- bind_rows(curve_blocks)
df_facet_overlay <- bind_rows(overlay_blocks)
df_facet_prior   <- bind_rows(prior_blocks)
present_levels   <- intersect(panel_levels, unique(df_facet_curve$panel))
df_facet_curve$panel   <- factor(df_facet_curve$panel,   levels = present_levels)
df_facet_overlay$panel <- factor(df_facet_overlay$panel, levels = present_levels)
df_facet_prior$panel   <- factor(df_facet_prior$panel,   levels = present_levels)

p_facet <- ggplot(df_facet_curve, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = panel),
              alpha = 0.25, colour = NA) +
  geom_line(aes(y = median, colour = panel), linewidth = 0.8) +
  geom_line(
    data = df_facet_prior, aes(x = x, y = lower, colour = panel),
    linetype = "dashed", linewidth = prior_lw, alpha = prior_alpha,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_line(
    data = df_facet_prior, aes(x = x, y = upper, colour = panel),
    linetype = "dashed", linewidth = prior_lw, alpha = prior_alpha,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_errorbar(
    data = df_facet_overlay,
    aes(x = x, ymin = 100 * b_lo, ymax = 100 * b_hi),
    colour = point_colour, width = 0.1, alpha = 0.5, linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = df_facet_overlay,
    # aes(x = x, y = 100 * b_est, colour = panel),
    aes(x = x, y = 100 * b_est),
    colour = point_colour, size = 2.0, alpha = 0.85, stroke = NA,
    # colour = point_colour, size = 2.0, alpha = 0.85, stroke = 0.8,
    inherit.aes = FALSE
  ) +
  scale_color_okabe_ito() +
  scale_fill_okabe_ito() +
  facet_wrap(~ panel, ncol = 1, strip.position = "left",
             labeller = labeller(panel = panel_labels)) +
  make_ann_layers() +
  common_x + labs(y = NULL) + common_y +
  theme_minimal() + theme_transparent +
  theme(legend.position   = "none",
        strip.placement   = "outside",
        strip.background  = element_blank(),
        strip.text.y.left = element_text(angle = 90, face = "bold"))

ggsave("atn_main_bidir_splitnH_three_curve_comparison_with_data.png", p_facet,
       width = 8, height = 11, dpi = 450, bg = "white")




p_facet <- ggplot(df_facet_curve, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = panel),
              alpha = 0.25, colour = NA) +
  geom_line(aes(y = median, colour = panel), linewidth = 0.8) +
  #geom_line(
  #  data = df_facet_prior, aes(x = x, y = lower, colour = panel),
  #  linetype = "dashed", linewidth = prior_lw, alpha = prior_alpha,
  #  inherit.aes = FALSE, show.legend = FALSE
  #) +
  #geom_line(
  #  data = df_facet_prior, aes(x = x, y = upper, colour = panel),
  #  linetype = "dashed", linewidth = prior_lw, alpha = prior_alpha,
  #  inherit.aes = FALSE, show.legend = FALSE
  #) +
  # scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  # scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_colour_manual(values = c("#C89B2B", "#EE7733", "#9A4D5E")) +
  scale_fill_manual(values = c("#C89B2B", "#EE7733", "#9A4D5E")) +

  ggnewscale::new_scale_colour() +

  geom_errorbar(
    data = df_facet_overlay,
    aes(x = x, ymin = 100 * b_lo, ymax = 100 * b_hi, colour = panel),
    width = 0.1, alpha = 0.5, linewidth = 0.4,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_point(
    data = df_facet_overlay,
    aes(x = x, y = 100 * b_est, colour = panel),
    size = 2.0, alpha = 0.85, stroke = NA,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  #scale_colour_manual(values = c("#A66F00", "#2F7FAE", "#00664D")) +
  scale_colour_manual(values = c("#A77816", "#B85A26", "#7A2F3E")) +

  facet_wrap(~ panel, ncol = 1, strip.position = "left",
             labeller = labeller(panel = panel_labels)) +
  make_ann_layers() +
  common_x + labs(y = NULL) + common_y +
  theme_minimal() + theme_transparent +
  scale_x_continuous(
    breaks = seq(-4, 4, 1),
    limits = c(-4, 4)
  ) +
  theme(
    legend.position   = "none",
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.y.left = element_text(angle = 90, face = "bold")
  )

ggsave("test_atn_main_bidir_splitnH_three_curve_comparison_with_data.png", p_facet,
       width = 8, height = 11, dpi = 450, bg = "white")



# Comparison test

p_compare_data <- p_compare +
  geom_errorbar(
    data = df_facet_overlay,
    aes(x = x, ymin = 100 * b_lo, ymax = 100 * b_hi, colour = panel),
    width = 0.1, alpha = 0.5, linewidth = 0.4,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_point(
    data = df_facet_overlay,
    aes(x = x, y = 100 * b_est, colour = panel),
    size = 2.0, alpha = 0.85, stroke = NA,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(-4, 4, 1),
    limits = c(-4, 4)
  ) +
  #scale_colour_manual(values = c("#A66F00", "#2F7FAE", "#00664D")) +
  scale_colour_manual(values = c("#A77816", "#B85A26", "#7A2F3E"))

ggsave("test_three_curve_comparison_with_data.png", p_compare_data,
       width = 8, height = 5, dpi = 450, bg = "white")

# ---- Summary table at key |s| values, both sides ----
key_s <- c(0, 0.5, 1, 2, 3, 5)

get_at_x <- function(b_mat, x_val) {
  idx  <- which.min(abs(x_full - x_val))
  vals <- b_mat[, idx]
  c(median = median(vals),
    lower  = quantile(vals, 0.025, names = FALSE),
    upper  = quantile(vals, 0.975, names = FALSE))
}

build_summary_side <- function(side_label, sign) {
  do.call(rbind, lapply(key_s, function(s_val) {
    x_val  <- sign * s_val
    q_tra  <- get_at_x(b_tra_all,        x_val)
    q_fld  <- get_at_x(b_tba_field,      x_val)
    q_lo   <- get_at_x(b_tba_field_low,  x_val)
    q_hi   <- get_at_x(b_tba_field_high, x_val)
    data.frame(
      side                 = side_label,
      s_days               = s_val,
      lab_TRA_med_pct      = round(100 * q_tra["median"], 2),
      lab_TRA_lo_pct       = round(100 * q_tra["lower"],  2),
      lab_TRA_hi_pct       = round(100 * q_tra["upper"],  2),
      field_TBA_med_pct    = round(100 * q_fld["median"], 2),
      field_TBA_lo_pct     = round(100 * q_fld["lower"],  2),
      field_TBA_hi_pct     = round(100 * q_fld["upper"],  2),
      field_TBA_low_m_med  = round(100 * q_lo["median"],  2),
      field_TBA_high_m_med = round(100 * q_hi["median"],  2)
    )
  }))
}

summary_tbl <- rbind(
  build_summary_side("pre",  -1),
  build_summary_side("post", +1)
)
rownames(summary_tbl) <- NULL

write.csv(summary_tbl, "atn_main_bidir_splitnH_field_tba_summary.csv",
          row.names = FALSE)
cat("\nSummary table:\n")
print(summary_tbl, row.names = FALSE)
cat("\nSaved 6 plots (3 base, 3 with empirical overlays) and 1 CSV.\n")
cat("All done.\n")
