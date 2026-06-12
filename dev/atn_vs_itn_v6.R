library(dplyr)
library(tidyr)
library(ggplot2)
library(viridisLite)
library(viridis)
library(RColorBrewer)

# User inputs

read_itn_pars <- TRUE

n_steps <- 13 * 365
old_distribution_times <- c(1, 4) * 365
new_distribution_times <- c(7, 10) * 365

retention_time <- 588

draws <- 50

seed <- 123

eir_levels <- c(0.7, 4.2, 26)
res_levels <- c(0.5, 0.9, 0.99)

old_coverage <- 0.9
new_coverage <- 0.9

use_sampled_atn_halflife <- TRUE

eir_labels <- c(
  "0.7" = "Low transmission",
  "4.2" = "Moderate transmission",
  "26"  = "High transmission"
)

# -----------------------------------------------

distribution_times <- c(old_distribution_times, new_distribution_times)
n_old_dist <- length(old_distribution_times)
n_new_dist <- length(new_distribution_times)
n_dist <- length(distribution_times)
old_covs <- rep(old_coverage, n_old_dist)
new_covs <- rep(new_coverage, n_new_dist)
all_covs <- c(old_covs, new_covs)
eir_labels <- eir_labels[as.character(eir_levels)]

# ------------------------------------------------

if (read_itn_pars) {
  only_pars <- readRDS("./dev/net_params/dat_res_only.rds") |>
    mutate(gamman = gamman * 365)
  pbo_pars  <- readRDS("./dev/net_params/dat_res_pbo.rds") |>
    mutate(gamman = gamman * 365)
  cfp_pars  <- readRDS("./dev/net_params/dat_res_cfp.rds") |>
    mutate(gamman = gamman * 365)
}

n_itn_pars <- max(only_pars$draw)

set.seed(seed); itn_draws <- floor(runif(draws, 1, n_itn_pars + 1))

# ------------------------------------------------

eip_rds = "./dev/eip_fit_hill_result.rds"
tra_rds = "./dev/atn_fit_tra_main_bidir_splitnH.rds"
tba_rds = "./dev/atn_fit_tba_main_bidir_splitnH.rds"

eip_samples <- rstan::extract(readRDS(eip_rds))
tra_samples <- rstan::extract(readRDS(tra_rds))
tba_samples <- rstan::extract(readRDS(tba_rds))

n_eip_samples <- length(eip_samples$s_half)
n_tra_samples <- length(tra_samples$s_half_pre)
n_tba_samples <- length(tba_samples$s_half_pre)

set.seed(seed); eip_draws <- floor(runif(draws, 1, n_eip_samples + 1))
set.seed(seed); tra_draws <- floor(runif(draws, 1, n_tra_samples + 1))
set.seed(seed); tba_draws <- floor(runif(draws, 1, n_tba_samples + 1))

gamman_samples <- c(only_pars$gamman, pbo_pars$gamman, cfp_pars$gamman)
n_gamman_samples <- length(gamman_samples)
set.seed(seed); gamma_atn_draws <- floor(runif(draws, 1, n_gamman_samples + 1))
atn_halflife <- 2.64 * 365

# ------------------------------------------------
# Storage
# ------------------------------------------------
all_runs <- list()

# ------------------------------------------------
# Loop over EIR and resistance levels
# ------------------------------------------------

ptm <- proc.time()

for (eir in eir_levels) {

  cat("\n==============================\n")
  cat("Starting EIR =", eir, "\n")
  cat("==============================\n")

  for (draw in 1:draws) {

    cat("\n==============================\n")
    cat("Starting draw =", draw, "of ", draws, "\n")
    cat("==============================\n")

    itn_draw <- itn_draws[draw]
    eip_draw <- eip_draws[draw]
    tra_draw <- tra_draws[draw]
    tba_draw <- tba_draws[draw]

    only_sample <- only_pars |> filter(draw == itn_draw)
    pbo_sample <- pbo_pars |> filter(draw == itn_draw)
    cfp_sample <- cfp_pars |> filter(draw == itn_draw)

    # EIP lengthening draws
    s_half_eip <- eip_samples$s_half[eip_draw]
    nH_eip <- eip_samples$nH[eip_draw]
    rho_frac <- eip_samples$r0_frac[eip_draw]

    # ATN pre-infection blocking (lab TRA -> field TBA assumed)
    s_half_pre <- tra_samples$s_half_pre[tra_draw]
    nH_pre <- tra_samples$nH_pre[tra_draw]
    B_max_post <- tra_samples$b_max[tra_draw]
    s_half_post <- tra_samples$s_half_post[tra_draw]
    nH_post <- tra_samples$nH_post[tra_draw]

    # ATN efficacy decay
    if (use_sampled_atn_halflife) {
      gamma_atn_draw <- gamma_atn_draws[draw]
      atn_halflife <- gamman_samples[gamma_atn_draw]
    }

    # ----- ITN loop -----
    for (res in res_levels) {

      cat("  -> Resistance =", res, "\n")

      only_res_sample <- filter(only_sample, resistance == res)
      pbo_res_sample <- filter(pbo_sample, resistance == res)
      cfp_res_sample <- filter(cfp_sample, resistance == res)

      # Baseline
      df_base <- get_parameters(
        n_days = n_steps
      ) |>
        set_bednets(
          days      = distribution_times,
          coverages = c(old_covs, rep(0, n_new_dist)),
          retention = retention_time,
          dn0 = matrix(rep(cfp_res_sample$dn0, n_dist), nrow = n_dist, ncol = 1),
          rn = matrix(rep(cfp_res_sample$rn0, n_dist), nrow = n_dist, ncol = 1),
          rnm = matrix(rep(.24, n_dist), nrow = n_dist, ncol = 1),
          gamman = rep(cfp_res_sample$gamman, n_dist)
        ) |>
        set_equilibrium(init_EIR = eir) |>
        run_aitn_simulation_v3() |>
        as.data.frame() |>
        mutate(
          itn_type            = "none",
          resistance          = res,
          EIR                 = eir,
          rep                 = draw,
          pfpr2to10           = n_detect_730_3650 / n_730_3650,
          base_clin_inc_0_Inf = n_clin_inc_0_Inf,
          cfp_clin_inc_0_Inf  = NA
        )

      # CFP
      df_cfp <- get_parameters(n_days = n_steps) |>
        set_bednets(
          days      = distribution_times,
          coverages = all_covs,
          retention = retention_time,
          dn0 = matrix(
            c(
              rep(cfp_res_sample$dn0, n_old_dist),
              rep(cfp_res_sample$dn0, n_new_dist)
            ),
            nrow = n_dist,
            ncol = 1
          ),
          rn = matrix(
            c(
              rep(cfp_res_sample$rn0, n_old_dist),
              rep(cfp_res_sample$rn0, n_new_dist)
            ),
            nrow = n_dist,
            ncol = 1
          ),
          rnm = matrix(rep(.24, n_dist), nrow = n_dist, ncol = 1),
          gamman = c(
            rep(cfp_res_sample$gamman, n_old_dist),
            rep(cfp_res_sample$gamman, n_new_dist)
          )
        ) |>
        set_equilibrium(init_EIR = eir)  |>
        run_aitn_simulation_v3() |>
        as.data.frame() |>
        mutate(
          itn_type            = "Pyr-CFP",
          resistance          = res,
          EIR                 = eir,
          rep                 = draw,
          pfpr2to10           = n_detect_730_3650 / n_730_3650,
          base_clin_inc_0_Inf = df_base |> pull(n_clin_inc_0_Inf),
          cfp_clin_inc_0_Inf  = n_clin_inc_0_Inf
        )

      # ONLY
      df_only <- get_parameters(n_days = n_steps) |>
        set_bednets(
          days      = distribution_times,
          coverages = all_covs,
          retention = retention_time,
          dn0 = matrix(
            c(
              rep(cfp_res_sample$dn0, n_old_dist),
              rep(only_res_sample$dn0, n_new_dist)
            ),
            nrow = n_dist,
            ncol = 1
          ),
          rn = matrix(
            c(
              rep(cfp_res_sample$rn0, n_old_dist),
              rep(only_res_sample$rn0, n_new_dist)
            ),
            nrow = n_dist,
            ncol = 1
          ),
          rnm = matrix(rep(.24, n_dist), nrow = n_dist, ncol = 1),
          gamman = c(
            rep(cfp_res_sample$gamman, n_old_dist),
            rep(only_res_sample$gamman, n_new_dist)
          )
        ) |>
        set_equilibrium(init_EIR = eir) |>
        run_aitn_simulation_v3() |>
        as.data.frame() |>
        mutate(
          itn_type            = "Pyr-only",
          resistance          = res,
          EIR                 = eir,
          rep                 = draw,
          pfpr2to10           = n_detect_730_3650 / n_730_3650,
          base_clin_inc_0_Inf = df_base |> pull(n_clin_inc_0_Inf),
          cfp_clin_inc_0_Inf  = df_cfp |> pull(n_clin_inc_0_Inf)
        )

      # PBO
      df_pbo <- get_parameters(n_days = n_steps) |>
        set_bednets(
          days      = distribution_times,
          coverages = all_covs,
          retention = retention_time,
          dn0 = matrix(
            c(
              rep(cfp_res_sample$dn0, n_old_dist),
              rep(pbo_res_sample$dn0, n_new_dist)
            ),
            nrow = n_dist,
            ncol = 1
          ),
          rn = matrix(
            c(
              rep(cfp_res_sample$rn0, n_old_dist),
              rep(pbo_res_sample$rn0, n_new_dist)
            ),
            nrow = n_dist,
            ncol = 1
          ),
          rnm = matrix(rep(.24, n_dist), nrow = n_dist, ncol = 1),
          gamman = c(
            rep(cfp_res_sample$gamman, n_old_dist),
            rep(pbo_res_sample$gamman, n_new_dist)
          )
        ) |>
        set_equilibrium(init_EIR = eir) |>
        run_aitn_simulation_v3() |>
        as.data.frame() |>
        mutate(
          itn_type            = "Pyr-PBO",
          resistance          = res,
          EIR                 = eir,
          rep                 = draw,
          pfpr2to10           = n_detect_730_3650 / n_730_3650,
          base_clin_inc_0_Inf = df_base |> pull(n_clin_inc_0_Inf),
          cfp_clin_inc_0_Inf  = df_cfp |> pull(n_clin_inc_0_Inf)
        )

      # ATN
      df_atn <- get_parameters(
        n_days      = n_steps,
        gamma_atn   = log(2) / atn_halflife,
        s_half_eip  = s_half_eip,
        nH_eip      = nH_eip,
        rho_frac    = rho_frac,
        s_half_pre  = s_half_pre,
        nH_pre      = nH_pre,
        B_max_post  = B_max_post,
        s_half_post = s_half_post,
        nH_post     = nH_post,
        Q0_atn      = new_covs,
        t0_atn      = new_distribution_times
      ) |>
        set_bednets(
          days      = distribution_times,
          coverages = all_covs,
          retention = retention_time,
          dn0 = matrix(
            c(
              rep(cfp_res_sample$dn0, n_old_dist),
              rep(0, n_new_dist)
            ),
            nrow = n_dist,
            ncol = 1
          ),
          rn = matrix(
            c(
              rep(cfp_res_sample$rn0, n_old_dist),
              rep(.24, n_new_dist)
            ),
            nrow = n_dist,
            ncol = 1
          ),
          rnm = matrix(rep(.24, n_dist), nrow = n_dist, ncol = 1),
          gamman = c(
            rep(cfp_res_sample$gamman, n_old_dist),
            rep(atn_halflife, n_new_dist)
          )
        ) |>
        set_equilibrium(init_EIR = eir) |>
        run_aitn_simulation_v3() |>
        as.data.frame() |>
        mutate(
          itn_type            = "ATN",
          resistance          = res,
          EIR                 = eir,
          rep                 = draw,
          pfpr2to10           = n_detect_730_3650 / n_730_3650,
          base_clin_inc_0_Inf = df_base |> pull(n_clin_inc_0_Inf),
          cfp_clin_inc_0_Inf  = df_cfp |> pull(n_clin_inc_0_Inf)
        )

      # ----- Pyr-ATN -----
      df_pyratn <- get_parameters(
        n_days      = n_steps,
        gamma_atn   = log(2) / atn_halflife,
        s_half_eip  = s_half_eip,
        nH_eip      = nH_eip,
        rho_frac    = rho_frac,
        s_half_pre  = s_half_pre,
        nH_pre      = nH_pre,
        B_max_post  = B_max_post,
        s_half_post = s_half_post,
        nH_post     = nH_post,
        Q0_atn      = new_covs,
        dn0_atn     = only_res_sample$dn0,
        t0_atn      = new_distribution_times
      ) |>
        set_bednets(
          days      = distribution_times,
          coverages = all_covs,
          retention = retention_time,
          dn0 = matrix(
            c(
              rep(cfp_res_sample$dn0, n_old_dist),
              rep(only_res_sample$dn0, n_new_dist)
            ),
            nrow = n_dist,
            ncol = 1
          ),
          rn = matrix(
            c(
              rep(cfp_res_sample$rn0, n_old_dist),
              rep(only_res_sample$rn0, n_new_dist)
            ),
            nrow = n_dist,
            ncol = 1
          ),
          rnm = matrix(rep(.24, n_dist), nrow = n_dist, ncol = 1),
          gamman = c(
            rep(cfp_res_sample$gamman, n_old_dist),
            rep(only_res_sample$gamman, n_new_dist)
          )
        ) |>
        set_equilibrium(init_EIR = eir)  |>
        run_aitn_simulation_v3() |>
        as.data.frame() |>
        mutate(
          itn_type            = "Pyr-ATN",
          resistance          = res,
          EIR                 = eir,
          rep                 = draw,
          pfpr2to10           = n_detect_730_3650 / n_730_3650,
          base_clin_inc_0_Inf  = df_base |> pull(n_clin_inc_0_Inf),
          cfp_clin_inc_0_Inf  = df_cfp |> pull(n_clin_inc_0_Inf)
        )



      # ----- Pyr-CFP-ATN -----
      df_cfpatn <- get_parameters(
        n_days      = n_steps,
        gamma_atn   = log(2) / atn_halflife,
        s_half_eip  = s_half_eip,
        nH_eip      = nH_eip,
        rho_frac    = rho_frac,
        s_half_pre  = s_half_pre,
        nH_pre      = nH_pre,
        B_max_post  = B_max_post,
        s_half_post = s_half_post,
        nH_post     = nH_post,
        Q0_atn      = new_covs,
        dn0_atn     = cfp_res_sample$dn0,
        t0_atn      = new_distribution_times
      ) |>
        set_bednets(
          days      = distribution_times,
          coverages = all_covs,
          retention = retention_time,
          dn0 = matrix(
            c(
              rep(cfp_res_sample$dn0, n_old_dist),
              rep(cfp_res_sample$dn0, n_new_dist)
            ),
            nrow = n_dist,
            ncol = 1
          ),
          rn = matrix(
            c(
              rep(cfp_res_sample$rn0, n_old_dist),
              rep(cfp_res_sample$rn0, n_new_dist)
            ),
            nrow = n_dist,
            ncol = 1
          ),
          rnm = matrix(rep(.24, n_dist), nrow = n_dist, ncol = 1),
          gamman = c(
            rep(cfp_res_sample$gamman, n_old_dist),
            rep(cfp_res_sample$gamman, n_new_dist)
          )
        ) |>
        set_equilibrium(init_EIR = eir)  |>
        run_aitn_simulation_v3() |>
        as.data.frame() |>
        mutate(
          itn_type            = "Pyr-CFP-ATN",
          resistance          = res,
          EIR                 = eir,
          rep                 = draw,
          pfpr2to10           = n_detect_730_3650 / n_730_3650,
          base_clin_inc_0_Inf  = df_base |> pull(n_clin_inc_0_Inf),
          cfp_clin_inc_0_Inf  = df_cfp |> pull(n_clin_inc_0_Inf)
        )

      all_runs[[length(all_runs) + 1]] <- df_only
      all_runs[[length(all_runs) + 1]] <- df_pbo
      all_runs[[length(all_runs) + 1]] <- df_cfp
      all_runs[[length(all_runs) + 1]] <- df_atn
      all_runs[[length(all_runs) + 1]] <- df_pyratn
      all_runs[[length(all_runs) + 1]] <- df_cfpatn

      cat("     Completed resistance =", res, "for EIR =", eir, "\n")

    }

    cat("Finished resistance scenarios for draw = ", draw , " and EIR =", eir, "\n")

  }

}

runtime <- proc.time() - ptm

cat("Runtime = ", unname(runtime[3]), " seconds")

# ------------------------------------------------
# Final combined dataframe
# ------------------------------------------------
df_full <- bind_rows(all_runs) |>
  mutate(
    avert_clin_inc_0_Inf = base_clin_inc_0_Inf - n_clin_inc_0_Inf,
    add_avert_clin_inc_0_Inf = cfp_clin_inc_0_Inf - n_clin_inc_0_Inf,
    itn_type = factor(
      itn_type, levels = c(
        "ATN", "Pyr-ATN", "Pyr-CFP-ATN", "Pyr-only", "Pyr-PBO", "Pyr-CFP"
      )
    ),
    year = (time - new_distribution_times[1]) / 365
  ) |>
  filter(year <= 6) |>
  filter(year >= 0)

saveRDS(df_full, file = "atn_sims_02JUN26_long.RDS")

#saveRDS(df_full, file = "atn_sims_MAY26_long.RDS")
#df_full <- readRDS(file = "atn_sims_MAY26_long.RDS")

df_full <- df_full |>
  mutate(
    itn_type = factor(
      itn_type, levels = c(
        #"ATN", "Pyr-ATN", "Pyr-CFP-ATN", "Pyr-only", "Pyr-PBO", "Pyr-CFP"
        "Pyr-only", "Pyr-PBO", "Pyr-CFP", "ATN", "Pyr-ATN", "Pyr-CFP-ATN"
      )
    )
  )

res_pct_labeller <- function(values) {
  vals_num <- as.numeric(values)
  paste0("Pyrethroid resistance = ", vals_num * 100, "%")
}


df_pfpr <- df_full |>
  group_by(itn_type, resistance, EIR, time, year) |>
  summarise(
    pfpr2to10_med = median(pfpr2to10),
    pfpr2to10_lo  = quantile(pfpr2to10, 0.025),
    pfpr2to10_hi  = quantile(pfpr2to10, 0.975),
    .groups = "drop"
  )

itn_cols <- c(
  "Pyr-only"    = "#0077BB",
  "Pyr-PBO"     = "#33BBEE",
  "Pyr-CFP"     = "#009988",
  "ATN"         = "#EE7733",
  "Pyr-ATN"     = "#CC3311",
  "Pyr-CFP-ATN" = "#EE3377"
)

itn_cols_dark <- c(
  "Pyr-only"    = "#005A8D",
  "Pyr-PBO"     = "#2486A8",
  "Pyr-CFP"     = "#006F63",
  "ATN"         = "#B85A26",
  "Pyr-ATN"     = "#992609",
  "Pyr-CFP-ATN" = "#B32659"
)

# itn_cols <- c(
#   "Pyr-only"    = "#0077BB",
#   "Pyr-PBO"     = "#33BBEE",
#   "Pyr-CFP"     = "#009988",
#   "ATN"         = "#EE3377",
#   "Pyr-ATN"     = "#EE7733",
#   "Pyr-CFP-ATN" = "#CC3311"
# )
#
# itn_cols_dark <- c(
#   "Pyr-only"    = "#005A8D",
#   "Pyr-PBO"     = "#2486A8",
#   "Pyr-CFP"     = "#006F63",
#   "ATN"         = "#B32659",
#   "Pyr-ATN"     = "#B85A26",
#   "Pyr-CFP-ATN" = "#992609"
# )

p_pfpr <- df_pfpr |>
  filter(itn_type %in% c("Pyr-only", "Pyr-CFP", "ATN", "Pyr-ATN")) |>
  ggplot(
    aes(
      x      = year,
      y      = pfpr2to10_med * 100,
      ymin   = pfpr2to10_lo * 100,
      ymax   = pfpr2to10_hi * 100,
      colour = itn_type,
      fill   = itn_type
    )
  ) +
  geom_ribbon(alpha = 0.15, colour = NA) +
  geom_line() +
  scale_color_manual(values = itn_cols) +
  scale_fill_manual(values = itn_cols) +
  facet_grid(
    EIR ~ resistance,
    labeller = labeller(
      EIR = eir_labels,
      resistance = res_pct_labeller
    ),
    scales = "free_y"
  ) +
  theme_minimal(base_size = 13) +
  labs(
    x      = "Years post distribution",
    y      = expression(italic(Pf) * PR[2-10] * " (%)"),
    colour = "",
    fill   = ""
  )

p_pfpr

ggsave("atn_pfpr_grid_fixed.png",   p_pfpr,
       width = 9, height = 7, dpi = 450, bg = "white")




df_clin <- df_full |>
  group_by(itn_type, resistance, EIR, rep) |>
  summarise(
    total_base = sum(base_clin_inc_0_Inf, na.rm = TRUE),
    total_averted = sum(avert_clin_inc_0_Inf, na.rm = TRUE),
    total_ref = sum(cfp_clin_inc_0_Inf, na.rm = TRUE),
    total_add_averted = sum(add_avert_clin_inc_0_Inf, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    avert_prop = total_averted / total_base,
    add_avert_prop = total_add_averted / total_ref
    ) |>
  group_by(itn_type, resistance, EIR) |>
  summarise(
    avert_med = median(total_averted),
    avert_lo  = quantile(total_averted, 0.025),
    avert_hi  = quantile(total_averted, 0.975),
    avert_prop_med = median(avert_prop),
    avert_prop_lo  = quantile(avert_prop, 0.025),
    avert_prop_hi  = quantile(avert_prop, 0.975),
    add_avert_med = median(total_add_averted),
    add_avert_lo  = quantile(total_add_averted, 0.025),
    add_avert_hi  = quantile(total_add_averted, 0.975),
    add_avert_prop_med = median(add_avert_prop),
    add_avert_prop_lo  = quantile(add_avert_prop, 0.025),
    add_avert_prop_hi  = quantile(add_avert_prop, 0.975),
    .groups = "drop"
  ) |>
  mutate(
    avert_text = sprintf(
      "%.1f%%\n(%.1f, %.1f)",
      avert_prop_med * -100,
      avert_prop_lo  * -100,
      avert_prop_hi  * -100
    ),
    add_avert_text = sprintf(
      "%.1f%%\n(%.1f, %.1f)",
      add_avert_prop_med * -100,
      add_avert_prop_lo  * -100,
      add_avert_prop_hi  * -100
    )
  )

p_clin_abs <- df_clin |>
  ggplot(
    aes(
      x       = itn_type,
      y       = avert_med / 100 / 3,
      ymin    = avert_lo / 100 / 3,
      ymax    = avert_hi / 100 / 3,
      fill    = itn_type,
      colour  = itn_type
    )
  ) +
  geom_col(width = 0.75, colour = NA, alpha = 0.6) +
  geom_errorbar(width = 0.18, linewidth = 0.5) +
  scale_fill_manual(values = itn_cols) +
  scale_colour_manual(values = itn_cols_dark) +
  facet_grid(
    EIR ~ resistance,
    labeller = labeller(
      EIR = eir_labels,
      resistance = res_pct_labeller
    )#,
    #scales = "free_y"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x    = "",
    y    = "Annual average clinical cases averted vs no future distribution per 1,000",
    fill = ""
  )

p_clin_abs
# ggsave("atn_clin_abs_grid_fixed.png",   p_clin_abs,
#        width = 9, height = 7, dpi = 450, bg = "white")

p_clin_rel <- df_clin |>
  ggplot(
    aes(
      x       = itn_type,
      y       = avert_prop_med * 100,
      ymin    = avert_prop_lo * 100,
      ymax    = avert_prop_hi * 100,
      fill    = itn_type,
      colour  = itn_type
    )
  ) +
  geom_col(width = 0.75, colour = NA, alpha = 0.6) +
  geom_errorbar(width = 0.18, linewidth = 0.5) +
  scale_fill_manual(values = itn_cols) +
  scale_colour_manual(values = itn_cols_dark) +
  facet_grid(
    EIR ~ resistance,
    labeller = labeller(
      EIR = eir_labels,
      resistance = res_pct_labeller
    )#,
    #scales = "free_y"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x    = "",
    y    = "Reduction in clinical cases vs no future distribution (%)",
    fill = ""
  )

p_clin_rel
# ggsave("atn_clin_rel_grid_fixed.png",   p_clin_rel,
#        width = 9, height = 7, dpi = 450, bg = "white")



p_clin_add_abs <- df_clin |>
  filter(itn_type %in% c("ATN", "Pyr-ATN", "Pyr-CFP-ATN")) |>
  ggplot(
    aes(
      x       = itn_type,
      y       = add_avert_med / 100 / 3,
      ymin    = add_avert_lo / 100 / 3,
      ymax    = add_avert_hi / 100 / 3,
      fill    = itn_type,
      colour  = itn_type
    )
  ) +
  geom_col(width = 0.75, colour = NA, alpha = 0.6) +
  geom_errorbar(width = 0.18, linewidth = 0.5) +
  scale_fill_manual(values = itn_cols) +
  scale_colour_manual(values = itn_cols_dark) +
  scale_y_continuous(
    breaks = seq(-200, 400, by = 50),       # major ticks
    minor_breaks = seq(-200, 400, by = 25),  # minor ticks
    limits = c(-105, 305)
  ) +
  facet_grid(
    EIR ~ resistance,
    labeller = labeller(
      EIR = eir_labels,
      resistance = res_pct_labeller
    )#,
    #scales = "free_y"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x    = "",
    y    = "Annual average clinical cases averted vs Pyr-CFP per 1,000",
    fill = ""
  )

p_clin_add_abs
# ggsave("atn_clin_add_abs_grid_fixed_v2.png",   p_clin_add_abs,
#        width = 9, height = 7, dpi = 450, bg = "white")

p_clin_add_rel <- df_clin |>
  filter(itn_type %in% c("ATN", "Pyr-ATN", "Pyr-CFP-ATN")) |>
  ggplot(
    aes(
      x       = itn_type,
      y       = add_avert_prop_med * 100,
      ymin    = add_avert_prop_lo * 100,
      ymax    = add_avert_prop_hi * 100,
      fill    = itn_type,
      colour  = itn_type
    )
  ) +
  geom_col(width = 0.75, colour = NA, alpha = 0.6) +
  geom_errorbar(width = 0.18, linewidth = 0.5) +
  scale_fill_manual(values = itn_cols) +
  scale_colour_manual(values = itn_cols_dark) +
  scale_y_continuous(
   breaks = seq(-200, 100, by = 50),       # major ticks
   minor_breaks = seq(-200, 100, by = 25),  # minor ticks
   limits = c(-155, 105)
  ) +
  facet_grid(
    EIR ~ resistance,
    labeller = labeller(
      EIR = eir_labels,
      resistance = res_pct_labeller
    ),
    scales = "free_y"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x    = "",
    y    = "Reduction in clinical cases vs Pyr-CFP (%)",
    fill = ""
  )

p_clin_add_rel
# ggsave("atn_clin_add_rel_grid_fixed_v2.png",   p_clin_add_rel,
#        width = 9, height = 7, dpi = 450, bg = "white")


p_clin_add_abs_rel <- df_clin |>
  filter(itn_type %in% c("ATN", "Pyr-ATN", "Pyr-only")) |>
  ggplot(
    aes(
      x       = itn_type,
      y       = add_avert_med / 100 / 6,
      ymin    = add_avert_lo / 100 / 6,
      ymax    = add_avert_hi / 100 / 6,
      fill    = itn_type,
      colour  = itn_type
    )
  ) +
  geom_text(
    aes(
      y = pmax(0, add_avert_hi / 100 / 6) + 50,
      label = add_avert_text
    ),
    #vjust = -0.3,
    size = 2
  ) +
  geom_col(width = 0.75, colour = NA, alpha = 0.6) +
  geom_errorbar(width = 0.18, linewidth = 0.5) +
  scale_fill_manual(values = itn_cols) +
  scale_colour_manual(values = itn_cols_dark) +
  scale_y_continuous(
    breaks = seq(-600, 1000, by = 100),       # major ticks
    minor_breaks = seq(-600, 1000, by = 25),  # minor ticks
    limits = c(-250, 300)
  ) +
  facet_grid(
    EIR ~ resistance,
    labeller = labeller(
      EIR = eir_labels,
      resistance = res_pct_labeller
    ),
    scales = "free_y"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x    = "",
    y    = "Annual average clinical cases averted vs Pyr-CFP per 1,000",
    fill = ""
  )# +
  #ylim(-160,110)

p_clin_add_abs_rel
ggsave("atn_clin_add_abs_rel_grid_v2.png",   p_clin_add_abs_rel,
       width = 9, height = 7, dpi = 450, bg = "white")







p_clin_abs_rel <- df_clin |>
  filter(itn_type %in% c("ATN", "Pyr-ATN", "Pyr-only", "Pyr-CFP")) |>
  ggplot(
    aes(
      x       = itn_type,
      y       = avert_med / 100 / 6,
      ymin    = avert_lo / 100 / 6,
      ymax    = avert_hi / 100 / 6,
      fill    = itn_type,
      colour  = itn_type
    )
  ) +
  # geom_text(
  #   aes(
  #     y = pmax(0, avert_hi / 100 / 3),# + 15,
  #     label = avert_text
  #   ),
  #   vjust = -0.4,
  #   size = 2
  # ) +
  geom_text(
    aes(
      y = pmax(0, avert_hi / 100 / 6) + 50,
      label = avert_text
    ),
    #vjust = -0.4,
    size = 2
  ) +
  geom_col(width = 0.75, colour = NA, alpha = 0.6) +
  geom_errorbar(width = 0.18, linewidth = 0.5) +
  scale_fill_manual(values = itn_cols) +
  scale_colour_manual(values = itn_cols_dark) +
  scale_y_continuous(
    breaks = seq(-500, 1500, by = 100),       # major ticks
    minor_breaks = seq(-500, 1500, by = 25),  # minor ticks
    limits = c(0, 650)
  ) +
  facet_grid(
    EIR ~ resistance,
    labeller = labeller(
      EIR = eir_labels,
      resistance = res_pct_labeller
    ),
    scales = "free_y"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x    = "",
    y    = "Annual average clinical cases averted vs no future distribution per 1,000",
    fill = ""
  )# +
#ylim(-160,110)

p_clin_abs_rel
ggsave("p_clin_abs_rel.png",   p_clin_abs_rel,
      width = 9, height = 7, dpi = 450, bg = "white")


