library(dplyr)
library(tidyr)
library(ggplot2)
library(viridisLite)
library(viridis)

read_itn_pars <- TRUE

n_steps <- 365 * 9

# ------------------------------------------------
# User-editable settings
# ------------------------------------------------
eir_levels <- c(0.5, 4, 30)
res_levels <- c(0, 0.5, 0.8, 0.9, 1) #c(0, 0.5, 0.8, 0.9, 1)
# ------------------------------------------------

if (read_itn_pars) {
  only_pars <- readRDS("./dev/net_params/dat_res_only.rds")
  pbo_pars  <- readRDS("./dev/net_params/dat_res_pbo.rds")
  cfp_pars  <- readRDS("./dev/net_params/dat_res_cfp.rds")
}

# ------------------------------------------------
# Median parameter sets (precomputed once)
# ------------------------------------------------
only_med <- only_pars |> group_by(resistance) |> summarise(
  dn0 = median(dn0), rn0 = median(rn0), gamman = median(gamman), .groups = "drop"
)

pbo_med <- pbo_pars |> group_by(resistance) |> summarise(
  dn0 = median(dn0), rn0 = median(rn0), gamman = median(gamman), .groups = "drop"
)

cfp_med <- cfp_pars |> group_by(resistance) |> summarise(
  dn0 = median(dn0), rn0 = median(rn0), gamman = median(gamman), .groups = "drop"
)

# ------------------------------------------------
# Storage
# ------------------------------------------------
all_runs <- list()

# ------------------------------------------------
# Loop over EIR and resistance levels
# ------------------------------------------------
for (eir in eir_levels) {

  # ------------------------------------------------
  # Baseline scenario
  # ------------------------------------------------
  params_base <- get_parameters(n_days = n_steps,
                                spor_len = 12,
                                Q0_atn = 0) |>
    set_equilibrium(init_EIR = eir)
  params_base$max_atn_cov <- 0
  df_base <- run_atn_simulation(params_base) |>
    as.data.frame() |>
    mutate(
      model     = "basline",
      itn_type  = "none",
      pfpr2to10 = n_detect_730_3650 / n_730_3650,
      EIR       = eir
    )

  # ----- ATN (run once per EIR) -----
  params_atn <- get_parameters(n_days = n_steps,
                               spor_len = 12,
                               gamma_atn = (log(2) / (2.64 * 365))) |>
    set_bednets(days = 1095,
                coverages = 0.9,
                gamman =   2.64 * 365,
                retention =  588,
                dn0 = 0,
                rn  = 0.24,#56,
                rnm = 0.24) |>
    set_equilibrium(init_EIR = eir)

  params_atn$max_atn_cov <- 0

  df_atn <- run_atn_simulation(params_atn) |>
    as.data.frame() |>
    mutate(
      model     = "ATN",
      itn_type  = "ATN",
      pfpr2to10 = n_detect_730_3650 / n_730_3650,
      EIR       = eir
    )
  df_atn$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

  # Repeat ATN for each resistance
  df_atn_rep <- bind_rows(lapply(res_levels, function(r) {
    df_atn |> mutate(resistance = r)
  }))

  # ----- ITN loop -----
  for (res in res_levels) {

    # ONLY
    pars_only <- filter(only_med, resistance == res)
    df_only <- run_atn_simulation(
      get_parameters(n_days = n_steps, spor_len = 12, Q0_atn = 0) |>
        set_bednets(
          days = 1095, coverages = 0.9,
          gamman = pars_only$gamman * 365,
          retention = 588, dn0 = pars_only$dn0,
          rn = pars_only$rn0, rnm = 0.24
        ) |>
        set_equilibrium(init_EIR = eir)
    ) |>
      as.data.frame() |>
      mutate(
        model = "ITN", itn_type = "Pyr-only",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_only$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

    # PBO
    pars_pbo <- filter(pbo_med, resistance == res)
    df_pbo <- run_atn_simulation(
      get_parameters(n_days = n_steps, spor_len = 12, Q0_atn = 0) |>
        set_bednets(
          days = 1095, coverages = 0.9,
          gamman = pars_pbo$gamman * 365,
          retention = 588, dn0 = pars_pbo$dn0,
          rn = pars_pbo$rn0, rnm = 0.24
        ) |>
        set_equilibrium(init_EIR = eir)
    ) |>
      as.data.frame() |>
      mutate(
        model = "ITN", itn_type = "Pyr-PBO",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_pbo$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

    # CFP
    pars_cfp <- filter(cfp_med, resistance == res)
    df_cfp <- run_atn_simulation(
      get_parameters(n_days = n_steps, spor_len = 12, Q0_atn = 0) |>
        set_bednets(
          days = 1095, coverages = 0.9,
          gamman = pars_cfp$gamman * 365,
          retention = 588, dn0 = pars_cfp$dn0,
          rn = pars_cfp$rn0, rnm = 0.24
        ) |>
        set_equilibrium(init_EIR = eir)
    ) |>
      as.data.frame() |>
      mutate(
        model = "ITN", itn_type = "Pyr-CFP",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_cfp$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

    # ----- Pyr-ATN -----
    #pars_pyratn1 <- filter(only_med, resistance == res)
    df_pyratn1 <- run_aitn_simulation(
      get_parameters(n_days = n_steps, spor_len = 12,
                     gamma_atn = (log(2) / (2.64 * 365)),
                     dn0_atn = pars_only$dn0) |>
        set_bednets(
          days = 1095, coverages = 0.9,
          gamman = pars_only$gamman * 365,
          retention = 588, dn0 = 0, #pars_only$dn0,
          rn = pars_only$rn0, rnm = 0.24
        ) |>
        set_equilibrium(init_EIR = eir)
    ) |>
      as.data.frame() |>
      mutate(
        model = "ITN", itn_type = "Pyr-ATNv1",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_pyratn1$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf
    #
    # pars_pyratn <- get_parameters(n_days = n_steps,
    #                                spor_len = 12,
    #                                gamma_atn = (log(2) / (2.64 * 365))) |>
    #   set_bednets(
    #     days = 1095, coverages = 0.9,
    #     gamman = pars_only$gamman * 365,
    #     retention = 588, dn0 = pars_only$dn0,
    #     rn = pars_only$rn0, rnm = 0.24
    #   ) |>
    #   set_equilibrium(init_EIR = eir)
    #
    # pars_pyratn$max_atn_cov <- 0
    #
    # df_pyratn <- run_atn_simulation(pars_pyratn) |>
    #   as.data.frame() |>
    #   mutate(
    #     model     = "ATN-only",
    #     itn_type  = "Pyr-ATN",
    #     pfpr2to10 = n_detect_730_3650 / n_730_3650,
    #     EIR       = eir
    #   )

    # ----- Pyr-ATN -----
    #pars_indeppyratn <- filter(only_med, resistance == res)
    df_pyratn0 <- run_atn_simulation(
      get_parameters(n_days = n_steps, spor_len = 12,
                     gamma_atn = (log(2) / (2.64 * 365))) |>
        set_bednets(
          days = 1095, coverages = 0.9,
          gamman = pars_only$gamman * 365,
          retention = 588, dn0 = pars_only$dn0,
          rn = pars_only$rn0, rnm = 0.24
        ) |>
        set_equilibrium(init_EIR = eir)
    ) |>
      as.data.frame() |>
      mutate(
        model = "ITN", itn_type = "Pyr-ATNv0",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_pyratn0$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

    all_runs[[length(all_runs) + 1]] <- df_only
    all_runs[[length(all_runs) + 1]] <- df_pbo
    all_runs[[length(all_runs) + 1]] <- df_cfp
    all_runs[[length(all_runs) + 1]] <- df_pyratn0
    all_runs[[length(all_runs) + 1]] <- df_pyratn1
  }

  all_runs[[length(all_runs) + 1]] <- df_atn_rep
}

# ------------------------------------------------
# Final combined dataframe
# ------------------------------------------------
df_plot <- bind_rows(all_runs)
df_plot$avert_clin_inc_0_Inf <- df_plot$base_clin_inc_0_Inf - df_plot$n_clin_inc_0_Inf


# ------------------------------------------------
# Facet plot: resistance × EIR
# ------------------------------------------------
ggplot(df_plot, aes(x = time, y = pfpr2to10, color = itn_type)) +
  geom_line(linewidth = 0.8) +
  facet_grid(EIR ~ resistance) + #, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    x = "Time (days)",
    y = "PfPR2–10",
    color = "Type",
    title = "ATN vs ITN performance across resistance and EIR levels"
  )

eir_labels <- c(
  `0.5`   = "Low transmission",
  `4`  = "Moderate transmission",
  `30` = "High transmission"
)

res_pct_labeller <- function(values) {
  vals_num <- as.numeric(values)
  paste0("Pyr res = ", vals_num * 100, "%")
}

df_plot |>
  mutate(
    itn_type = factor(
      itn_type,
      levels = c("Pyr-only", "Pyr-PBO", "Pyr-CFP", "Pyr-ATNv0", "Pyr-ATNv1", "ATN")
    ),
    year = (time - 1095) / 365
  ) |>
  filter(year >= -1) |>
  ggplot(aes(x = year, y = pfpr2to10, color = itn_type)) +
  geom_line() +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, direction = -1) +
  #scale_color_viridis_d(option = "magma", end = 0.8, direction = 1) +
  facet_grid(EIR ~ resistance,
             labeller = labeller(
               EIR = eir_labels,
               resistance = res_pct_labeller
             )
  ) + #, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    x = "Years since distribution",
    y = "PfPR2–10",
    color = "Net type"
    #title = "Pyrethroid resistance (columns) vs baseline EIR (rows)"
  )



df_plot |>
  mutate(
    itn_type = factor(
      itn_type,
      levels = c("Pyr-only", "Pyr-PBO", "Pyr-CFP", "Pyr-ATNv0", "Pyr-ATNv1", "ATN")
    ),
    year = (time - 1095) / 365
  ) |>
  filter(year >= -1) |>
  ggplot(aes(x = year, y = avert_clin_inc_0_Inf, color = itn_type)) +
  geom_line(linewidth = 0.8) +
  scale_color_viridis_d(option = "turbo", begin = 0, end = 0.5, direction = -1) +
  #scale_color_viridis_d(option = "plasma", end = 0.95, direction = -1) +
  facet_grid(EIR ~ resistance,
             labeller = labeller(
               EIR = eir_labels,
               resistance = res_pct_labeller
             )
  ) + #, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    x = "Years since distribution",
    y = "Clinical cases averted",
    color = "Net type"
    #title = "Pyrethroid resistance (columns) vs baseline EIR (rows)"
  )

df_summary <- df_plot |>
  mutate(
    itn_type = factor(
      itn_type,
      levels = c("Pyr-only", "Pyr-PBO", "Pyr-CFP", "Pyr-ATNv0", "Pyr-ATNv1", "ATN")
    ),
    year = (time - 1095) / 365
  ) |>
  filter(year <= 3) |>
  filter(year >= 0) |>
  group_by(itn_type, resistance, EIR) |>
  summarise(
    total_base = sum(base_clin_inc_0_Inf, na.rm = TRUE),
    total_averted = sum(avert_clin_inc_0_Inf, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    avert_prop = total_averted / total_base
  )

df_summary |>
  ggplot(aes(x = itn_type, y = total_averted/100, fill = itn_type)) +
  geom_col() +
  geom_text(
    aes(label = scales::percent(avert_prop, accuracy = 0.1)),
    vjust = -0.5
  ) +
  #scale_fill_viridis_d(option = "turbo", begin = 0, end = 0.5, direction = -1) +
  facet_grid(EIR ~ resistance,
             labeller = labeller(
               EIR = eir_labels,
               resistance = res_pct_labeller
             )
  ) + #, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    y = "Clinical cases averted per 1,000 over 3 years",
    fill = "Net type"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  ylim(0, 1200)

