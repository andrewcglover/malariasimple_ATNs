library(dplyr)
library(tidyr)
library(ggplot2)
library(viridisLite)
library(viridis)
library(RColorBrewer)

read_itn_pars <- TRUE

n_steps <- 365 * 9

# ------------------------------------------------
# User-editable settings
# ------------------------------------------------
eir_levels <- c(0.7, 4.2, 26)
res_levels <- c(0.5, 0.9, 1)
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
  dn0 = median(dn0), rn0 = median(rn0), gamman = median(gamman),
  .groups = "drop"
)

pbo_med <- pbo_pars |> group_by(resistance) |> summarise(
  dn0 = median(dn0), rn0 = median(rn0), gamman = median(gamman),
  .groups = "drop"
)

cfp_med <- cfp_pars |> group_by(resistance) |> summarise(
  dn0 = median(dn0), rn0 = median(rn0), gamman = median(gamman),
  .groups = "drop"
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
  df_base <- get_parameters(
    n_days = n_steps,
    spor_len = 12,
    Q0_atn = 0
  ) |>
    set_equilibrium(init_EIR = eir) |>
    run_atn_simulation() |>
    as.data.frame() |>
    mutate(
      itn_type  = "none",
      pfpr2to10 = n_detect_730_3650 / n_730_3650,
      EIR       = eir
    )

  # ----- ATN (run once per EIR) -----
  # optimistic version (default of xi = 0.3)
  df_atn <- get_parameters(
    n_days    = n_steps,
    spor_len  = 12,
    gamma_atn = (log(2) / (2.64 * 365))
  ) |>
    set_bednets(
      days      = 1095,
      coverages = 0.9,
      gamman    = 2.64 * 365,
      retention = 588,
      dn0       = 0,
      rn        = 0.24, #56,
      rnm       = 0.24
    ) |>
    set_equilibrium(init_EIR = eir) |>
    run_atn_simulation() |>
    as.data.frame() |>
    mutate(
      itn_type            = "ATN",
      atn_action          = "Both",
      inhib_action        = "Long",
      pfpr2to10           = n_detect_730_3650 / n_730_3650,
      EIR                 = eir,
      base_clin_inc_0_Inf = df_base |> pull(n_clin_inc_0_Inf)
    )
  # Repeat ATN for each resistance
  df_atn_rep <- bind_rows(
    lapply(
      res_levels, function(r) {
        df_atn |> mutate(resistance = r)
      }
    )
  )

  # pessimistic version (default of xi = 0.5)
  df_atn_fasti <- get_parameters(
    n_days = n_steps,
    spor_len = 12,
    gamma_atn = (log(2) / (2.64 * 365)),
    xi = 0.6
  ) |>
    set_bednets(days = 1095,
                coverages = 0.9,
                gamman =   2.64 * 365,
                retention =  588,
                dn0 = 0,
                rn  = 0.24, #56,
                rnm = 0.24) |>
    set_equilibrium(init_EIR = eir) |>
    run_atn_simulation() |>
    as.data.frame() |>
    mutate(
      itn_type            = "ATN",
      atn_action          = "Both",
      inhib_action        = "Short",
      pfpr2to10 = n_detect_730_3650 / n_730_3650,
      EIR       = eir,
      base_clin_inc_0_Inf = df_base |> pull(n_clin_inc_0_Inf)
    )
  # Repeat ATN for each resistance
  df_atn_fasti_rep <- bind_rows(
    lapply(
      res_levels, function(r) {
        df_atn_fasti |> mutate(resistance = r)
      }
    )
  )

  # pessimistic version (default of xi = 0.5)
  df_atn_vfasti <- get_parameters(
    n_days = n_steps,
    spor_len = 12,
    gamma_atn = (log(2) / (2.64 * 365)),
    xi = 1e3
  ) |>
    set_bednets(days = 1095,
                coverages = 0.9,
                gamman =   2.64 * 365,
                retention =  588,
                dn0 = 0,
                rn  = 0.24, #56,
                rnm = 0.24) |>
    set_equilibrium(init_EIR = eir) |>
    run_atn_simulation() |>
    as.data.frame() |>
    mutate(
      itn_type            = "ATN",
      atn_action          = "Both",
      inhib_action        = "Concurrent",
      pfpr2to10 = n_detect_730_3650 / n_730_3650,
      EIR       = eir,
      base_clin_inc_0_Inf = df_base |> pull(n_clin_inc_0_Inf)
    )
  # Repeat ATN for each resistance
  df_atn_vfasti_rep <- bind_rows(
    lapply(
      res_levels, function(r) {
        df_atn_vfasti |> mutate(resistance = r)
      }
    )
  )

  # ATN EIP lengthening only
  df_atn_eip <- get_parameters(
    n_days      = n_steps,
    spor_len    = 12,
    gamma_atn   = (log(2) / (2.64 * 365)),
    Lambda00sf  = 1
  ) |>
    set_bednets(
      days = 1095,
      coverages = 0.9,
      gamman =   2.64 * 365,
      retention =  588,
      dn0 = 0,
      rn  = 0.24, #56,
      rnm = 0.24
    ) |>
    set_equilibrium(init_EIR = eir) |>
    run_atn_simulation() |>
    as.data.frame() |>
    mutate(
      itn_type            = "ATN",
      atn_action          = "EIP",
      inhib_action        = "Long",
      pfpr2to10 = n_detect_730_3650 / n_730_3650,
      EIR       = eir,
      base_clin_inc_0_Inf = df_base |> pull(n_clin_inc_0_Inf)
    )
  df_atn_eip_rep <- bind_rows(
    lapply(
      res_levels, function(r) {
        df_atn_eip |> mutate(resistance = r)
      }
    )
  )

  # ATN blocking only
  df_atn_block <- get_parameters(
    n_days = n_steps,
    spor_len = 12,
    gamma_atn = (log(2) / (2.64 * 365)),
    rho00 = 1
  ) |>
    set_bednets(
      days = 1095,
      coverages = 0.9,
      gamman =   2.64 * 365,
      retention =  588,
      dn0 = 0,
      rn  = 0.24, #56,
      rnm = 0.24
    ) |>
    set_equilibrium(init_EIR = eir) |>
    run_atn_simulation() |>
    as.data.frame() |>
    mutate(
      itn_type            = "ATN",
      atn_action          = "Inhibitory",
      inhib_action        = "Long",
      pfpr2to10           = n_detect_730_3650 / n_730_3650,
      EIR                 = eir,
      base_clin_inc_0_Inf = df_base |> pull(n_clin_inc_0_Inf)
    )
  df_atn_block_rep <- bind_rows(
    lapply(
      res_levels, function(r) {
        df_atn_block |> mutate(resistance = r)
      }
    )
  )

  # ATN blocking only
  df_atn_block_fasti <- get_parameters(
    n_days = n_steps,
    spor_len = 12,
    gamma_atn = (log(2) / (2.64 * 365)),
    rho00 = 1,
    xi = 0.6
  ) |>
    set_bednets(
      days = 1095,
      coverages = 0.9,
      gamman =   2.64 * 365,
      retention =  588,
      dn0 = 0,
      rn  = 0.24, #56,
      rnm = 0.24
    ) |>
    set_equilibrium(init_EIR = eir) |>
    run_atn_simulation() |>
    as.data.frame() |>
    mutate(
      itn_type            = "ATN",
      atn_action          = "Inhibitory",
      inhib_action        = "Short",
      pfpr2to10           = n_detect_730_3650 / n_730_3650,
      EIR                 = eir,
      base_clin_inc_0_Inf = df_base |> pull(n_clin_inc_0_Inf)
    )
  df_atn_block_fasti_rep <- bind_rows(
    lapply(
      res_levels, function(r) {
        df_atn_block_fasti |> mutate(resistance = r)
      }
    )
  )

  df_atn_block_vfasti <- get_parameters(
    n_days = n_steps,
    spor_len = 12,
    gamma_atn = (log(2) / (2.64 * 365)),
    rho00 = 1,
    xi = 1e3
  ) |>
    set_bednets(
      days = 1095,
      coverages = 0.9,
      gamman =   2.64 * 365,
      retention =  588,
      dn0 = 0,
      rn  = 0.24, #56,
      rnm = 0.24
    ) |>
    set_equilibrium(init_EIR = eir) |>
    run_atn_simulation() |>
    as.data.frame() |>
    mutate(
      itn_type            = "ATN",
      atn_action          = "Inhibitory",
      inhib_action        = "Concurrent",
      pfpr2to10           = n_detect_730_3650 / n_730_3650,
      EIR                 = eir,
      base_clin_inc_0_Inf = df_base |> pull(n_clin_inc_0_Inf)
    )
  df_atn_block_vfasti_rep <- bind_rows(
    lapply(
      res_levels, function(r) {
        df_atn_block_vfasti |> mutate(resistance = r)
      }
    )
  )

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
        itn_type            = "Pyr-only",
        atn_action          = "Both",
        inhib_action        = "Long",
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
        itn_type            = "Pyr-PBO",
        atn_action          = "Both",
        inhib_action        = "Long",
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
        itn_type            = "Pyr-CFP",
        atn_action          = "Both",
        inhib_action        = "Long",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_cfp$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

    # ----- Pyr-ATN -----
    #pars_pyratn1 <- filter(only_med, resistance == res)
    df_pyratn <- run_aitn_simulation(
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
        itn_type            = "Pyr-ATN",
        atn_action          = "Both",
        inhib_action        = "Long",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_pyratn$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

    df_pyratn_fasti <- run_aitn_simulation(
      get_parameters(n_days = n_steps, spor_len = 12,
                     gamma_atn = (log(2) / (2.64 * 365)),
                     dn0_atn = pars_only$dn0, xi = 0.6) |>
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
        itn_type            = "Pyr-ATN",
        atn_action          = "Both",
        inhib_action        = "Short",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_pyratn_fasti$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

    df_pyratn_vfasti <- run_aitn_simulation(
      get_parameters(n_days = n_steps, spor_len = 12,
                     gamma_atn = (log(2) / (2.64 * 365)),
                     dn0_atn = pars_only$dn0, xi = 1e3) |>
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
        itn_type            = "Pyr-ATN",
        atn_action          = "Both",
        inhib_action        = "Concurrent",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_pyratn_vfasti$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

    df_pyratn_eip <- run_aitn_simulation(
      get_parameters(n_days = n_steps, spor_len = 12,
                     gamma_atn = (log(2) / (2.64 * 365)),
                     dn0_atn = pars_only$dn0,
                     Lambda00sf = 1) |>
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
        itn_type            = "Pyr-ATN",
        atn_action          = "EIP",
        inhib_action        = "Long",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_pyratn_eip$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

    df_pyratn_block <- run_aitn_simulation(
      get_parameters(n_days = n_steps, spor_len = 12,
                     gamma_atn = (log(2) / (2.64 * 365)),
                     dn0_atn = pars_only$dn0,
                     rho00 = 1) |>
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
        itn_type            = "Pyr-ATN",
        atn_action          = "Inhibitory",
        inhib_action        = "Long",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_pyratn_block$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

    df_pyratn_block_fasti <- run_aitn_simulation(
      get_parameters(n_days = n_steps, spor_len = 12,
                     gamma_atn = (log(2) / (2.64 * 365)),
                     dn0_atn = pars_only$dn0,
                     rho00 = 1,
                     xi = 0.6) |>
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
        itn_type            = "Pyr-ATN",
        atn_action          = "Inhibitory",
        inhib_action        = "Short",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_pyratn_block_fasti$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

    df_pyratn_block_vfasti <- run_aitn_simulation(
      get_parameters(n_days = n_steps, spor_len = 12,
                     gamma_atn = (log(2) / (2.64 * 365)),
                     dn0_atn = pars_only$dn0,
                     rho00 = 1,
                     xi = 1e3) |>
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
        itn_type            = "Pyr-ATN",
        atn_action          = "Inhibitory",
        inhib_action        = "Concurrent",
        resistance = res, EIR = eir,
        pfpr2to10 = n_detect_730_3650 / n_730_3650
      )
    df_pyratn_block_vfasti$base_clin_inc_0_Inf <- df_base$n_clin_inc_0_Inf

    all_runs[[length(all_runs) + 1]] <- df_only
    all_runs[[length(all_runs) + 1]] <- df_pbo
    all_runs[[length(all_runs) + 1]] <- df_cfp
    all_runs[[length(all_runs) + 1]] <- df_pyratn
    all_runs[[length(all_runs) + 1]] <- df_pyratn_fasti
    all_runs[[length(all_runs) + 1]] <- df_pyratn_vfasti
    all_runs[[length(all_runs) + 1]] <- df_pyratn_eip
    all_runs[[length(all_runs) + 1]] <- df_pyratn_block
    all_runs[[length(all_runs) + 1]] <- df_pyratn_block_fasti
    all_runs[[length(all_runs) + 1]] <- df_pyratn_block_vfasti
  }

  all_runs[[length(all_runs) + 1]] <- df_atn_rep
  all_runs[[length(all_runs) + 1]] <- df_atn_fasti_rep
  all_runs[[length(all_runs) + 1]] <- df_atn_vfasti_rep
  all_runs[[length(all_runs) + 1]] <- df_atn_eip_rep
  all_runs[[length(all_runs) + 1]] <- df_atn_block_rep
  all_runs[[length(all_runs) + 1]] <- df_atn_block_fasti_rep
  all_runs[[length(all_runs) + 1]] <- df_atn_block_vfasti_rep
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

label_names <- paste(
  c("Low", "Moderate", "High", "Very high", "Extreme")[seq_along(eir_levels)],
  "transmission"
)

eir_labels <- setNames(label_names, eir_levels)

res_pct_labeller <- function(values) {
  vals_num <- as.numeric(values)
  paste0("Pyrethroid resistance = ", vals_num * 100, "%")
}

df_plot |>
  mutate(
    itn_type = factor(
      itn_type,
      levels = c("Pyr-only", "Pyr-PBO", "Pyr-CFP", "Pyr-ATN", "ATN")
    ),
    inhib_action = factor(
      inhib_action,
      levels = c("Concurrent", "Short", "Long")
    ),
    year = (time - 1095) / 365
  ) |>
  filter(year >= -1) |>
  filter(atn_action == "Both") |>
  filter(inhib_action == "Long") |>
  #filter(itn_type != "ATN") |>
  #filter(itn_type == "Pyr-ATN") |>
  ggplot(
    aes(
      x = year, y = pfpr2to10,
      color = itn_type,
      #linetype = atn_action,
      #alpha = inhib_action
      )
    ) +
  #geom_line(colour = "#e7298a") +
  geom_line() +
  scale_colour_brewer(palette = "Dark2") +
  #scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, direction = -1) +
  #scale_color_viridis_d(option = "magma", end = 0.8, direction = 1) +
  scale_alpha_discrete(range = c(0.4, 1)) +
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
    color = "Net type",
    #linetype = "ATN action",
    #alpha = "Inhibition\nduration"
    #title = "Pyrethroid resistance (columns) vs baseline EIR (rows)"
  )




df_plot |>
  mutate(
    itn_type = factor(
      itn_type,
      levels = c("Pyr-only", "Pyr-PBO", "Pyr-CFP", "Pyr-ATN", "ATN")
    ),
    inhib_action = factor(
      inhib_action,
      levels = c("Concurrent", "Short", "Long")
    ),
    year = (time - 1095) / 365
  ) |>
  filter(year >= -1) |>
  #filter(atn_action == "Both") |>
  #filter(inhib_action == "Long") |>
  filter(itn_type %in% c("Pyr-ATN", "ATN")) |>
  filter(itn_type == "Pyr-ATN") |>
  ggplot(
    aes(
      x = year, y = pfpr2to10,
      color = itn_type,
      linetype = atn_action,
      alpha = inhib_action
    )
  ) +
  #geom_line(colour = "#e7298a") +
  geom_line() +
  #scale_colour_brewer(palette = "Dark2") +
  scale_color_manual(values = c("Pyr-ATN" = "#e7298a", "ATN" = "#66a61e")) +
  #scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, direction = -1) +
  #scale_color_viridis_d(option = "magma", end = 0.8, direction = 1) +
  scale_alpha_discrete(range = c(0.4, 1)) +
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
    color = "Net type",
    linetype = "ATN action",
    alpha = "Inhibition\nduration"
    #title = "Pyrethroid resistance (columns) vs baseline EIR (rows)"
  )




df_plot |>
  mutate(
    itn_type = factor(
      itn_type,
      levels = c("Pyr-only", "Pyr-PBO", "Pyr-CFP", "Pyr-ATN", "ATN")
    ),
    inhib_action = factor(
      inhib_action,
      levels = c("Concurrent", "Short", "Long")
    ),
    year = (time - 1095) / 365
  ) |>
  filter(year >= -1) |>
  #filter(atn_action == "Both") |>
  #filter(inhib_action == "Long") |>
  filter(itn_type %in% c("Pyr-ATN", "ATN")) |>
  filter(itn_type == "ATN") |>
  ggplot(
    aes(
      x = year, y = pfpr2to10,
      color = itn_type,
      linetype = atn_action,
      alpha = inhib_action
    )
  ) +
  #geom_line(colour = "#e7298a") +
  geom_line() +
  #scale_colour_brewer(palette = "Dark2") +
  scale_color_manual(values = c("Pyr-ATN" = "#e7298a", "ATN" = "#66a61e")) +
  #scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, direction = -1) +
  #scale_color_viridis_d(option = "magma", end = 0.8, direction = 1) +
  scale_alpha_discrete(range = c(0.4, 1)) +
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
    color = "Net type",
    linetype = "ATN action",
    alpha = "Inhibition\nduration"
  )




df_summary <- df_plot |>
  mutate(
    itn_type = factor(
      itn_type, levels = c("Pyr-only", "Pyr-PBO", "Pyr-CFP", "Pyr-ATN", "ATN")
    ),
    inhib_action = factor(
      inhib_action, levels = c("Concurrent", "Short", "Long")
    ),
    atn_action = factor(
      atn_action, levels = c("Both", "Inhibitory", "EIP")
    ),
    year = (time - 1095) / 365
  ) |>
  filter(year <= 3) |>
  filter(year >= 0) |>
  group_by(itn_type, resistance, EIR, atn_action, inhib_action) |>
  summarise(
    total_base = sum(base_clin_inc_0_Inf, na.rm = TRUE),
    total_averted = sum(avert_clin_inc_0_Inf, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(avert_prop = total_averted / total_base) |>
  mutate(scenario = paste(itn_type, atn_action, inhib_action)) |>
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "Pyr-only Both Long"            , "Pyr-PBO Both Long"       ,
        "Pyr-CFP Both Long"             , "Pyr-ATN Both Long"       ,
        "Pyr-ATN Both Short"            , "Pyr-ATN Both Concurrent" ,
        "Pyr-ATN Inhibitory Long"       , "Pyr-ATN Inhibitory Short",
        "Pyr-ATN Inhibitory Concurrent" , "Pyr-ATN EIP Long"        ,
        "ATN Both Long"             ,
        "ATN Both Short"            , "ATN Both Concurrent" ,
        "ATN Inhibitory Long"       , "ATN Inhibitory Short",
        "ATN Inhibitory Concurrent" , "ATN EIP Long"
      )
    )
  )

df_summary |>
  filter(atn_action == "Both") |>
  ggplot(aes(x = scenario,
             y = total_averted/100,
             fill = itn_type)) +
  geom_col(aes(alpha = inhib_action)) +
  geom_text(
    aes(label = scales::percent(avert_prop, accuracy = 0.1)),
    vjust = -0.5, size = 2
  ) +
  scale_fill_brewer(palette = "Dark2") +
  scale_alpha_discrete(range = c(0.4, 1)) +
  facet_grid(EIR ~ resistance,
             labeller = labeller(
               EIR = eir_labels,
               resistance = res_pct_labeller
             )
  ) + #, scales = "free_y") +
  theme_minimal() +
  labs(
    y = "Clinical cases averted per 1,000 over 3 years",
    fill = "Net type",
    alpha = "Inhibition\nduration"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  ylim(0, 1250)






library(ggpattern)

df_summary |>
  #filter(atn_action == "Both") |>
  ggplot(aes(x = scenario,
             y = total_averted/100,
             fill = itn_type)) +
  geom_col_pattern(aes(alpha = inhib_action, pattern = atn_action),
                   pattern_fill = "black", pattern_colour = NA) +
  geom_text(
    aes(label = scales::percent(avert_prop, accuracy = 1)),#accuracy = 0.1)),
    vjust = -0.25, hjust = -0.25, size = 2, angle = 60
  ) +
  scale_pattern_manual(values=c(NA, 'circle', 'stripe')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_alpha_discrete(range = c(0.4, 1)) +
  facet_grid(EIR ~ resistance,
             labeller = labeller(
               EIR = eir_labels,
               resistance = res_pct_labeller
             )
  ) +
  theme_minimal() +
  labs(
    y = "Clinical cases averted per 1,000 over 3 years",
    fill = "Net type",
    alpha = "Inhibition\nduration",
    pattern = "ATN action"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  guides(
    alpha = guide_legend(
      override.aes = list(
        pattern = "none",
        pattern_fill = NA,
        pattern_colour = NA
      )
    ),
    fill = guide_legend(
      override.aes = list(
        pattern = "none",
        pattern_fill = NA,
        pattern_colour = NA
      )
    ),
    pattern = guide_legend(
      override.aes = list(
        alpha = 0.7
      )
    )
  ) +
  ylim(0, 1250)



















# df_plot |>
#   mutate(
#     itn_type = factor(
#       itn_type,
#       levels = c("Pyr-only", "Pyr-PBO", "Pyr-CFP", "Pyr-ATNv0", "Pyr-ATNv1", "ATN")
#     ),
#     year = (time - 1095) / 365
#   ) |>
#   filter(year >= -1) |>
#   ggplot(aes(x = year, y = avert_clin_inc_0_Inf, color = itn_type)) +
#   geom_line(linewidth = 0.8) +
#   scale_color_viridis_d(option = "turbo", begin = 0, end = 0.5, direction = -1) +
#   #scale_color_viridis_d(option = "plasma", end = 0.95, direction = -1) +
#   facet_grid(EIR ~ resistance,
#              labeller = labeller(
#                EIR = eir_labels,
#                resistance = res_pct_labeller
#              )
#   ) + #, scales = "free_y") +
#   theme_minimal(base_size = 13) +
#   labs(
#     x = "Years since distribution",
#     y = "Clinical cases averted",
#     color = "Net type"
#     #title = "Pyrethroid resistance (columns) vs baseline EIR (rows)"
#   )

df_summary <- df_plot |>
  mutate(
    itn_type = factor(
      itn_type,
      levels = c("Pyr-only", "Pyr-PBO", "Pyr-CFP",
                 "Pyr-ATN", "Pyr-ATN(EIP)", "Pyr-ATN(inhib)",
                 "ATN", "ATN(EIP)", "ATN(inhib)")
    ),
    year = (time - 1095) / 365
  ) |>
  filter(year <= 3) |>
  filter(year >= 0) |>
  filter(itn_type != "Pyr-ATN(EIP)") |>
  filter(itn_type != "Pyr-ATN(inhib)") |>
  filter(itn_type != "ATN(EIP)") |>
  filter(itn_type != "ATN(inhib)") |>
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
    vjust = -0.5, size = 2
  ) +
  scale_fill_viridis_d(option = "turbo", begin = 0.1, end = 0.9, direction = -1) +
  #scale_fill_viridis_d(option = "turbo", begin = 0, end = 0.5, direction = -1) +
  facet_grid(EIR ~ resistance,
             labeller = labeller(
               EIR = eir_labels,
               resistance = res_pct_labeller
             )
  ) + #, scales = "free_y") +
  theme_minimal() +
  labs(
    y = "Clinical cases averted per 1,000 over 3 years",
    fill = "Net type"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  ylim(0, 1250)



df_plot |>
  mutate(
    itn_type = factor(
      itn_type,
      levels = c("Pyr-only", "Pyr-PBO", "Pyr-CFP",
                 "Pyr-ATN", "Pyr-ATN(inhib)", "Pyr-ATN(EIP)",
                 "ATN", "ATN(inhib)", "ATN(EIP)")
    ),
    year = (time - 1095) / 365
  ) |>
  filter(year >= -1) |>
  filter(itn_type != "Pyr-only") |>
  filter(itn_type != "Pyr-PBO") |>
  filter(itn_type != "Pyr-CFP") |>
  filter(itn_type != "ATN") |>
  filter(itn_type != "ATN(EIP)") |>
  filter(itn_type != "ATN(inhib)") |>
  ggplot(aes(x = year, y = pfpr2to10, color = itn_type)) +
  geom_line() +
  scale_color_viridis_d(option = "mako", begin = 0.1, end = 0.8, direction = -1) +
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


df_summary <- df_plot |>
  mutate(
    itn_type = factor(
      itn_type,
      levels = c("Pyr-only", "Pyr-PBO", "Pyr-CFP",
                 "Pyr-ATN", "Pyr-ATN(inhib)", "Pyr-ATN(EIP)",
                 "ATN", "ATN(inhib)", "ATN(EIP)")
    ),
    year = (time - 1095) / 365
  ) |>
  filter(year <= 3) |>
  filter(year >= 0) |>
  filter(itn_type != "Pyr-only") |>
  filter(itn_type != "Pyr-PBO") |>
  filter(itn_type != "Pyr-CFP") |>
  filter(itn_type != "ATN") |>
  filter(itn_type != "ATN(EIP)") |>
  filter(itn_type != "ATN(inhib)") |>
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
    vjust = -0.5, size = 2
  ) +
  scale_fill_viridis_d(option = "mako", begin = 0.1, end = 0.8, direction = -1) +
  #scale_fill_viridis_d(option = "turbo", begin = 0, end = 0.5, direction = -1) +
  facet_grid(EIR ~ resistance,
             labeller = labeller(
               EIR = eir_labels,
               resistance = res_pct_labeller
             )
  ) + #, scales = "free_y") +
  theme_minimal() +
  labs(
    y = "Clinical cases averted per 1,000 over 3 years",
    fill = "Net type"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  ylim(0, 1250)
