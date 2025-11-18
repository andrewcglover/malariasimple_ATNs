library(dplyr)
library(tidyr)
library(ggplot2)
library(future)
library(future.apply)

# ------------------------------------------------
# USER SETTINGS (edit easily)
# ------------------------------------------------
read_itn_pars <- TRUE
n_steps <- 365 * 9

eir_levels <- c(5, 40, 100)   # YOU CAN CHANGE THIS
res_levels <- c(0, 0.5, 1)    # YOU CAN CHANGE THIS

workers <- 6                  # number of CPU cores
plan(multisession, workers = workers)
# ------------------------------------------------


# ------------------------------------------------
# LOAD PARAMETER TABLES
# ------------------------------------------------
if (read_itn_pars) {
  only_pars <- readRDS("./dev/net_params/dat_res_only.rds")
  pbo_pars  <- readRDS("./dev/net_params/dat_res_pbo.rds")
  cfp_pars  <- readRDS("./dev/net_params/dat_res_cfp.rds")
}


# ------------------------------------------------
# PRECOMPUTE MEDIAN PARAMETERS PER RESISTANCE
# ------------------------------------------------
only_med <- only_pars |>
  group_by(resistance) |>
  summarise(
    dn0 = median(dn0),
    rn0 = median(rn0),
    gamman = median(gamman),
    .groups = "drop"
  )

pbo_med <- pbo_pars |>
  group_by(resistance) |>
  summarise(
    dn0 = median(dn0),
    rn0 = median(rn0),
    gamman = median(gamman),
    .groups = "drop"
  )

cfp_med <- cfp_pars |>
  group_by(resistance) |>
  summarise(
    dn0 = median(dn0),
    rn0 = median(rn0),
    gamman = median(gamman),
    .groups = "drop"
  )


# ------------------------------------------------
# RUN ATN ONCE PER EIR, THEN DUPLICATE FOR RESISTANCE
# ------------------------------------------------
atn_runs <- lapply(eir_levels, function(eir) {

  params_atn <- get_parameters(
    n_days = n_steps,
    spor_len = 10,
    gamma_atn = (log(2) / (2.64 * 365))
  ) |>
    set_bednets(
      days = 1095, coverages = 0.9,
      gamman = 2.64 * 365,
      retention = 588,
      dn0 = 0,
      rn  = 0.56,
      rnm = 0.24
    ) |>
    set_equilibrium(init_EIR = eir)

  params_atn$max_atn_cov <- 0

  base_atn <- run_atn_simulation(params_atn) |>
    as.data.frame() |>
    mutate(
      model = "ATN",
      itn_type = "ATN",
      pfpr2to10 = n_detect_730_3650 / n_730_3650,
      EIR = eir
    )

  # duplicate across all resistance levels
  bind_rows(lapply(res_levels, function(r) {
    base_atn |> mutate(resistance = r)
  }))
})

df_atn_all <- bind_rows(atn_runs)


# ------------------------------------------------
# CREATE TASK GRID FOR PARALLEL ITN RUNS
# ------------------------------------------------
task_grid <- expand.grid(
  EIR        = eir_levels,
  resistance = res_levels,
  itn_type   = c("Pyr-only", "Pyr-PBO", "Pyr-CFP"),
  stringsAsFactors = FALSE
)


# ------------------------------------------------
# PARALLEL ITN RUNS
# ------------------------------------------------
itn_runs <- future_lapply(seq_len(nrow(task_grid)), function(i) {

  task <- task_grid[i, ]

  # Select correct median parameter row
  pars <- switch(task$itn_type,
                 "Pyr-only" = only_med |> filter(resistance == task$resistance),
                 "Pyr-PBO"  = pbo_med  |> filter(resistance == task$resistance),
                 "Pyr-CFP"  = cfp_med  |> filter(resistance == task$resistance)
  )

  # Build parameter object
  params_itn <- get_parameters(n_days = n_steps, Q0_atn = 0) |>
    set_bednets(
      days      = 1095,
      coverages = 0.9,
      gamman    = pars$gamman * 365,   # converting from years if needed
      retention = 588,
      dn0       = pars$dn0,
      rn        = pars$rn0,
      rnm       = 0.24
    ) |>
    set_equilibrium(init_EIR = task$EIR)

  # Run model
  out <- run_atn_simulation(params_itn) |> as.data.frame()

  # Annotate
  out |> mutate(
    pfpr2to10 = n_detect_730_3650 / n_730_3650,
    model     = "ITN",
    itn_type  = task$itn_type,
    resistance = task$resistance,
    EIR        = task$EIR
  )
})

df_itn_all <- bind_rows(itn_runs)


# ------------------------------------------------
# COMBINE ATN + ITN
# ------------------------------------------------
df_plot <- bind_rows(df_atn_all, df_itn_all)


# ------------------------------------------------
# FACET PLOT: resistance × EIR
# ------------------------------------------------
ggplot(df_plot, aes(x = time, y = pfpr2to10, color = itn_type)) +
  geom_line(linewidth = 0.7) +
  facet_grid(EIR ~ resistance) + #, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    x = "Time (days)",
    y = "PfPR2–10",
    color = "ITN type",
    title = "ATN vs ITN performance across resistance and EIR levels"
  )
