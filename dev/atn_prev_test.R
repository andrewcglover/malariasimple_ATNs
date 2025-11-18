# atn_prev_test.R
library(dplyr)
library(ggplot2)

read_itn_pars <- FALSE

n_steps <- 365 * 9

if (read_itn_pars) {
  only_pars <- readRDS("./dev/net_params/pyrethroid_uncertainty2.rds")
  pbo_pars <- readRDS("./dev/net_params/pbo_uncertainty_using_pyrethroid_dn0_for_mn_durability2.rds")
  cfp_pars <- readRDS("./dev/net_params/pyrrole_uncertainty_using_pyrethroid_dn0_for_mn_durability2.rds")
}

params_atn <- get_parameters(n_days = n_steps,
                             spor_len = 10,
                             gamma_atn = (log(2) / (2.64 * 365))) |>
  set_bednets(days = 1095,
              coverages = 0.9,
              gamman =   2.64 * 365,
              retention =   588,
              dn0 = 0,
              rn = 0.56,
              rnm = 0.24) |>
  set_equilibrium(init_EIR = 10)
params_atn$max_atn_cov <- 0
out_atn <- run_atn_simulation(params_atn)
df_atn <- as.data.frame(out_atn)

# params_base <- get_parameters(n_days = n_steps) |>
#   set_bednets(days = 1095,
#               coverages = 0.9,
#               gamman =   2.64 * 365,
#               retention =   588,
#               dn0 = 0.1,
#               rn = 0.56,
#               rnm = 0.24) |>
#   set_equilibrium(init_EIR = 10)
# out_base <- run_simulation(params_base)
# df_base <- as.data.frame(out_base)

params_base <- get_parameters(n_days = n_steps, Q0_atn = 0) |>
  set_bednets(days = 1095,
              coverages = 0.9,
              gamman =   2.64 * 365,
              retention =   588,
              dn0 = 0.41,
              rn = 0.56,
              rnm = 0.24) |>
  set_equilibrium(init_EIR = 10)
out_base <- run_atn_simulation(params_base)
df_base <- as.data.frame(out_base)

df_combined <- bind_rows(
  df_base |> mutate(model = "Baseline"),
  df_atn  |> mutate(model = "ATN")
)

df_combined$pfpr2to10 <- df_combined$n_detect_730_3650 / df_combined$n_730_3650

ggplot(df_combined, aes(x = time, y = pfpr2to10, color = model)) +
  geom_line(size = 1) +
  theme_minimal(base_size = 13) +
  labs(
    x = "Time (days)",
    y = "PfPR2–10",
    color = "Model"
    )

