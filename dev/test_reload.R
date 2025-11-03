# test script
params <- get_parameters() |> set_equilibrium(init_EIR = 10)
params$n_days <- 100
out <- run_atn_simulation(params)
plot(out[, "time"], out[, "Iv"], type = "l")
