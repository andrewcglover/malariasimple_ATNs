# --- Run both models for 1000 days ---

params_atn <- get_parameters(n_days = 1000, spor_len = 10) |>
  set_equilibrium(init_EIR = 10)
out_atn <- run_atn_simulation(params_atn)

params_det <- get_parameters(n_days = 1000) |>
  set_equilibrium(init_EIR = 10)
out_det <- run_simulation(params_det)

# --- Extract and align variables ---
t <- out_atn[, "time"]

Evtot <- out_atn[, "Sv"]
Pv <- out_det[, "Sv"]
Evtot <- out_atn[, "Iv"]
Pv <- out_det[, "Iv"]
Evtot <- out_atn[, "Evtot"]
Pv <- out_det[, "Pv"]

# --- Plot setup: leave room for legend on the right ---
#op <- par(no.readonly = TRUE)
#par(xpd = TRUE, mar = c(5, 5, 4, 8))  # more right margin

# --- Plot ---
matplot(
  t,
  cbind(Evtot, Pv),
  type = "l", lwd = 2, lty = 1,
  col = c("red", "black"),
  xlab = "Time (days)",
  ylab = "Latent mosquito stage",
  main = "Comparison of Ev vs Pv (ATN vs Legacy)"
)

# --- Legend placed outside to the right ---
legend("bottomright",
       #inset = c(-0.25, 0),
       legend = c("Evtot (ATN)", "Pv (Legacy)"),
       col = c("red", "black"),
       lty = 1, lwd = 2, bty = "n")

# --- Restore plotting parameters ---
#par(op)

# --- Numeric comparison ---
mean_rel_diff <- mean(abs(Evtot - Pv) / Pv)
cat("Mean relative difference:", mean_rel_diff, "\n")
