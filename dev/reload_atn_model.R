# Helper script to recompile and reload ATN models during development
# Safe to run manually in RStudio or source("dev/reload_atn_model.R")

reload_atn_model <- function() {
  message(" Rebuilding malariasimple_ATNs models...")

  # Step 1: regenerate Dust2 + C++ code
  odin2::odin_package(".")

  # Step 2: recompile & reload package
  devtools::clean_dll()
  devtools::load_all()

  message("Rebuild complete. Ready to run run_atn_simulation(params)")
}
