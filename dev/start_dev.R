#start_dev.R

library(devtools)
library(roxygen2)
library(odin2)
library(dust2)

devtools::load_all()

source("dev/reload_atn_model.R")

reload_atn_model()
