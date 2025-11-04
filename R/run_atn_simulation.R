#' @title Run deterministic malariasimple
#' @description Generates and runs the deterministic malariasimple model
#' @param params List of model parameters generated using get_parameters() and other helper functions
#' @param n_particles Number of simulations to perform
#' @param full_output Boolean variable stating whether model should output model outputs, or a selected few.
#' @examples
#' params <- get_parameters() |>
#'           set_equilibrium(init_EIR = 10)
#' simulation_output <- run_simulation(params)
#' @export

run_atn_simulation <- function(params, n_particles = 1, full_output = FALSE){
  if(!params$stochastic){
    gen <- malariasimple_atn_deterministic
  } else if(params$stochastic){
    gen <- malariasimple_atn_stochastic
  }
  sys <- dust2::dust_system_create(gen(), params, n_particles = n_particles, dt = 1/params$tsd)

  dust2::dust_system_set_state_initial(sys)
  time <- 0:params$n_days
  out <- dust2::dust_system_simulate(sys, time)

  if(n_particles > 1){
    out <- aperm(out, perm = c(3, 1, 2))
    out <- out[-1,,]
  } else{
    out <- aperm(out,c(2,1))
    out <- out[-1,]
  }
  #Add time to output
  time <- 1:params$n_days
  time_slice <- array(time, dim = c(1, length(time), n_particles)) |> drop()
  out <- abind::abind(out, time_slice, along = 2)  # Add "time" as the last column of the second dimension

  #Add colnames to output
  colnames <- get_output_colnames(sys, params)
  dimnames(out)[[2]] <- c(colnames,"time")

  #Only output select variables (unless requested otherwise)
  if(!full_output){
    # selected_cols <- c("time","EIR_mean","natural_deaths","mu_mosq",
    #                    grep("_count$", colnames(out), value = TRUE),
    #                    "ica_mean","icm_mean","ib_mean", "id_mean",
    #                    grep("^n_", colnames(out), value = TRUE),
    #                    "EL","LL","PL","Sv","Pv","Iv","mv")
    selected_cols <- c(
      "time", "EIR_mean", "natural_deaths", "mu_mosq",
      grep("_count$", colnames(out), value = TRUE),
      "ica_mean", "icm_mean", "ib_mean", "id_mean",
      grep("^n_", colnames(out), value = TRUE),
      "EL", "LL", "PL", "Svtot", "Evtot", "Ivtot", "mv", "Q_atn"
    )
    if(n_particles == 1){
      out <- out[,selected_cols]
    } else {
      out <- out[,selected_cols,]
    }

  }
  return(out)
}

#' @title Get column names for simulation output
#' @param sys Dust system. Output from dust2::dust_system_create()
#' @param params Simulation parameters
#' @export
get_output_colnames <- function(sys, params){
  ## Get column names
  colname_ls <- dust2::dust_unpack_index(sys)
  colname_length <- sum(sapply(colname_ls, length))
  index <- vector("character", colname_length)
  pos <- 1

  for (i in seq_along(colname_ls)) {
    len <- length(colname_ls[[i]])
    index[pos:(pos + len - 1)] <- names(colname_ls)[i]
    pos <- pos + len
  }

  #Improve naming of user-defined outputs
  index[index == "n_ud_prev"] <- paste("n",params$age_vector[params$min_age_prev], params$prevalence_rendering_max_age,sep="_")
  index[index == "n_ud_detect_prev"] <- paste("n_detect",params$age_vector[params$min_age_prev], params$prevalence_rendering_max_ages,sep="_")
  index[index == "n_ud_inc"] <- paste("n_clin_inc",params$age_vector[params$min_age_inc], params$clin_inc_rendering_max_ages,sep="_")
  return(index)
}

