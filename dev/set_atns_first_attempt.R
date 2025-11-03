#' Set ATN (Antimalarial-Treated Net) parameters
#' @export
#' @param params malariasimple parameters
#' @param days Vector of days on which ATN distribution events occur
#' @param coverages Vector of coverage proportions for each ATN distribution (must be same length as `days`)
#' @param lambda_atn Mean rate of loss of ATNs (1 / mean duration of use, in days)
#' @param p_atn Probability that a mosquito attempting to bite under an ATN is exposed to antimalarial drugs
#' @param phi_atn Probability a mosquito attempts to bite a human in bed (usually equal to `phi_bednets`)
#'
#' @returns Updates the input parameter list to include ATN parameters
#'
#' @examples
#' params <- get_parameters(n_days = 500)
#' params <- set_atns(params,
#'                    days = c(100),
#'                    coverages = c(0.4),
#'                    lambda_atn = 1/(2*365),
#'                    p_atn = 0.8)
set_atns <- function(params,
                     days = NULL,
                     coverages = NULL,
                     lambda_atn = 1 / (3 * 365), # mean 3 years
                     p_atn = 0.8,
                     phi_atn = NULL) {

  #---------------------- Sanity checks ---------------------------------
  if (params$equilibrium_set == 1) warning("Equilibrium must be set last")

  n_dists <- length(days)
  if (length(coverages) != n_dists)
    stop("`days` and `coverages` must have the same length")

  if (any(days %% 1 != 0))
    stop("'days' must be integer values")

  if (any(diff(days) < 1))
    stop("'days' must be unique and in chronological order")

  if (max(coverages) > 1)
    stop("'coverages' cannot exceed 1")

  if (is.null(phi_atn)) {
    if (!is.null(params$phi_bednets)) {
      phi_atn <- params$phi_bednets
    } else {
      stop("`phi_atn` not provided and no phi_bednets found in params")
    }
  }

  #---------------------- Define helper function ------------------------
  get_atn_coverage <- function(n_days, days, coverages, lambda_atn) {
    cov_daily <- rep(0, n_days)
    for (i in seq_along(days)) {
      day0 <- days[i]
      cov_daily[day0:n_days] <- cov_daily[day0:n_days] +
        coverages[i] * exp(-lambda_atn * (0:(n_days - day0)))
    }
    pmin(cov_daily, 1) # cap at 1
  }

  #---------------------- Set parameters --------------------------------
  atn_cov <- get_atn_coverage(params$n_days, days, coverages, lambda_atn)

  params$atn_days <- days
  params$atn_coverages <- coverages
  params$lambda_atn <- lambda_atn
  params$p_atn <- p_atn
  params$phi_atn <- phi_atn
  params$Q0_atn <- max(coverages)
  params$atn_cov_daily <- c(0, atn_cov)

  # Update flags
  params$atn_set <- 1
  params$equilibrium_set <- 0

  message("ATN parameters set successfully for ", n_dists, " distribution(s).")
  return(params)
}
