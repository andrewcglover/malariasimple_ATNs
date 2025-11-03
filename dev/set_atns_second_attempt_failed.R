#' #' Set ATN (Antimalarial-Treated Net) parameters
#' #'
#' #' @description
#' #' Defines coverage and entomological parameters for antimalarial-treated nets (ATNs),
#' #' analogous to how ITNs and SMC are handled in the model. Supports one or more discrete
#' #' distribution campaigns with exponential decay of coverage (defined by `lambda_atn`).
#' #'
#' #' @param params malariasimple parameters object
#' #' @param days Vector of days on which ATN distributions occur
#' #' @param coverages Vector of coverage proportions for each distribution
#' #' @param lambda_atn Mean rate of loss of ATNs (1 / mean duration of use, in days)
#' #' @param p_atn Probability that a mosquito feeding under an ATN is exposed to antimalarial drugs
#' #' @param phi_atn Probability a mosquito attempts to bite a human in bed (defaults to `phi_bednets`)
#' #'
#' #' @returns
#' #' The updated parameter list including ATN-specific quantities:
#' #' - `atn_days`, `atn_coverages`
#' #' - `lambda_atn`, `p_atn`, `phi_atn`
#' #' - `atn_cov_daily`, `atn_eff_cov_daily`
#' #' - `max_atn_cov`
#' #' - flags `atn_set`, `equilibrium_set`
#' #'
#' #' @export
#' #'
#' #' @examples
#' #' params <- get_parameters(n_days = 500)
#' #' params <- set_atns(params,
#' #'                    days = c(100, 300),
#' #'                    coverages = c(0.4, 0.3),
#' #'                    lambda_atn = 1 / (3 * 365),
#' #'                    p_atn = 0.8)
#' set_atns <- function(params,
#'                      days = NULL,
#'                      coverages = NULL,
#'                      lambda_atn = 1 / (3 * 365),  # mean 3 years
#'                      p_atn = 0.8,
#'                      phi_atn = NULL) {
#'
#'   ## ------------------ Sanity checks ------------------
#'   if (params$equilibrium_set == 1)
#'     warning("Equilibrium must be set last")
#'
#'   if (is.null(days) || is.null(coverages))
#'     stop("Both 'days' and 'coverages' must be specified for ATN distributions")
#'
#'   n_dists <- length(days)
#'   if (length(coverages) != n_dists)
#'     stop("'days' and 'coverages' must be the same length")
#'
#'   if (any(days %% 1 != 0))
#'     stop("'days' must be integer values")
#'
#'   if (any(diff(days) < 1))
#'     stop("'days' must be unique and in chronological order")
#'
#'   if (max(coverages) > 1)
#'     stop("'coverages' cannot exceed 1")
#'
#'   if (is.null(phi_atn)) {
#'     if (!is.null(params$phi_bednets)) {
#'       phi_atn <- params$phi_bednets
#'     } else {
#'       stop("`phi_atn` not provided and no `phi_bednets` found in params")
#'     }
#'   }
#'
#'   ## ------------------ Helper: exponential decay of coverage ------------------
#'   get_atn_coverage <- function(n_days, days, coverages, lambda_atn) {
#'     cov_daily <- rep(0, n_days)
#'     for (i in seq_along(days)) {
#'       day0 <- days[i]
#'       cov_daily[day0:n_days] <- cov_daily[day0:n_days] +
#'         coverages[i] * exp(-lambda_atn * (0:(n_days - day0)))
#'     }
#'     pmin(cov_daily, 1) # cap at 1
#'   }
#'
#'   ## ------------------ Compute daily coverage trajectory ------------------
#'   atn_cov <- get_atn_coverage(params$n_days, days, coverages, lambda_atn)
#'   max_atn_cov <- max(atn_cov)
#'
#'   ## ------------------ Store parameters in params list ------------------
#'   params$atn_days <- days
#'   params$atn_coverages <- coverages
#'   params$lambda_atn <- lambda_atn
#'   params$p_atn <- p_atn
#'   params$phi_atn <- phi_atn
#'   params$Q0_atn <- max(coverages)
#'
#'   params$atn_cov_daily <- c(0, atn_cov)
#'   params$atn_eff_cov_daily <- if (max_atn_cov > 0) c(0, atn_cov / max_atn_cov) else rep(0, params$n_days + 1)
#'   params$max_atn_cov <- max_atn_cov
#'
#'   ## ------------------ Integration with intervention system ------------------
#'   # num_int = 8 when ITNs + SMC + ATNs are all present
#'   params$num_int <- 8
#'
#'   ## ------------------ Flags ------------------
#'   params$atn_set <- 1
#'   params$equilibrium_set <- 0
#'
#'   message("ATN parameters set successfully for ", n_dists, " distribution(s).")
#'   return(params)
#' }

#' Set ATN (Antimalarial-Treated Net) parameters
#' @export
set_atns <- function(params,
                     days = NULL,
                     coverages = NULL,
                     lambda_atn = 1 / (3 * 365), # mean 3 years
                     p_atn = 0.8,
                     phi_atn = NULL) {

  #---------------------- Sanity checks ---------------------------------
  if (params$equilibrium_set == 1)
    warning("Equilibrium must be set last")

  n_dists <- length(days)
  if (length(coverages) != n_dists)
    stop("`days` and `coverages` must have the same length")
  if (any(days %% 1 != 0))
    stop("'days' must be integer values")
  if (any(diff(days) < 1))
    stop("'days' must be unique and in chronological order'")
  if (max(coverages) > 1)
    stop("'coverages' cannot exceed 1")

  # if phi_atn not given, fall back to phi_bednets
  if (is.null(phi_atn)) {
    if (!is.null(params$phi_bednets)) {
      phi_atn <- params$phi_bednets
    } else {
      stop("`phi_atn` not provided and no phi_bednets found in params")
    }
  }

  #---------------------- Define daily coverage -------------------------
  get_atn_coverage <- function(n_days, days, coverages, lambda_atn) {
    cov_daily <- rep(0, n_days)
    for (i in seq_along(days)) {
      day0 <- days[i]
      cov_daily[day0:n_days] <- cov_daily[day0:n_days] +
        coverages[i] * exp(-lambda_atn * (0:(n_days - day0)))
    }
    pmin(cov_daily, 1)
  }

  atn_cov <- get_atn_coverage(params$n_days, days, coverages, lambda_atn)

  #---------------------- Set parameters for odin -----------------------
  params$atn_days <- days
  params$atn_coverages <- coverages
  params$lambda_atn <- lambda_atn
  params$p_atn <- p_atn
  params$phi_atn <- phi_atn
  params$max_atn_cov <- max(atn_cov)
  params$atn_cov_daily <- c(0, atn_cov)  # length = n_days + 1 for odin2

  # flags
  params$atn_set <- 1
  params$equilibrium_set <- 0

  message("ATN parameters set successfully for ", n_dists, " distribution(s).")
  return(params)
}
