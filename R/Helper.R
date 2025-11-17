# Helper function: poverty headcounts -----------------------
#' Compute Model-Based Poverty Headcounts for Multiple Distributions
#'
#' This function computes poverty headcount ratios for one or more
#' parametric income distributions (Lognormal, Fisk, and New Pareto).
#' It evaluates the cumulative distribution function (CDF) at each
#' poverty line provided in `PL_vals` and returns model-specific
#' headcounts (in percent) alongside the observed headcount.
#'
#' Users may specify which models to compute using the `models` argument.
#' Required distribution parameters must be provided only for the
#' requested models.
#'
#' @param PL_vals A data frame containing poverty lines and observed
#'   headcounts. Must include:
#'   \describe{
#'     \item{pl}{Numeric vector of poverty lines.}
#'     \item{obs}{Observed headcounts associated with each poverty line.}
#'   }
#'
#' @param models Character vector specifying which models to compute.
#'   Options are `"LN"`, `"FISK"`, `"NP"`, or `"all"` (default).
#'   The matching is case-insensitive.
#'
#' @param Average Numeric. Mean income parameter required for the
#'   Lognormal (LN) model.
#'
#' @param Fisk_scale,Fisk_shape Numerics. Shape parameters for the
#'   Fisk distribution.
#'
#' @param LN_sigma Numeric. Standard deviation of log-income for the
#'   Lognormal model.
#'
#' @param NP_alpha,NP_beta Numerics. Parameters for the New Pareto (NP)
#'   distribution.
#'
#' @details
#' For each poverty line `p`, the function computes:
#' \itemize{
#'   \item \strong{FISK:} \eqn{F(p | \alpha, \beta)}
#'   \item \strong{LN:} \eqn{ \Phi( ( \log(p/Average) / \sigma ) + \sigma/2 ) }
#'   \item \strong{NP:} \eqn{F_{\text{NP}}(p | \alpha, \beta)}
#' }
#'
#' Required parameters depend on the selected models.
#' Missing parameters for requested models result in an error.
#'
#' @return A data frame with:
#'   \describe{
#'     \item{Observed}{Observed headcounts from `PL_vals`.}
#'     \item{FISK_H}{Headcount from the Fisk model (if requested).}
#'     \item{LN_H}{Headcount from the Lognormal model (if requested).}
#'     \item{NP_H}{Headcount from the New Pareto model (if requested).}
#'   }
#'   Row names correspond to the poverty lines in `PL_vals$pl`.
#'
#' @examples
#' # Example poverty lines and observed values
#' PL_vals <- data.frame(
#'   pl = c(2, 3, 5),
#'   obs = c(10, 20, 35)
#' )
#'
#' # Lognormal model only
#' compute_headcounts(
#'   PL_vals,
#'   models = "LN",
#'   Average = 50,
#'   LN_sigma = 0.6
#' )
#'
#' # Fisk and NP models together
#' compute_headcounts(
#'   PL_vals,
#'   models = c("FISK", "NP"),
#'   Fisk_scale = 2, Fisk_shape = 1.5,
#'   NP_alpha = 1.2, NP_beta = 0.8
#' )
#'
#' # All three models
#' compute_headcounts(
#'   PL_vals,
#'   models = "all",
#'   Average = 50, LN_sigma = 0.6,
#'   Fisk_scale = 2, Fisk_shape = 1.5,
#'   NP_alpha = 1.2, NP_beta = 0.8
#' )
#' compute_headcounts(PL_vals, models = "LN", Average = avg, LN_sigma = sigma)
#' compute_headcounts(PL_vals, models = c("FISK", "NP"), Fisk_scale = a, Fisk_shape = b, NP_alpha = na, NP_beta = nb)
#' compute_headcounts(PL_vals, models = "all", ...)
#' @export
compute_headcounts <- function(
    PL_vals,
    models = c("all", "LN", "FISK", "NP"),
    Average = NULL,
    Fisk_scale = NULL,
    Fisk_shape  = NULL,
    LN_sigma   = NULL,
    NP_alpha   = NULL,
    NP_beta    = NULL
) {

  # --- Normalize and validate models input ---
  models <- tolower(models)
  valid_models <- c("all", "ln", "fisk", "np", "lognormal", "Fisk", "NewPareto", "newpareto")
  if (!all(models %in% valid_models)) {
    stop("models must be one of: all, LN, FISK, NP")
  }
  if ("all" %in% models) {
    models <- c("ln", "fisk", "np")
  }

  # Prepare result container:
  HC <- data.frame(
    Observed = PL_vals$obs,
    row.names = PL_vals$pl
  )

  # --------------------------------------------------------
  # FISK headcounts
  # --------------------------------------------------------
  if ("fisk" %in% models) {
    if (is.null(Fisk_scale) || is.null(Fisk_shape)) {
      stop("FISK parameters (Fisk_scale, Fisk_shape) must be provided when models include 'FISK'.")
    }

    HC$FISK_H <- sapply(
      PL_vals$pl,
      function(p) CDF_Fisk(y = p, scale = Fisk_scale, shape = Fisk_shape) * 100
    )
  }

  # --------------------------------------------------------
  # Lognormal headcounts
  # --------------------------------------------------------
  if ("ln" %in% models) {
    if (is.null(Average) || is.null(LN_sigma)) {
      stop("LN parameters (Average, LN_sigma) must be provided when models include 'LN'.")
    }

    HC$LN_H <- sapply(
      PL_vals$pl,
      function(p) pnorm( (log(p / Average) / LN_sigma) + LN_sigma / 2 ) * 100
    )
  }

  # --------------------------------------------------------
  # New Pareto headcounts
  # --------------------------------------------------------
  if ("np" %in% models) {
    if (is.null(NP_alpha) || is.null(NP_beta)) {
      stop("NP parameters (NP_alpha, NP_beta) must be provided when models include 'NP'.")
    }

    HC$NP_H <- sapply(
      PL_vals$pl,
      function(p) CDF_NP(alpha = NP_alpha, beta = NP_beta, q = p) * 100
    )
  }

  rownames(HC) = PL_vals$ratio
  return(HC)
}

# Helper function: parameters summary -------------------------
#' Summarize Estimated Distribution Parameters and Standard Errors
#'
#' This function creates a compact summary table containing the estimated
#' parameters and their standard errors for three income distribution models:
#' Fisk, Lognormal, and New Pareto.
#'
#' The function outputs a two-row table:
#' \itemize{
#'   \item Row 1 contains the first parameter for each model
#'         (Fisk \eqn{\alpha}, Lognormal \eqn{\mu}, NP \eqn{\alpha})
#'   \item Row 2 contains the second parameter for each model
#'         (Fisk \eqn{\beta}, Lognormal \eqn{\sigma}, NP \eqn{\beta})
#' }
#'
#' @param Fisk_scale,Fisk_shape Numeric. Estimated parameters of the Fisk distribution.
#' @param LN_mu,LN_sigma Numeric. Estimated parameters of the Lognormal distribution.
#' @param NP_alpha,NP_beta Numeric. Estimated parameters of the New Pareto distribution.
#'
#' @param se_Fisk_alpha,se_Fisk_beta Numeric. Standard errors for the Fisk parameters.
#' @param se_LN_mu,se_LN_sigma Numeric. Standard errors for Lognormal parameters.
#' @param se_NP_alpha,se_NP_beta Numeric. Standard errors for New Pareto parameters.
#'   All SEs default to `NA` if not provided.
#'
#' @return A data frame with two rows and columns for:
#'   \describe{
#'     \item{Parameter}{Labels: `"alpha / μ"` and `"beta / σ"`}
#'     \item{Fisk}{Estimated Fisk parameters}
#'     \item{Fisk_SE}{Standard errors for Fisk parameters}
#'     \item{LN}{Estimated Lognormal parameters}
#'     \item{LN_SE}{Standard errors for Lognormal parameters}
#'     \item{NP}{Estimated New Pareto parameters}
#'     \item{NP_SE}{Standard errors for NP parameters}
#'   }
#'
#' @examples
#' compute_param_summary(
#'   Fisk_scale = 2.1, Fisk_shape = 1.4,
#'   LN_mu = 3.2, LN_sigma = 0.55,
#'   NP_alpha = 1.7, NP_beta = 0.9,
#'   se_Fisk_alpha = 0.05, se_Fisk_beta = 0.03
#' )
#'
#' @export
compute_param_summary <- function(
    Fisk_scale, Fisk_shape,
    LN_mu, LN_sigma,
    NP_alpha, NP_beta,
    se_Fisk_alpha = NA, se_Fisk_beta = NA,
    se_LN_mu = NA, se_LN_sigma = NA,
    se_NP_alpha = NA, se_NP_beta = NA
) {
  data.frame(
    Parameter = c("alpha / μ", "beta / σ"),
    Fisk     = c(Fisk_scale, Fisk_shape),
    Fisk_SE  = c(se_Fisk_alpha, se_Fisk_beta),
    LN       = c(LN_mu, LN_sigma),
    LN_SE    = c(se_LN_mu, se_LN_sigma),
    NP       = c(NP_alpha, NP_beta),
    NP_SE    = c(se_NP_alpha, se_NP_beta),
    row.names = NULL
  )
}

# Helper function: Headcount Ratio Variance-------------------------

headcounts_SE = function(H, N){
  if (H>1 | H <0) stop("H should be a proportion between 0 and 1")
  sqrt(H*(1-H)/N)
}

compute_headcounts_SE = function(H,N){
  data.frame(
FISK_SE = sapply(H, FUN = function(p) headcounts_SE(H=p, N=N)),
LN_SE   = NA,
NP_SE   = NA
)
}

library(dplyr)
##Function to expand data according to HH weight vector -----------
replicate_rows <- function(df) {
  df = df %>%
    slice(rep(1:n(), hwgt))  %>%
    ungroup()
  return(as.data.frame(df))
}

##Function to expand data according to HH weight vector -------------
replicate_rows_p <- function(df) {
  df = df %>%
    slice(rep(1:n(), pwgt))  %>%
    ungroup()
  return(as.data.frame(df))
}

# Helper function: Confidence interval-------------------------
#' Compute Confidence Bounds for Proportions or Estimates
#'
#' This function computes a two-sided (lower/upper) confidence interval
#' for an estimate `x` given its standard error `se`. It supports vectorized
#' inputs, meaning `x` and `se` can be numeric vectors of the same length.
#'
#' @param x Numeric vector. The point estimates.
#' @param se Numeric vector. The standard errors associated with each estimate.
#' @param significance Numeric value between 0 and 1 (default = 0.95).
#'                     The confidence level.
#'
#' @return A data.frame with two columns:
#'         - lower: lower confidence bound
#'         - upper: upper confidence bound
#'
#' @examples
#' conf_bound(0.5, 0.1)
#' conf_bound(c(0.5, 0.7), c(0.1, 0.05))
#'
conf_bound <- function(x, se, significance = 0.95) {

  # --- Input checks ---
  if (!is.numeric(x) || !is.numeric(se)) {
    stop("x and se must be numeric.")
  }
  if (length(x) != length(se)) {
    stop("x and se must have the same length.")
  }
  if (significance <= 0 || significance >= 1) {
    stop("significance must be between 0 and 1.")
  }

  # --- Compute z-value ---
  z <- abs(qnorm((1 - significance) / 2))

  # --- Confidence interval ---
  lower <- x - z * se
  upper <- x + z * se

  # Return a clean data frame
  return(data.frame(lower = lower, upper = upper))
}

