# -------------------------
# Helper function: poverty headcounts
# -------------------------
compute_headcounts <- function(PL_vals, Average, a, b, sigma, alpha, beta) {
  data.frame(
    Observed = PL_vals$obs,
    FISK_H   = sapply(PL_vals$pl,  function(p) CDF_Fisk(y=p, alpha = a, beta = b) * 100),
    LN_H     = sapply(PL_vals$pl,  function(p) pnorm((log(p / Average)/sigma) + sigma/2) * 100),
    NP_H     = sapply(PL_vals$pl,  function(p) CDF_NP(alpha, beta, q=p) * 100),
  )
}

# -------------------------
# Helper function: parameters summary
# -------------------------
compute_param_summary <- function(a_FISK, b_FISK, mu_LN, sigma_LN, alpha_NP, beta_NP,
                                  se_a = NA, se_b = NA, se_mu = NA, se_sigma = NA, se_alpha = NA, se_beta = NA) {
  data.frame(
    row = c("Shape/μ", "Scale/σ"),
    Observed_FISK = c(a_FISK, b_FISK),
    SE_FISK = c(se_a, se_b),
    Observed_LN = c(mu_LN, sigma_LN),
    SE_LN = c(se_mu, se_sigma),
    Observed_NP = c(alpha_NP, beta_NP),
    SE_NP = c(se_alpha, se_beta)
  )
}

# -------------------------
# Helper function: Headcount Ratio Variance
# -------------------------
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

# -------------------------
# Helper function: Confidence interval
# -------------------------
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

