#' Confidence Intervals for Poverty Headcount Estimates
#'
#' A dataset providing 95\% confidence intervals for observed and model-based
#' headcount poverty ratios. Each row corresponds to a country-year and a
#' poverty threshold.
#'
#' @format A data frame with one row per threshold per country-year:
#'
#' \describe{
#'   \item{Country}{Character. Country code and survey year (e.g., `"at94"`).}
#'
#'   \item{threshold}{Numeric. Poverty threshold expressed as a fraction
#'         of mean income (e.g., 0.2 for 20\%).}
#'
#'   \item{Obs_lower}{Numeric. Lower 95\% CI for observed headcount.}
#'   \item{Obs_upper}{Numeric. Upper 95\% CI for observed headcount.}
#'
#'   \item{FISK_H_lower}{Numeric. Lower 95\% CI for Fisk-model headcount.}
#'   \item{FISK_H_upper}{Numeric. Upper 95\% CI for Fisk-model headcount.}
#'
#'   \item{LN_H_lower}{Numeric. Lower 95\% CI for Lognormal-model headcount.}
#'   \item{LN_H_upper}{Numeric. Upper 95\% CI for Lognormal-model headcount.}
#'
#'   \item{NP_H_lower}{Numeric. Lower 95\% CI for NP-model headcount.}
#'   \item{NP_H_upper}{Numeric. Upper 95\% CI for NP-model headcount.}
#' }
#'
#' @details
#' Confidence intervals are computed using appropriate standard errors
#' for each distribution and are useful for inference and model comparison.
#'
#' Poverty thresholds typically include 20\%, 50\%, and 80\% of the mean,
#' but the dataset can accommodate arbitrary values.
#'
#' @usage data("CI_limited_data")
#'
#' @examples
#' data(CI_limited_data)
#' head(CI_limited_data)
#'
"CI_limited_data"
