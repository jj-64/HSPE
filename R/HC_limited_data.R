#' Headcount Poverty Indicators from Limited Data
#'
#' A dataset containing observed and model-based headcount poverty indicators
#' for multiple countries and survey years. Headcount ratios are computed using
#' observed microdata, and fitted using three parametric income distributions:
#' Fisk (log-logistic), Lognormal, and Pareto type-NP.
#'
#' @format A data frame with one row per country-year and the following columns:
#'
#' \describe{
#'
#'   \item{Observed}{Numeric. Observed poverty headcount (share of population).}
#'   \item{FISK_H}{Numeric. Headcount predicted by the Fisk (log-logistic) distribution.}
#'   \item{LN_H}{Numeric. Headcount predicted by the Lognormal distribution.}
#'   \item{NP_H}{Numeric. Headcount predicted by the Pareto-NP model.}
#'
#'   \item{Country}{Character. Country code and survey year (e.g., `"at94"`).}
#'
#'   \item{FISK_H_SE}{Numeric. Standard error of the Fisk-model headcount.}
#'   \item{LN_H_SE}{Numeric. Standard error of the Lognormal-model headcount.}
#'   \item{NP_H_SE}{Numeric. Standard error of the NP-model headcount.}
#' }
#'
#' @details
#' This dataset contains point estimates and standard errors
#' of poverty headcounts computed from:
#' \itemize{
#'   \item observed microdata,
#'   \item fitted Fisk (log-logistic) distribution,
#'   \item fitted Lognormal distribution,
#'   \item fitted Pareto-NP distribution.
#' }
#'
#' These estimates are useful for comparing parametric and nonparametric
#' representations of income distributions.
#'
#' @usage data("HC_limited_data")
#'
#' @examples
#' data(HC_limited_data)
#' head(HC_limited_data)
#'
"HC_limited_data"
