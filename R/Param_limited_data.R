#' Parametric Income Distribution Estimates
#'
#' A dataset containing estimated parameters and standard errors for three
#' fitted income distributions (Fisk/log-logistic, Lognormal, and NP-Pareto)
#' for each country-year.
#'
#' @format A data frame with the following columns:
#'
#' \describe{
#'
#'   \item{Parameter}{Character. The name of the parameter (e.g.,
#'       `"shape/μ"`, `"scale"`, `"σ"`, depending on distribution).}
#'
#'   \item{Fisk}{Numeric. Estimated parameter value for the Fisk distribution.}
#'   \item{Fisk_SE}{Numeric. Standard error of the Fisk parameter estimate.}
#'
#'   \item{LN}{Numeric. Estimated parameter value for the Lognormal distribution.}
#'   \item{LN_SE}{Numeric. Standard error for the Lognormal parameter estimate.}
#'
#'   \item{NP}{Numeric. Estimated parameter value for the NP-Pareto model.}
#'   \item{NP_SE}{Numeric. Standard error of the NP-Pareto parameter estimate.}
#'
#'   \item{Country}{Character. Country identifier (e.g., `"at94"`).}
#' }
#'
#' @details
#' These parameters are obtained by maximum likelihood or equivalent
#' estimation methods and are used for computing fitted headcounts,
#' inequality measures, and income distribution summaries.
#'
#' @usage data("Param_limited_data")
#'
#' @examples
#' data(Param_limited_data)
#' head(Param_limited_data)
#'
"Param_limited_data"
