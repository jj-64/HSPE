#' Household Summary Statistics from LIS Microdata
#'
#' A dataset containing summary income distribution statistics computed
#' from LIS (Luxembourg Income Study) household microdata for multiple
#' countries and survey years. Each row corresponds to one country-year
#' dataset and includes measures of central tendency, inequality, distributional
#' shares, and poverty indicators.
#'
#' @format A data frame with one row per country-year and the following columns:
#'
#' \describe{
#'
#'   \item{Country}{Character. Country code and survey year (e.g., `"at94"`).}
#'
#'   \item{N}{Numeric. Number of observations in the micro dataset.}
#'
#'   \item{Mean}{Numeric. Mean disposable household income.}
#'   \item{Mean_SE}{Numeric. Standard error of the mean.}
#'
#'   \item{median}{Numeric. Median income.}
#'
#'   \item{Gini}{Numeric. Gini coefficient of inequality.}
#'   \item{Gini_SE}{Numeric. Standard error of the Gini index.}
#'
#'   \item{d1}{Numeric. Share of total income held by the poorest 10\% of households.}
#'   \item{d2}{Numeric. Share of total income held by the second decile.}
#'   \item{d3}{Numeric. Share of total income held by the third decile.}
#'   \item{d4}{Numeric. Share of total income held by the fourth decile.}
#'   \item{d5}{Numeric. Share of total income held by the fifth decile.}
#'   \item{d6}{Numeric. Share of total income held by the sixth decile.}
#'   \item{d7}{Numeric. Share of total income held by the seventh decile.}
#'   \item{d8}{Numeric. Share of total income held by the eighth decile.}
#'   \item{d9}{Numeric. Share of total income held by the ninth decile.}
#'   \item{d10}{Numeric. Share of total income held by the richest 10\% of households (should sum to 100 across all deciles).}
#'
#'   \item{HC_20}{Numeric. Headcount poverty ratio at 20\% of mean income.}
#'   \item{HC_50}{Numeric. Headcount poverty ratio at 50\% of mean income.}
#'   \item{HC_80}{Numeric. Headcount poverty ratio at 80\% of mean income.}
#'
#'   \item{HC_SE_20}{Numeric. Standard error of the 20\% poverty headcount (binomial).}
#'   \item{HC_SE_50}{Numeric. Standard error of the 50\% poverty headcount (binomial).}
#'   \item{HC_SE_80}{Numeric. Standard error of the 80\% poverty headcount (binomial).}
#'
#'   \item{CI_lower_20}{Numeric. Lower 95\% confidence interval bound for HC\_20.}
#'   \item{CI_upper_20}{Numeric. Upper 95\% confidence interval bound for HC\_20.}
#'
#'   \item{CI_lower_50}{Numeric. Lower 95\% confidence interval bound for HC\_50.}
#'   \item{CI_upper_50}{Numeric. Upper 95\% confidence interval bound for HC\_50.}
#'
#'   \item{CI_lower_80}{Numeric. Lower 95\% confidence interval bound for HC\_80.}
#'   \item{CI_upper_80}{Numeric. Upper 95\% confidence interval bound for HC\_80.}
#'
#' }
#'
#' @details
#' The dataset summarizes income distribution characteristics including:
#' \itemize{
#'   \item central tendency (mean, median),
#'   \item inequality (Gini and standard error),
#'   \item income distribution shares by decile,
#'   \item poverty headcounts at thresholds of 20\%, 50\%, and 80\% of the mean,
#'   \item standard errors and confidence intervals of poverty measures.
#' }
#'
#' This dataset is intended for inequality analysis, cross-country comparisons,
#' and testing distributional methods.
#'
#' @usage data("SumData")
#'
#' @examples
#' data(SumData)
#' head(SumData)
#'
"SumData"
