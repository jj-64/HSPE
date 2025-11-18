#' Fast Estimation of Parameters Using Mean and Gini
#'
#' Computes:
#' distribution parameters as in \code{\link{fit_Fisk_param}},
#' or \code{\link{fit_LN_param}} or \code{\link{fit_NP_param}}.
#'
#' Optionally computes delta-method standard errors.
#'
#' @param dist character vector, one or multiple of "LN", "NP", "FISK"
#' @param mean_y Mean of the data (must be > 0)
#' @param Gini Gini coefficient in (0,1)
#' @param se_mean Standard error of mean (optional)
#' @param se_Gini Standard error of Gini (optional)
#'
#' @return A list with:
#'   \describe{
#'     \item{par}{c(shape, scale)}
#'     \item{se}{c(se_shape, se_scale)}
#'   }
#' @export
#' @examples
#' res = estimate_param(dist = "NP", mean_y = 6000, Gini = 0.36, se_mean = 250, se_Gini = 0.01)
#' # res$NP
#' # $par
#' # shape       scale
#' # 2.128961 2523.457703
#' #
#' # $se
#' # se_shape se_scale
#' # 0.04810901 125.06197149
#'
#'
#' res = estimate_param(dist = c("LN","NP"), mean_y = 6000, Gini = 0.36, se_mean = 250, se_Gini = 0.01)
#' # $FISK
#' # [1] NA
#' #
#' # $LN
#' # $LN$par
#' # mu    sigma
#' # 8.480773 0.661426
#' #
#' # $LN$se
#' # se_mu       se_sigma
#' # 0.04367100 0.01977307
#' #
#' #
#' # $NP
#' # $NP$par
#' # shape       scale
#' # 2.128961 2523.457703
#' #
#' # $NP$se
#' # se_shape se_scale
#' # 0.04810901 125.06197149
estimate_param <- function(dist = c("LN", "NP", "FISK"), mean_y, Gini, se_mean = NA, se_Gini = NA){
  #dist <- match.arg(dist)
  #if(! (all(dist) %in% c("LN", "NP", "FISK"))) stop("dist should be one or multiple of LN, FISK, NP")
  results = list("FISK" = NA, "LN" =NA, "NP" = NA)

  if("FISK" %in% dist){
    results[["FISK"]] = fit_Fisk_param (mean_y, Gini, se_mean = se_mean, se_Gini = se_Gini)
  }

  if("LN" %in% dist){
    results[["LN"]] = fit_LN_param(mean_y, Gini, se_mean = se_mean, se_Gini = se_Gini)
  }

  if( "NP" %in% dist){
    results[["NP"]]  = fit_NP_param(mean_y, Gini, se_mean = se_mean, se_Gini = se_Gini)
  }

  return(results)
}
