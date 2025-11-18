# Helper: Fit a grouped model -----------------------------
#' Fit a Grouped-Income Distribution Model
#'
#' This function provides a unified wrapper to fit one of several
#' parametric income distributions (DA, SM, B2, GB2, FISK, LN, NP)
#' using grouped data (non-cumulative decile shares), mean income,
#' and Gini index.
#'
#' @param model Character scalar. One of:
#'   `"DA"`, `"SM"`, `"B2"`, `"GB2"`, `"FISK"`, `"LN"`, `"NP"`.
#' @param L_nonCum Numeric vector of non-cumulative income shares
#'   (e.g., decile shares).
#' @param mean_y Mean income.
#' @param Gini Observed Gini coefficient.
#' @param N Sample size (affects SEs of some fitters).
#'
#' @details
#' The function looks up the corresponding fitting routine from
#' `CDF_registry`, executes it, and returns parameter estimates
#' together with their standard errors.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{ok}{Logical indicating success.}
#'     \item{par}{Named vector of parameter estimates (if successful).}
#'     \item{se}{Standard errors corresponding to `par`.}
#'     \item{gini_est}{Model-implied Gini estimate.}
#'   }
#'
#' @examples
#' \dontrun{
#' fit_model_grouped("NP", L_nonCum = deciles, mean_y = 20000,
#'                   Gini = 0.32, N = 5000)
#' }
fit_model_grouped <- function(model = c("DA","SM","B2","GB2","FISK","LN","NP"),
                              L_nonCum,
                              mean_y,
                              Gini,
                              N)
{
  model <- match.arg(model)

  if (!model %in% names(CDF_registry)) {
    stop("Unknown model: ", model, ". Must be one of: ",
         paste(names(CDF_registry), collapse = ", "))
  }

  info <- CDF_registry[[model]]
  fitfun <- get(info$fitfun)

  tryCatch({
    fit <- fitfun(
      y       = L_nonCum,
      pc.inc  = mean_y,
      gini.e  = Gini,
      gini    = TRUE,
      N       = N,
      se.ewmd = TRUE,
      se.scale = TRUE,
      nrep    = 100
    )

    list(
      ok       = TRUE,
      par      = fit$ewmd.estimation["Coef.", info$params],
      se       = fit$ewmd.estimation["se",    info$params],
      gini_est = fit$gini.estimation["gini"]
    )

  }, error = function(e){
    message(model, " FAILED: ", e$message)
    list(ok = FALSE)
  })
}


