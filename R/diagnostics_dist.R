## Grouped data -----------------
#' Diagonistics for grouped data
#'
#' Compute multiple measures for diagnostics from grouped data
#'
#' @param L_obs observed Lorenz cumulative shares at p
#' @param Lorenz_fun Lorenz function(p, par): model Lorenz curve
#' @param par list of parameters estimated
#' @param N sample size (for grouped likelihood)
#' @export
#' @examples
#' observed = c(0.2,0.28,0.31,0.4,0.6)
#' Lorenz = function (p, pars) pnorm(qnorm(p) - pars$s) # cdf_registry[["LN"]]$lorenzfun
#' param_list = list(mu = 1, s =2)
#' diagnostics_grouped(L_obs = observed, Lorenz_fun = Lorenz, par = param_list, N=100)
#' # $KS
#' # [1] 0.3637596
#' #
#' # $MSE
#' # [1] 0.08773977
#' #
#' # $RMSE
#' # [1] 0.296209
#' #
#' # $logLik
#' # [1] -434.7308
#' #
#' # $AIC
#' # [1] 873.4617
#' #
#' # $BIC
#' # [1] 878.672
#' #
#' # $L_theo
#' # [1] 0.0005161882 0.0057947902 0.0227501319 0.0700257214 0.2362404159
#' #
#' # $L_obs
#' # [1] 0.20 0.28 0.31 0.40 0.60
diagnostics_grouped <- function(L_obs, Lorenz_fun, par, N) {

  # Model Lorenz curve
  p <- seq(0.1 , 1, by = (1/length(L_obs)) )

  L_obs = as.numeric(L_obs) ## remove any vector names
  if(L_obs[length(L_obs)]>=99) L_obs = L_obs/100  ## Make sure L_obs are not in %

  L_theo <- Lorenz_fun(p, par)

  # --- 1. KS-like statistic (Lorenz curve)
  KS <- max(abs(L_obs - L_theo))

  # --- 2. MSE on cumulative shares
  MSE <- mean((L_obs - L_theo)^2)
  RMSE <- sqrt(MSE)

  # --- 3. Grouped log-likelihood
  # Prob. mass in each bin = Δcdf
  # But we approximate Δcdf using ΔLorenz * (mean * population)
  # Simpler alternative: use ΔLorenz only (scale-free)
  dL <- diff(c(0, L_theo))
  dL[dL <= 0] <- 1e-12

  nj <- rep(N/length(p), length(p))  # equal decile sample count

  logLik <- sum(nj * log(dL))

  K <- length(par)
  AIC <- -2 * logLik + 2 * K
  BIC <- -2 * logLik + log(N) * K

  return(list(
    KS = KS,
    MSE = MSE,
    RMSE = RMSE,
    logLik = logLik,
    AIC = AIC,
    BIC = BIC,
    L_theo = L_theo,
    L_obs = L_obs
                )
    )
}

## Micro data----------------------
#' Diagnostics for micro data
#'
#' Compute multiple measures for diagnostics from micro data
#'
#' @param y observed individual income values
#' @param cdf function(y, param) returning vector of cdf values
#' @param pdf function(y, param) returning vector of PDF values
#' @param param list of named vector of parameters
#' @export
#' @examples
#' y = c(100, 200, 250, 300, 400)
#' param_list = list(a=1, b=2)
#' pdf_FISK_param = function(y, param){pdf_FISK(y=y, scale= param$a, shape = param$b)}
#' cdf_FISK_param = function(y, param){cdf_FISK(y=y, scale= param$a, shape = param$b)}
#' diagnostics_micro(y = y, cdf = cdf_FISK_param, pdf = pdf_FISK_param, param = param_list)
#' # $KS_stat
#' # [1] 0.7999
#' #
#' # $KS_pvalue
#' # [1] 1
#' #
#' # $logLik
#' # [1] -77.89517
#' #
#' # $AIC
#' # [1] 159.7903
#' #
#' # $BIC
#' # [1] 159.0092
diagnostics_micro <- function(y, cdf, pdf, param) {

  # -- Sort data for KS test
  y_sorted <- sort(y)
  Fn <- ecdf(y_sorted)(y_sorted)
  Ftheta <- cdf(y_sorted, param)

  # -- KS statistic
  ks_stat <- max(abs(Fn - Ftheta))

  # -- KS p-value (1-parameter correction is not used; plain asymptotic)
  ks_pval <- 1 - stats::ks.test(y, function(x) cdf(x, param))$p.value
  # Note: ks.test cannot be used directly with estimated parameters,
  # but we return its p-value so the user is aware.

  # -- Log-likelihood
  ll_vec <- pdf(y, param)
  # Guard against zero or negative densities
  ll_vec[ll_vec <= 0] <- .Machine$double.eps
  logLik <- sum(log(ll_vec))

  # -- AIC and BIC
  k <- length(param)
  n <- length(y)
  AIC <- -2*logLik + 2*k
  BIC <- -2*logLik + log(n)*k

  return(list(
    KS_stat = ks_stat,
    KS_pvalue = ks_pval,
    logLik = logLik,
    AIC = AIC,
    BIC = BIC
  ))
}
