## Grouped data -----------------
#diagnostics_grouped(L_obs = L_values, Lorenz_fun = Lorenz_FISK, par = 2, N=N)
diagnostics_grouped <- function(L_obs, Lorenz_fun, par, N) {
  # L_obs: observed Lorenz cumulative shares at p
  # Lorenz_fun(p, par): model Lorenz curve
  # par: parameter list
  # N: sample size (for grouped likelihood)

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
  # Prob. mass in each bin = ΔCDF
  # But we approximate ΔCDF using ΔLorenz * (mean * population)
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
#pdf_FISK_param = function(y, param){pdf_FISK(y=y, scale= param$a, shape = param$b)}
#diagonistics_dist(y = y, cdf = CDF_FISK_param, pdf = pdf_FISK_param, param = as.list(c(a=1, b=2)))
diagnostics_micro <- function(y, cdf, pdf, param) {
  # y: data vector
  # cdf: function(y, param) returning vector of CDF values
  # pdf: function(y, param) returning vector of PDF values
  # param: list of named vector of parameters

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
