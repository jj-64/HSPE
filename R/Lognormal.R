
# Lognormal sigma from Gini --------------
Sigma_LN<- function(Gini) {
  sqrt(2) * qnorm((1 + Gini) / 2)
}

## Lognormal meanlog (mu) -------------
mu_LN = function(mean_y, LN_sigma){
  log(mean_y) - LN_sigma^2/2
}

## Standard Error for Sigma parameter ------------
se_sigma_LN = function(mean_y, Gini, se_Gini){

  #LN_sigma = Sigma_LN(Gini)
  #mu_LN    <- mu_LN(mean_y, LN_sigma)

  # Derivative dσ/dG
  z <- (1 + Gini) / 2
  d_sigma_dG <- sqrt(2) * (1/dnorm(qnorm(z))) * (1/2)

  se_sigma_LN <- abs(d_sigma_dG) * se_Gini
  return(se_sigma_LN)
}

## Standard error for mu parameter -----------------
se_mu_LN = function(mean_y, Gini, se_mean, se_Gini){

  LN_sigma = Sigma_LN(Gini)

  # Derivative dμ/dmean
  d_mu_dmean <- 1 / mean_y

  # Derivative dσ/dG
  z <- (1 + Gini) / 2
  d_sigma_dG <- sqrt(2) * (1 / dnorm(qnorm(z))) * (1/2)

  # Derivative dμ/dG = -sigma * dσ/dG
  d_mu_dG <- -LN_sigma * d_sigma_dG

  se_mu_LN <- sqrt(
    (d_mu_dmean * se_mean)^2 +
      (d_mu_dG * se_Gini)^2
  )
  return(se_mu_LN)
}

## Compute the Lorenz Curve - Lopez ----------------
Lorenz_LN = function(LN_sigma, P) {

  Lorenz=0
  for(i in 1:length(P)) {
    Lorenz[i] = pnorm(qnorm(P[i]) - LN_sigma)
  }
  return(Lorenz)
}

## PDF of Lognormal, same as dlnorm -----------
pdf_LN = function(y,mu_LN,LN_sigma) { ## lmu is mean in logscale
  exp(-(log(y) - mu_LN )^2/( 2*LN_sigma^2) ) / (y*LN_sigma*sqrt(2*pi))
}

## CDF of Lognormal, same as plnorm ------------
CDF_LN = function(y, mean_y, LN_sigma){
  #pnorm((log(y) - mu_LN)/LN_sigma)
  pnorm( (log(y / mean_y) / LN_sigma) + LN_sigma / 2 )
}

## Function for the Lorenz Curve for the country --------------
Lorenz_LN = function(y, mu_LN, LN_sigma) { ## the income group y, mean income and sigma

  ## vector for the intergrals

  integrand <- function(z) {z * dlnorm(z, mu_LN, LN_sigma) }
  int = integrate(integrand, lower = 0, upper = y)$value

  ## Formula
  Lorenz = exp(-(mu_LN + LN_sigma^2/2)) * int
  return(Lorenz)
}

## Function for the Lorenz Curve for the region ---------------
R_Lorenz_LN = function(y, R_mu, popw, mu_LN, sigma) { ## the income group y, regional mean income and vector of population weight, vector of mu_LN, of sigma

  ## vector for the intergrals
  int = 0

  for( i in 1:length(mu_LN)) {
    integrand <- function(z) {z * dlnorm(z, mu_LN[i], sigma[i]) }
    int[i] = integrate(integrand, lower = 0, upper = y)$value
  }

  ## Formula
  Lorenz = (1/R_mu) * sum (popw * int)
  return(Lorenz)
}

## Headcount Standard Error ---------------
HC_se_LN <- function(p, mean_y, LN_sigma, se_mean, se_sigma, se_mu, cov_mu_sigma = 0) {

  #LN_loc = mu_LN(mean_y, LN_sigma)

  # compute z value
  z <- (log(p / mean_y) / LN_sigma) + LN_sigma / 2

  # normal pdf
  phi_z <- dnorm(z)

  # derivatives
  dH_dmean    <- phi_z * (-1 / (mean_y * LN_sigma))
  dH_dsigma <- phi_z * ( -log(p / mean_y) / (LN_sigma^2) + 0.5 )
  #dH_dmu <- phi_z * (-1 / (mean_y * LN_sigma)) * exp(LN_sigma^2/2 + LN_loc)

  # delta-method variance
  var_H <- (dH_dmean^2)    * se_Average^2 +
    (dH_dsigma^2) * se_sigma^2 +
    2 * dH_dmean * dH_dsigma * cov_mu_sigma

  # var_H <- (dH_dmu^2)    * se_mu^2 +
  #   (dH_dsigma^2) * se_sigma^2 +
  #   2 * dH_dmean * dH_dsigma * cov_mu_sigma
  sqrt(var_H)  # return SE
}
