
# Lognormal sigma from Gini --------------
Sigma_LN<- function(Gini) {
  sqrt(2) * qnorm((1 + Gini) / 2)
}

## Lognormal meanlog (mu) -------------
mu_LN = function(mean_y, Sigma){
  log(mean_y) - Sigma^2/2
}

## Standard Error for Sigma parameter ------------
se_sigma_LN = function(mean_y, Gini, se_Gini){

  sigma_LN = Sigma_LN(Gini)
  mu_LN    <- mu_LN(mean_y, sigma_LN)

  # Derivative dσ/dG
  z <- (1 + Gini) / 2
  d_sigma_dG <- sqrt(2) * (1 / dnorm(qnorm(z))) * (1/2)

  se_sigma_LN <- abs(d_sigma_dG) * se_Gini
  return(se_sigma_LN)
}

## Standard error for mu parameter -----------------
se_mu_LN = function(mean_y, Gini, se_mean, se_Gini){

  sigma_LN = Sigma_LN(Gini)

  # Derivative dμ/dmean
  d_mu_dmean <- 1 / mean_y

  # Derivative dσ/dG
  z <- (1 + Gini) / 2
  d_sigma_dG <- sqrt(2) * (1 / dnorm(qnorm(z))) * (1/2)

  # Derivative dμ/dG = -sigma * dσ/dG
  d_mu_dG <- -sigma_LN * d_sigma_dG

  se_mu_LN <- sqrt(
    (d_mu_dmean * se_mean)^2 +
      (d_mu_dG * se_Gini)^2
  )
  return(se_mu_LN)
}

## Compute the Lorenz Curve - Lopez ----------------
Lorenz_LN = function(Sigma, P) {

  Lorenz=0
  for(i in 1:length(P)) {
    Lorenz[i] = pnorm(qnorm(P[i]) - Sigma)
  }
  return(Lorenz)
}

## PDF of Lognormal, same as dlnorm -----------
pdf_LN = function(y,lmu,std) { ## lmu is mean in logscale
  exp(-(log(y) - lmu )^2/( 2*std^2) ) / (y*std*sqrt(2*pi))
}

## Function for the Lorenz Curve for the country --------------
Lorenz_LN = function(y, mu, Sigma) { ## the income group y, mean income and sigma

  ## vector for the intergrals

  integrand <- function(z) {z * dlnorm(z, mu, Sigma) }
  int = integrate(integrand, lower = 0, upper = y)$value

  ## Formula
  Lorenz = exp(-(mu + Sigma^2/2)) * int
  return(Lorenz)
}

## Function for the Lorenz Curve for the region ---------------
R_Lorenz_LN = function(y, R_mu, popw, mu, sigma) { ## the income group y, regional mean income and vector of population weight, vector of mu, of sigma

  ## vector for the intergrals
  int = 0

  for( i in 1:length(mu)) {
    integrand <- function(z) {z * dlnorm(z, mu[i], sigma[i]) }
    int[i] = integrate(integrand, lower = 0, upper = y)$value
  }

  ## Formula
  Lorenz = (1/R_mu) * sum (popw * int)
  return(Lorenz)
}
