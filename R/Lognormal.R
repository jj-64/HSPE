
# Lognormal sigma from Gini --------------
#' Lognormal Sigma From Gini
#'
#' Computes the lognormal sigma parameter from the Gini coefficient using:
#'   sigma = sqrt(2) * qnorm((1 + Gini)/2)
#'
#' @param Gini Numeric in (0,1). The Gini coefficient.
#'
#' @return The lognormal sigma.
#' @export
#' @examples
#' Sigma_LN(0.36)
#' [1] 0.661426
Sigma_LN <- function(Gini) {
  sqrt(2) * qnorm((1 + Gini) / 2)
}


## Lognormal meanlog (mu) -------------
#' Lognormal Meanlog (mu)
#'
#' For a lognormal distribution:
#'   mean = exp(mu + sigma^2 / 2)
#' so
#'   mu = log(mean) - sigma^2/2
#'
#' @param mean_y numerical value, observed mean of the variable
#' @param sigma numerical value, lognormal scale positive parameter
#'
#' @return mu (meanlog)
#' @export
#' @examples
#' mu_LN(1000, 0.36)
#' [1] 6.8429
mu_LN <- function(mean_y, sigma) {
  log(mean_y) - sigma^2 / 2
}


## Standard Error for Sigma parameter ------------
#' Standard Error of Lognormal Sigma From Gini
#'
#' Uses delta method: SE(sigma) = |d sigma / d G| * SE(Gini)
#'
#' @param mean_y numerical value, observed mean of the variable. Not used (kept for consistency).
#' @param Gini Numeric in (0,1). The Gini coefficient.
#' @param se_Gini Standard error of the Gini index.
#'
#' @return positive value, Standard error of sigma
#' @export
#' @examples
#' se_sigma_LN(1000,0.36, 0.01)
#' [1] 0.01977
se_sigma_LN <- function(mean_y, Gini, se_Gini) {
  if( Gini <=0 | Gini >=1){stop("Gini should be between 0 and 1")}
  if (mean_y <=0) stop ("mean_y should be positive")
   z <- (1 + Gini) / 2
  d_sigma_dG <- sqrt(2) * (1 / dnorm(qnorm(z))) * 0.5

  se_sigma_LN <- abs(d_sigma_dG) * se_Gini
  return(se_sigma_LN)
}

## Standard error for mu parameter -----------------
#' Standard Error of Lognormal mu
#'
#' Uses delta method based on partial derivatives wrt mean and Gini.
#'
#' @param mean_y numerical value, Mean of the variable.
#' @param Gini Numeric in (0,1). The Gini coefficient.
#' @param se_mean standard error of the mean
#' @param se_Gini standard error of the Gini index
#'
#' @return positive value, standard error of mu
#' @export
#' @examples
#' se_mu_LN(1000,0.36, 250, 0.01)
#' [1] 0.2503419
se_mu_LN <- function(mean_y, Gini, se_mean, se_Gini) {
  if( Gini <=0 | Gini >=1){stop("Gini should be between 0 and 1")}
  if (mean_y <=0) stop ("mean_y should be positive")
  LN_sigma <- Sigma_LN(Gini)

  # d mu / d mean
  d_mu_dmean <- 1 / mean_y

  # d sigma / d G
  z <- (1 + Gini) / 2
  d_sigma_dG <- sqrt(2) * (1 / dnorm(qnorm(z))) * 0.5

  # d mu / d G
  d_mu_dG <- -LN_sigma * d_sigma_dG

  se_mu_LN <- sqrt(
    (d_mu_dmean * se_mean)^2 +
      (d_mu_dG * se_Gini)^2
  )
  return(se_mu_LN)
}

## Compute the Lorenz Curve - Lopez ----------------
#' Lognormal Lorenz Curve
#'
#' Computes L(p) = Φ(Φ⁻¹(p) − sigma)
#'
#' @param P Quantiles in [0,1]
#' @param sigma numerical value, lognormal scale positive parameter
#'
#' @return Lorenz curve values.
#' @export
#' @examples
#' Lorenz_LN(Gini = 0.36)
#' # [1] 0.02600944 0.06641343 0.11784541 0.18015538 0.25416958 0.34160788 0.44550533 0.57150036
#' # [9] 0.73241244 1.00000000
Lorenz_LN <- function(p = seq(0.1,1,0.1), Gini) {
  if(Gini <=0 | Gini >=1) stop("sigma should be between 0 and 1")
  sigma = Sigma_LN(Gini)

  pnorm(qnorm(p) - sigma)
}


## PDF of Lognormal, same as dlnorm -----------
#' Lognormal PDF (fast)
#'
#' @param y numerical vector of the variable
#' @param mu numerical value, lognormal location parameter
#' @param sigma numerical value, lognormal scale positive parameter
#'
#' @return Density, probability between 0 and 1
#' @export
#' @examples
#' pdf_LN(100, 6, 12)
#' [1] 0.0003302136
#' pdf_LN(c(100, 200), 6, 12)
#' [1] 0.0003302136 0.0001659420
pdf_LN <- function(y, mu, sigma) {
  exp(-(log(y) - mu)^2 / (2 * sigma^2)) / (y * sigma * sqrt(2*pi))
}


## CDF of Lognormal, same as plnorm ------------
#' Lognormal CDF Expressed Using (mean, sigma)
#'
#' @param y numerical vector for the mean income
#' @param mean_y numerical value, average income
#' @param sigma numerical value, lognormal scale positive parameter
#'
#' @return CDF, probability between 0 and 1
#' @export
#' @examples
#' CDF_LN(1500, mean_y = 6000, 0.6)
#' [1] 0.02218964
#' CDF_LN(c(1500,7500), mean_y = 6000, 0.6)
#' [1] 0.02218964 0.74917820
CDF_LN <- function(y, mean_y, sigma) {
  if(sigma <=0) stop("sigma should be positive")
  if(mean_y <=0) stop("mean_y should be positive")
  pnorm((log(y / mean_y) / sigma) + sigma/2)
}

## Function for the Lorenz Curve for the country --------------
Lorenz_LN_fromIncome = function(y, mu, sigma) { ## the income group y, mean income and sigma
  if(sigma <=0) stop("sigma should be positive")
  if(mean_y <=0) stop("mean_y should be positive")
  ## vector for the intergrals

  integrand <- function(z) {z * dlnorm(z, mu, sigma) }
  int = integrate(integrand, lower = 0, upper = y)$value

  ## Formula
  Lorenz = exp(-(mu + sigma^2/2)) * int
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

## Headcount Standard Error ---------------
#' Standard Error of Lognormal Headcount Ratio
#'
#' @param p Poverty line
#' @param mean_y numerical value, observed mean of the variable
#' @param sigma numerical value, lognormal scale positive parameter
#' @param se_mean standard error of the mean
#' @param se_sigma standard error of sigma parameter
#' @param cov_mu_sigma Covariance of (mu, sigma), default 0
#'
#' @return standard error of the headcount
#' @export
#' @examples
#' HC_se_LN(900, mean_y = 2000, sigma = 0.36, se_mean = 250, se_sigma = 0.01)
#' 0.01767621 # 1.76%
HC_se_LN <- function(pov_line, mean_y, sigma, se_mean, se_sigma, cov_mu_sigma = 0) {

  z <- (log(pov_line / mean_y) / sigma) + sigma/2
  phi_z <- dnorm(z)

  dH_dmean  <- phi_z * (-1 / (mean_y * sigma))
  dH_dsigma <- phi_z * (-log(pov_line / mean_y) / sigma^2 + 0.5)

  var_H <- (dH_dmean^2)  * se_mean^2 +
    (dH_dsigma^2) * se_sigma^2 +
    2*dH_dmean*dH_dsigma*cov_mu_sigma

  sqrt(var_H)
}
