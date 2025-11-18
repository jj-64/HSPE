
## α > 0 scale
## β >  0  shape

## CDF of Fisk -----------
#' Fisk (Log-Logistic) CDF
#'
#' Computes the cumulative distribution function of the Fisk distribution
#' with scale \code{scale > 0} and shape \code{shape > 0}.
#'
#' @param y Numeric vector of evaluation points.
#' @param scale Positive scale parameter.
#' @param shape Positive shape parameter.
#'
#' @return Numeric vector of CDF values.
#'
#' @details
#' The Fisk / Log-logistic CDF is:
#' \deqn{ F(y) = \frac{1}{1 + (y/scale)^{-shape}} }
#'
#' @export
#' @examples
#' CDF_Fisk(1:5, scale = 2, shape = 3)
CDF_Fisk <- function(y, scale, shape) {
  if (scale <= 0 || shape <= 0)
    stop("scale and shape must be > 0")

  1 / (1 + (y / scale)^(-shape))
}

## Pdf of Fisk ------------
#' Fisk (Log-Logistic) PDF
#'
#' Computes the probability density function of a Fisk distribution.
#'
#' @param y Numeric vector.
#' @param scale Positive scale parameter.
#' @param shape Positive shape parameter.
#'
#' @return Numeric vector of density values.
#'
#' @details
#' The Fisk PDF is:
#' \deqn{
#' f(y) = \frac{shape}{scale}
#'        \left( \frac{y}{scale} \right)^{shape - 1}
#'        \left[ 1 + \left( \frac{y}{scale} \right)^{shape} \right]^{-2}
#' }
#'
#' @export
#' @examples
#' pdf_Fisk(1:5, scale = 2, shape = 3)
pdf_Fisk <- function(y, scale, shape) {
  if (scale <= 0 || shape <= 0)
    stop("scale and shape must be > 0")

  (shape / scale) *
    (y / scale)^(shape - 1) *
    (1 + (y / scale)^shape)^(-2)
}

## function to compute the mean --------------
#' Mean of the Fisk Distribution
#'
#' @param scale Scale parameter (> 0).
#' @param shape Shape parameter (> 1). Mean exists only for shape > 1.
#'
#' @return Mean value, or NA if the mean is undefined.
#'
#' @export
#' @examples
#' mean_Fisk(scale = 2, shape = 3)
mean_Fisk <- function(shape, scale) {
  if (shape <= 1)
    return(NA_real_)  # Mean does not exist

  scale * gamma(1 + 1/shape) * gamma(1 - 1/shape)
}

## Gini --------------
#' Gini Coefficient of the Fisk Distribution
#'
#' @param shape Shape parameter (> 0).
#'
#' @return Gini coefficient.
#'
#' @export
#' @examples
#' Gini_Fisk(3)
Gini_Fisk <- function(shape) {
  if (shape <= 0)
    stop("shape must be > 0")
  1 / shape
}

#### Fisk mean-matching scale (alpha, a) parameter --------------
#' Mean-Matching Scale Parameter for Fisk Distribution
#'
#' @param shape Shape parameter.
#' @param mean_y numerical value, observed mean of the variable.
#'
#' @return Scale parameter consistent with the specified mean.
#' @export
#' @examples
#' scale_Fisk (3, 4000)
scale_Fisk <- function(shape, mean_y) {
  if (shape <= 1)
    stop("Mean exists only if shape > 1")

  mean_y / (gamma(1 + 1/shape) * gamma(1 - 1/shape))
}

#### Fisk mean-matching shape (beta, b) parameter --------------
#' Mean-Matching Shape Parameter for Fisk Distribution
#'
#' @param Gini Numeric in (0,1). The Gini coefficient
#'
#' @return Shape parameter.
#' @export
#' @examples
#' shape_Fisk (0.367)
shape_Fisk <- function(Gini) {
  if (Gini <= 0 || Gini >= 1)
    stop("Gini must be in (0,1)")

  1 / Gini
}

#######  Standard Error of Fisk Shape Parameter--------------
#' Standard Error of Fisk Shape Parameter
#'
#' Uses the delta method: \eqn{shape = 1/Gini}.
#'
#' @param Gini Numeric in (0,1). The Gini coefficient
#' @param se_Gini Standard error of the Gini estimate.
#'
#' @return Standard error of the Fisk shape parameter.
#' @export
#' @examples
#' se_shape_Fisk (0.36, 0.012)
se_shape_Fisk <- function(Gini, se_Gini) {
  # Derivative dbeta/dG
  d_b_dG <- -1 / (Gini^2)

  se_shape_Fisk <- abs(d_b_dG) * se_Gini
  return(se_shape_Fisk)
}

#######  Standard Error of a point estimate --------------
#' Standard Error of Fisk Scale Parameter (Delta Method)
#'
#' @param mean_y numerical value, observed mean of the variable.
#' @param Gini Numeric in (0,1). The Gini coefficient
#' @param se_mean Standard error of mean.
#' @param se_Gini Standard error of Gini.
#' @param shape_Fisk Optional shape value (otherwise computed from Gini).
#'
#' @return Standard error of the Fisk scale parameter.
#' @export
#' @examples
#' se_scale_Fisk (4000,0.36, 750,0.012)
se_scale_Fisk <- function(mean_y, Gini, se_mean, se_Gini,
                          shape_Fisk = NA)
{
  if (is.na(shape_Fisk)) shape_Fisk <- 1 / Gini

  # derivative d(shape)/d(Gini)
  d_b_dG <- -1 / (Gini^2)

  # derivative d(scale)/d(mean)
  h_b <- gamma(1 + 1 / shape_Fisk) * gamma(1 - 1 / shape_Fisk)
  d_a_dmean <- 1 / h_b

  # derivative d(scale)/d(shape) (numeric)
  eps <- 1e-6
  a_plus  <- scale_Fisk(shape_Fisk + eps, mean_y)
  a_minus <- scale_Fisk(shape_Fisk - eps, mean_y)
  d_a_db  <- (a_plus - a_minus) / (2 * eps)

  # chain rule d(scale)/d(Gini)  db/dG = db/da * da/dG
  d_a_dG <- d_a_db * d_b_dG

  se_shape_Fisk <- sqrt((d_a_dmean * se_mean)^2 +
         (d_a_dG    * se_Gini)^2)
  return(se_shape_Fisk)
}


## Compute the Lorenz Curve - FISK ---------------
#' Lorenz Curve of the Fisk Distribution
#'
#' @param shape Shape parameter (> 1).
#' @param P Vector of cumulative population shares in (0,1).
#'
#' @return Vector of Lorenz curve values L(P).
#'
#' @details
#' The Lorenz curve is:
#' \deqn{
#' L(p) = B(p; 1 + 1/shape,\; 1 - 1/shape)
#' }
#' where B is the incomplete beta CDF.
#'
#' @export
#' @examples
#' Lorenz_Fisk(p= seq(0.1, 1, by=0.1), shape = 3)
Lorenz_Fisk <- function(p, shape) {
  if (shape <= 1)
    stop("Lorenz curve exists only if shape > 1")

  pbeta(p, 1 + 1/shape, 1 - 1/shape)
}

### Headcount SE ------------------
#' Standard Error of Headcount Ratio Under Fisk Model
#'
#' Uses the delta method for:
#' \deqn{ HC = 1 - F(x) }
#'
#' @param x Poverty line.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @param se_shape SE of shape.
#' @param se_scale SE of scale.
#' @param cov_shape_scale Optional covariance term.
#'
#' @return Standard error of the headcount ratio.
#'
#' @export
#' @examples
#' HC_se_Fisk(2, shape = 3, scale = 2, se_shape = 0.1, se_scale = 0.2)
HC_se_Fisk <- function(x, shape, scale, se_shape, se_scale,
                       cov_shape_scale = 0)
{
  # compute u = (p/scale)^shape
  u <- (x / scale)^shape

  # derivatives of HC = 1 - F(x)
  dH_db <- (1 + u)^(-2) * u * log(x / scale)         # wrt shape
  dH_da <- (1 + u)^(-2) * u * (-shape / scale)       # wrt scale

  # delta-method variance
  var_H <- dH_da^2 * se_scale^2 +
    dH_db^2 * se_shape^2 +
    2 * dH_da * dH_db * cov_shape_scale

  sqrt(var_H)
}

## fit_fisk_param function ------------------
#' Fast Parameter Estimation for Fisk Distribution Using Only (Mean, Gini)
#'
#' Computes the Fisk (log-logistic) parameters using moment-matching:
#'   shape = 1/Gini
#'   scale = mean / [Gamma(1+1/shape) * Gamma(1-1/shape)]
#'
#' Optionally propagates standard errors using the delta method.
#'
#' @param mean_y numerical value, observed mean of the variable.
#' @param Gini Numeric in (0,1). The Gini coefficient
#' @param se_mean Standard error of mean (optional).
#' @param se_Gini Standard error of Gini (optional).
#'
#' @return A list with:
#'   \describe{
#'     \item{par}{Named vector: c(shape, scale)}
#'     \item{se}{Named vector: c(se_shape, se_scale)}
#'   }
#' @export
#' @examples
#' fit <- fit_Fisk_param(mean_y= 20000, Gini = 0.32)
#' fit
#' $par
#' shape     scale
#' 3.125 16797.370
#'
#' $se
#' se_shape se_scale
#' NA       NA

fit_Fisk_param <- function(mean_y, Gini, se_mean = NA, se_Gini = NA) {

  # --- validation ---
  if (Gini <= 0 || Gini >= 1)
    stop("Gini must be in (0,1).")

  shape <- shape_Fisk(Gini)

  if (shape <= 1)
    stop("For Fisk, the mean exists only if shape > 1 (i.e. Gini < 1).")

  if (mean_y <= 0)
    stop("mean_y must be > 0.")

  # --- parameter estimates ---
  scale <- scale_Fisk(shape, mean_y)

  # --- standard errors (optional) ---
  if (!is.na(se_mean) && !is.na(se_Gini)) {

    se_shape <- se_shape_Fisk(Gini, se_Gini)
    se_scale <- se_scale_Fisk(
      mean_y = mean_y,
      Gini   = Gini,
      se_mean = se_mean,
      se_Gini = se_Gini,
      shape_Fisk = shape
    )

  } else {
    se_shape <- NA
    se_scale <- NA
  }

  list(
    par = c(shape = shape, scale = scale),
    se  = c(se_shape = se_shape, se_scale = se_scale)
  )
}


## fit fisk -------------------------
#' Loss function for Fisk (Log-Logistic) distribution
#'
#' @param params c(shape, scale)
#' @param pc.inc numerical value, average income
#' @param gini.e Numeric in (0,1). The Gini coefficient
#' @param y Observed non-cumulative income shares
#' @param p_dec Probabilities (e.g. seq(.1,.9,length.out=length(y)+1))
#'
#' @return Scalar loss value
loss_Fisk <- function(params, pc.inc, gini.e, y,
                      p_dec = seq(0.1, 0.9, length.out = length(y) + 1)) {

  shape <- params[1]
  scale <- params[2]

  # invalid parameter region
  if (shape <= 1 || scale <= 0) return(1e12)

  # theoretical mean
  mean_th <- mean_Fisk(shape, scale)
  mean_err <- (mean_th - pc.inc)^2

  # theoretical Gini
  gini_th <- Gini_Fisk(shape)
  gini_err <- (gini_th - gini.e)^2

  # theoretical Lorenz → deciles
  L_th <- Lorenz_Fisk(shape, P = c(0, p_dec))
  dec_th <- diff(L_th)
  dec_err <- sum((dec_th - y)^2)

  # weights
  w1 <- 1  # mean
  w2 <- 1  # gini
  w3 <- 10 # deciles

  w1*mean_err + w2*gini_err + w3*dec_err
}

#' Fit the Fisk Distribution from Grouped Data
#'
#' @param y Non-cumulative income shares (deciles, vigintiles, etc.)
#' @param pc.inc numerical value, average income
#' @param gini.e Numeric in (0,1). The Gini coefficient
#' @param N Sample size (optional; used only for future SE estimation)
#' @param se.ewmd Logical: return EWMD bootstrap SEs (placeholder)
#' @param se.scale Logical: return SEs for scale (placeholder)
#' @param nrep Number of bootstrap replications (future extension)
#' @param control A list of optimizer controls passed to \code{optim}
#'
#' @return A list:
#' \describe{
#'   \item{par}{c(shape, scale)}
#'   \item{value}{Loss function value}
#'   \item{convergence}{0 = ok}
#'   \item{gini}{Model-implied Gini}
#'   \item{mean}{Model-implied mean}
#' }
#'
#' @examples
#' fit <- fit_Fisk(y = deciles, pc.inc = 20000, gini.e = 0.32)
#' fit$par
fit_Fisk <- function(y, pc.inc, gini.e, N = NULL,
                     se.ewmd = TRUE, se.scale = TRUE,
                     nrep = 1000,
                     control = list(maxit = 500))
{
  p_dec <- seq(0.1, 0.9, length.out = length(y) + 1)

  # starting values based on moment-matching
  shape0 <- shape_Fisk(gini.e)              # 1/Gini
  scale0 <- scale_Fisk(shape0, pc.inc)

  out <- optim(
    par = c(shape = shape0, scale = scale0),
    fn = loss_Fisk,
    pc.inc = pc.inc,
    gini.e = gini.e,
    y = y,
    p_dec = p_dec,
    method = "L-BFGS-B",
    lower = c(1.01, 1e-6),
    upper = c(100, 1e6),
    control = control
  )

  # compute implied quantities
  shape_hat <- out$par[1]
  scale_hat <- out$par[2]

  list(
    par = c(shape = shape_hat, scale = scale_hat),
    value = out$value,
    convergence = out$convergence,
    mean = mean_Fisk(shape_hat, scale_hat),
    gini = Gini_Fisk(shape_hat)
  )
}
