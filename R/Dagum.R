## CDF ------------------
#' @title CDF of the Dagum Distribution
#'
#' @description
#' Computes the Dagum CDF:
#'
#' \deqn{
#' F(y \mid a, b, p)
#'   = \left( 1 + \left( \frac{y}{b} \right)^{-a} \right)^{-p}
#' }
#'
#' Equivalent to \code{VGAM::pdagum()} with
#' \code{shape1.a = a}, \code{shape2.p = p}, \code{scale = b}.
#'
#' @param y Numeric vector of evaluation points.
#' @param a Shape parameter \eqn{a > 0}.
#' @param b Scale parameter \eqn{b > 0}.
#' @param p Shape parameter \eqn{p > 0}.
#'
#' @return Numeric vector of CDF values.
#'
#' @examples
#' CDF_DA(1:5, a = 2, b = 1, p = 3)
#'
#' @export
CDF_DA <- function(y, a, b, p) {
  (1 + (y / b)^(-a))^(-p)
#VGAM::pdagum(q=y, shape1.a = a, shap2.p = p, scale=b)
}

## Lorenz-----------------------
#' @title Lorenz Curve of the Dagum Distribution
#'
#' @description
#' Computes the Lorenz curve of the Dagum distribution using the closed-form:
#'
#' \deqn{
#' L(u)
#'   = \frac{
#'       B\!\left( u^{1/p};\; p + 1/a,\; 1 - 1/a \right)
#'     }{
#'       B(p,\; 1 - 1/a)
#'     }
#' }
#'
#' where \eqn{B(\cdot,\cdot)} is the Beta function and
#' \eqn{u \in (0,1)} is the cumulative population share.
#'
#' @param u Numeric vector in \eqn{(0,1)}.
#' @param a Shape parameter of the Dagum distribution.
#' @param b Scale parameter (unused here but kept for consistency).
#' @param p Shape parameter.
#'
#' @return Numeric vector of Lorenz curve values.
#'
#' @examples
#' u <- seq(0.1, 0.9, length.out = 5)
#' Lorenz_DA(u, a = 2, b = 1, p = 3)
#'
#' @export
Lorenz_DA <- function(u = seq(0.1, 1, by = 0.1), a, b, p) {
  x <- u^(1 / p)
  num <- pbeta(x, p + 1 / a, 1 - 1 / a)
  den <- beta(p, 1 - 1 / a)
  L <- num / den
  return(L)
}

## HC Standard error -------------------
#' @title Standard Errors for Dagum Poverty Headcount via Delta Method
#'
#' @description
#' Computes **delta-method standard errors** for headcount poverty indices
#' using the Dagum CDF:
#'
#' \deqn{
#' H(z) = F(z \mid a, b, p)
#'   = \left( 1 + (z/b)^{-a} \right)^{-p}
#' }
#'
#' The delta-method approximation is:
#'
#' \deqn{
#' \mathrm{Var}(H)
#'   \approx
#'   \nabla_\theta H^\top
#'   \Sigma_\theta
#'   \nabla_\theta H,
#' }
#'
#' with parameters \eqn{\theta = (a, b, p)}.
#'
#' @param y Numeric vector of poverty lines.
#' @param a Estimated Dagum shape parameter.
#' @param b Estimated Dagum scale parameter.
#' @param p Estimated Dagum shape parameter.
#' @param se_a Standard error of \code{a}.
#' @param se_b Standard error of \code{b}.
#' @param se_p Standard error of \code{p}.
#'
#' @return Numeric vector of standard errors, same length as \code{y}.
#'
#' @details
#' - This implementation assumes **independent parameter estimates**
#'   (diagonal covariance matrix).
#' - If a full vcov matrix is available, the function can be extended easily.
#' - Computes the analytical gradient of the Dagum CDF wrt \eqn{a,b,p}.
#'
#' @examples
#' y <- c(1, 2, 3)
#' HC_se_DA(
#'   y,
#'   a = 2, b = 1, p = 3,
#'   se_a = 0.05, se_b = 0.02, se_p = 0.03
#' )
#'
#' @export
HC_se_DA <- function(y, a, b, p, se_a, se_b, se_p) {
  # y: poverty line(s)
  # a, b, p: Dagum parameters
  # vcov: 3x3 variance–covariance matrix from fitgroup.da

  y <- as.numeric(y)

  # compute u and F (CDF)
  u <- (y / b)^(-a)
  Fval <- (1 + u)^(-1 / p)

  # partial derivatives
  dF_da <- Fval * (u / (p * (1 + u))) * log(y / b)
  dF_db <- Fval * (a / (p * b)) * (u / (1 + u))
  dF_dp <- Fval * (log(1 + u) / p^2)

  # gradient matrix: each row corresponds to one y
  G <- cbind(dF_da, dF_db, dF_dp)

  # delta-method variance for each poverty line
  vcov = diag(as.numeric(c(se_a, se_b, se_p))^2)
  var_H <- rowSums((G %*% vcov) * G)

  # numerical safety
  var_H[var_H < 0 & var_H > -1e-12] <- 0
  var_H[var_H < 0] <- NA

  sqrt(var_H)
}
