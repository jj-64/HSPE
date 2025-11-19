
library(numDeriv)

# CDF_B2 ------
#' @title CDF of the Beta Prime (Beta-2) Distribution
#'
#' @description
#' Computes the cumulative distribution function (CDF) of a **Beta Prime
#' (Beta-2)** distribution:
#'
#' \deqn{F(y \mid p, q, b) = I_{\,y/(b+y)}(p, q)}
#'
#' where \eqn{I_z(p,q)} is the regularized incomplete beta function.
#'
#' @param y Numeric vector of evaluation points.
#' @param b Scale parameter \eqn{b > 0}.
#' @param p First shape parameter \eqn{p > 0}.
#' @param q Second shape parameter \eqn{q > 0}.
#'
#' @return Numeric vector of CDF values.
#'
#' @examples
#' CDF_B2(y = 1:5, b = 2, p = 3, q = 4)
#'
#' @export
CDF_B2 <- function(y, b, p, q) {
  # params: named vector or list with p, q, b
  z <- y / (b + y)
  # regularized incomplete beta
  pbeta(z, p, q)
}

## Lorenz -------------------
#' @title Lorenz Curve of the Beta Prime (Beta-2) Distribution
#'
#' @description
#' Computes the Lorenz curve for the **Beta-2 (Beta Prime)** distribution.
#'
#' The Lorenz curve is:
#'
#' \deqn{
#' L(u) = \frac{ B\!\left( Q(u; p, q);\, p+1,\, q-1 \right) }{ B(p, q) }
#' }
#'
#' where:
#' - \eqn{Q(u)} is the quantile function (inverse CDF)
#' - \eqn{B(\cdot,\cdot)} is the beta function
#' - \eqn{u \in (0,1)}
#'
#' @param u Numeric vector in (0,1); cumulative population shares.
#' @param b Scale parameter (unused but included for consistency).
#' @param p First shape parameter.
#' @param q Second shape parameter (must satisfy \eqn{q > 1}).
#'
#' @return Numeric vector of Lorenz curve values.
#'
#' @examples
#' u <- seq(0.01, 0.99, length.out = 5)
#' Lorenz_B2(u, b = 2, p = 3, q = 4)
#'
#' @export
Lorenz_B2 <- function(u, b, p, q) {
  x <- qbeta(u, p, q)
  num <- pbeta(x, p + 1, q - 1)
  den <- beta(p, q)
  L <- num / den
  return(L)
}

# Delta-method SE for Beta-2 headcount using numerical gradient -------------
#' @title Standard Errors for Beta-2 Poverty Headcount via Delta Method
#'
#' @description
#' Computes **delta-method standard errors** for poverty headcount measures
#' based on the CDF of the **Beta-2 (Beta Prime)** distribution.
#'
#' For a poverty line \eqn{z}, the headcount is:
#'
#' \deqn{
#' H(z) = F(z \mid p, q, b)
#' }
#'
#' and the delta-method variance is:
#'
#' \deqn{
#' \mathrm{Var}(H) \approx \nabla_\theta H^\top \Sigma_\theta \nabla_\theta H
#' }
#'
#' where \eqn{\theta = (p, q, b)}.
#'
#' @param y Numeric vector of poverty lines.
#' @param b Scale parameter estimate.
#' @param p First shape parameter estimate.
#' @param q Second shape parameter estimate.
#' @param se_b Standard error of \code{b}.
#' @param se_p Standard error of \code{p}.
#' @param se_q Standard error of \code{q}.
#' @param eps Step size for numerical gradient (default 1e-6).
#'
#' @return A numeric vector of standard errors (same length as \code{y}).
#'
#' @details
#' The function:
#' - constructs a diagonal covariance matrix if only SEs are provided
#' - uses Richardson numerical differentiation
#' - returns one SE per poverty line
#'
#' @examples
#' y <- c(1, 2, 3)
#' HC_se_B2(
#'   y, b = 2, p = 3, q = 4,
#'   se_b = 0.1, se_p = 0.05, se_q = 0.06
#' )
#'
#' @export
HC_se_B2 <- function(y, b, p, q, se_b, se_p, se_q, eps = 1e-6) {
  if (!is.numeric(y)) stop("y must be numeric.")

  # Build covariance matrix (independence assumed)
  Sigma <- diag(as.numeric(c(se_p, se_q, se_b))^2)
  rownames(Sigma) <- colnames(Sigma) <- c("p", "q", "b")

  scalar_cdf <- function(theta, y_scalar) {
    names(theta) <- c("p", "q", "b")
    pbeta(y_scalar / (theta["b"] + y_scalar), theta["p"], theta["q"])
  }

  se_out <- numeric(length(y))
  theta0 <- c(p = p, q = q, b = b)

  for (i in seq_along(y)) {
    yi <- y[i]
    g <- numDeriv::grad(
      func = function(th) scalar_cdf(th, yi),
      x = theta0,
      method = "Richardson"
    )
    var_h <- as.numeric(t(g) %*% Sigma %*% g)
    if (var_h < 1e-14) var_h <- 0
    se_out[i] <- sqrt(var_h)
  }

  se_out
}
