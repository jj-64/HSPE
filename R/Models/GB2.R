# ---- CDF for GB2 ----
#' @title CDF of the Generalized Beta of the Second Kind (GB2)
#'
#' @description
#' Computes the CDF of the GB2 distribution:
#'
#' \deqn{
#' F(y \mid a, b, p, q)
#'   = B\!\left( \frac{(y/b)^a}{1 + (y/b)^a};\; p,\; q \right)
#' }
#'
#' where \eqn{B(\cdot;\,p,q)} is the regularized incomplete beta function.
#'
#' @param y Numeric vector of evaluation points.
#' @param a Shape parameter \eqn{a > 0}.
#' @param b Scale parameter \eqn{b > 0}.
#' @param p Shape parameter \eqn{p > 0}.
#' @param q Shape parameter \eqn{q > 0}.
#'
#' @return Numeric vector of CDF values.
#'
#' @examples
#' CDF_GB2(1:5, a = 2, b = 1, p = 3, q = 4)
#'
#' @export
CDF_GB2 <- function(y, a, b, p, q) {
  z <- (y / b)^a
  u <- z / (1 + z)
  pbeta(u, p, q)
}

## Lorenz -----------------
#' @title Lorenz Curve of the GB2 Distribution
#'
#' @description
#' Computes the Lorenz curve of the GB2 distribution using the identity:
#'
#' \deqn{
#' L(u)
#'   = \frac{
#'     B\!\left( x;\; p + 1/a,\; q - 1/a \right)
#'   }{
#'     B(p,\; q)
#'   },
#' }
#'
#' where \eqn{x = Q_{\mathrm{Beta}}(u;\, p, q)} is the beta quantile.
#'
#' @param u Numeric vector in \eqn{(0,1)}.
#' @param a GB2 shape parameter.
#' @param b GB2 scale parameter (unused here but kept for consistency).
#' @param p GB2 shape parameter.
#' @param q GB2 shape parameter.
#'
#' @return Numeric vector of Lorenz values.
#'
#' @examples
#' u <- seq(0.1, 0.9, length.out = 5)
#' Lorenz_GB2(u, a = 2, b = 1, p = 3, q = 4)
#'
#' @export
Lorenz_GB2 <- function(u, a, b, p, q) {
  x <- qbeta(u, p, q)
  num <- pbeta(x, p + 1 / a, q - 1 / a)
  den <- beta(p, q)
  num / den
  L <- num / den
  return(L)
}

# ---- Gradient of CDF wrt parameters (numerically) ----
#' @title Numerical Gradient of GB2 CDF
#'
#' @description
#' Computes a numerical gradient of the GB2 CDF with respect to
#' \eqn{(a, b, p, q)} using forward finite differences.
#'
#' \deqn{
#' g_j(y)
#'   \approx
#'   \frac{F(y \mid \theta_j + \varepsilon) - F(y \mid \theta)}{\varepsilon}
#' }
#'
#' @param y Numeric vector of evaluation points.
#' @param a,b,p,q GB2 parameters.
#' @param eps Numeric finite difference step.
#'
#' @return A matrix of size \code{length(y) × 4}, one gradient per parameter.
#'
#' @examples
#' grad_CDF_GB2(1:3, a = 2, b = 1, p = 3, q = 4)
#'
#' @export
grad_CDF_GB2 <- function(y, a, b, p, q, eps = 1e-6) {
  param <- c(a = a, b = b, p = p, q = q)
  f0 <- CDF_GB2(y, a, b, p, q)

  g <- matrix(0, nrow = length(y), ncol = 4)
  colnames(g) <- names(param)

  for (nm in names(param)) {
    param_eps <- param
    param_eps[nm] <- param_eps[nm] + eps

    f1 <- CDF_GB2(
      y,
      a = param_eps["a"],
      b = param_eps["b"],
      p = param_eps["p"],
      q = param_eps["q"]
    )

    g[, nm] <- (f1 - f0) / eps
  }

  g
}

# ---- Delta-method SE of headcount ----
#' @title Delta-Method Standard Errors for GB2 Poverty Headcount
#'
#' @description
#' Computes delta-method standard errors for the poverty headcount
#' \eqn{H(z) = F(z \mid a, b, p, q)} under the GB2 distribution.
#'
#' With gradient \eqn{\nabla_\theta H} and parameter covariance \eqn{\Sigma}:
#'
#' \deqn{
#' \mathrm{Var}(H)
#'   \approx
#'   \nabla_\theta H^\top
#'   \Sigma
#'   \nabla_\theta H.
#' }
#'
#' @param y Numeric vector of poverty lines.
#' @param a,b,p,q GB2 parameter estimates.
#' @param se_a,se_b,se_p,se_q Standard errors of \eqn{a,b,p,q}.
#'
#' @return Numeric vector of standard errors (same length as \code{y}).
#'
#' @details
#' - Assumes independence of parameter estimates (diagonal covariance).
#' - Uses numerical gradient from \code{grad_CDF_GB2()}.
#'
#' @examples
#' HC_se_GB2(
#'   y = c(1, 2, 3),
#'   a = 2, b = 1, p = 3, q = 4,
#'   se_a = 0.05, se_b = 0.03, se_p = 0.04, se_q = 0.02
#' )
#'
#' @export
HC_se_GB2 <- function(y, a, b, p, q, se_a, se_b, se_p, se_q) {
  g <- grad_CDF_GB2(y, a = a, b = b, p = p, q = q)

  vcov <- diag(c(se_a, se_b, se_p, se_q)^2)

  # Compute variance for each y: diag(G Σ Gᵀ)
  var_vec <- rowSums((g %*% vcov) * g)

  # numeric safety
  var_vec[var_vec < 0 & var_vec > -1e-12] <- 0
  var_vec[var_vec < 0] <- NA_real_

  sqrt(var_vec)
}
