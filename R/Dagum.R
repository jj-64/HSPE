CDF_DA <- function(y, a, b, p)
  { (1 + (y / b)^(-a))^(-p) }
#VGAM::pdagum(q=y, shape1.a = a, shap2.p = p, scale=b)

Lorenz_DA <- function(u = seq(0.1, 1, by =0.1), a, b, p) {
  # u is between o and 1
  x  <- u^(1/p)
  num <- pbeta(x, p + 1/a, 1 - 1/a)
  den <- beta(p, 1 - 1/a)
  L <- num / den
  return(L)
}

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
