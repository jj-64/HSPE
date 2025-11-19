CDF_SM <- function(y, a, b, q)
{ VGAM::psinmad(q=y, scale = b , shape1.a = a, shape3.q = q) }

Lorenz_SM <- function(p = seq(0.1,1,by =0.1), a, b, q) {
  # Lorenz of Singh-Maddala
  term <- (1 - p)^(1 - 1/a)
  num  <- pbeta(1 - p, q - 1/a, 1 + 1/a)
  den  <- beta(q, 1 + 1/a)
  L <- 1 - term * num / den
  return(L)
}

HC_se_SM <- function(y, a, b, q, se_a, se_b, se_q) {
  # y: poverty line(s)
  # a, b, q: Singh-Maddala parameters
  # vcov: 3x3 variance–covariance matrix from fitgroup.sm

  y <- as.numeric(y)

  # components
  v <- (y / b)^a
  u <- 1 + v
  Fval <- 1 - u^(-q)

  # partial derivatives
  dF_da <-  q * u^(-q - 1) * v * log(y / b)
  dF_db <- -q * u^(-q - 1) * (a / b) * v
  dF_dq <-  u^(-q) * log(u)

  # gradient matrix
  vcov = diag(as.numeric(c(se_a, se_b, se_q))^2)
  G <- cbind(dF_da, dF_db, dF_dq)

  # delta-method variance
  var_H <- rowSums((G %*% vcov) * G)

  # numerical safety
  var_H[var_H < 0 & var_H > -1e-12] <- 0
  var_H[var_H < 0] <- NA

  sqrt(var_H)
}
