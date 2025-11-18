CDF_SM <- function(y, a, b, q)
{ VGAM::psinmad(q=y, scale = b , shape1.a = a, shape3.q = q) }


HC_se_SM <- function(pov_line, a, b, q, se_a, se_b, se_q) {
  # pov_line: poverty line(s)
  # a, b, q: Singh-Maddala parameters
  # vcov: 3x3 variance–covariance matrix from fitgroup.sm

  pov_line <- as.numeric(pov_line)

  # components
  v <- (pov_line / b)^a
  u <- 1 + v
  Fval <- 1 - u^(-q)

  # partial derivatives
  dF_da <-  q * u^(-q - 1) * v * log(pov_line / b)
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
