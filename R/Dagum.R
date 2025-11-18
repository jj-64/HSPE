HC_se_DA <- function(pov_line, a, b, p, se_a, se_b, se_p) {
  # pov_line: poverty line(s)
  # a, b, p: Dagum parameters
  # vcov: 3x3 variance–covariance matrix from fitgroup.da

  pov_line <- as.numeric(pov_line)

  # compute u and F (CDF)
  u <- (pov_line / b)^(-a)
  Fval <- (1 + u)^(-1 / p)

  # partial derivatives
  dF_da <- Fval * (u / (p * (1 + u))) * log(pov_line / b)
  dF_db <- Fval * (a / (p * b)) * (u / (1 + u))
  dF_dp <- Fval * (log(1 + u) / p^2)

  # gradient matrix: each row corresponds to one y
  G <- cbind(dF_da, dF_db, dF_dp)

  # delta-method variance for each poverty line
  vcov = diag(c(se_a, se_b, se_p))
  var_H <- rowSums((G %*% vcov) * G)

  # numerical safety
  var_H[var_H < 0 & var_H > -1e-12] <- 0
  var_H[var_H < 0] <- NA

  sqrt(var_H)
}
