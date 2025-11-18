# ---- CDF for GB2 ----
CDF_GB2 <- function(y, a, b, p ,q){
  z <- (y/b)^a
  u <- z / (1 + z)
  pbeta(u, p, q)
}

# ---- Gradient of CDF wrt parameters (numerically) ----
grad_CDF_GB2 <- function(y, a, b, p ,q, eps = 1e-6){
  f0 <- CDF_GB2(y, a=a, b=b, p=p, q=q)
  param = c(a, b, p, q)
  g  <- numeric(4)
  names(g) <- c("a","b","p","q")

  for (nm in names(param)) {
    param_eps <- param
    param_eps[[nm]] <- param[[nm]] + eps
    g[nm] <- (CDF_GB2(y, a=a, b=b, p=p, q=q) - f0) / eps
  }
  g
}

# ---- Delta-method SE of headcount ----
HC_se_GB2 <- function(y, a, b ,p ,q , se_a, se_b, se_p, se_q){
  g <- grad_CDF_GB2(y, a=a, b=b, p=p, q=q)          # 4×1 gradient
  vcov = diag(c(se_a, se_b, se_p, se_q)^2)
  sqrt( as.numeric(t(g) %*% vcov %*% g) )
}
