# ---- CDF for GB2 ----
CDF_GB2 <- function(y, a, b, p ,q){
  p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
  z <- (y/b)^a
  u <- z / (1 + z)
  pbeta(u, p, q)
}

# ---- Gradient of CDF wrt parameters (numerically) ----
grad_CDF_GB2 <- function(y, a, b, p ,q, eps = 1e-6) {
  param <- c(a=a, b=b, p=p, q=q)
  f0 <- CDF_GB2(y, a, b, p, q)

  g <- matrix(0, nrow = length(y), ncol = 4)
  colnames(g) <- names(param)

  for (nm in names(param)) {
    param_eps <- param
    param_eps[nm] <- param_eps[nm] + eps
    f1 <- CDF_GB2(y, a=as.numeric(param_eps["a"]), b=as.numeric(param_eps["b"]),
                  p=as.numeric(param_eps["p"]), q=as.numeric(param_eps["q"]))
    g[, nm] <- (f1 - f0) / eps
  }

  g
}

# ---- Delta-method SE of headcount ----
HC_se_GB2 <- function(y, a, b ,p ,q , se_a, se_b, se_p, se_q){
  g <- grad_CDF_GB2(y, a=a, b=b, p=p, q=q)          # 4×1 gradient
  vcov = diag(c(se_a, se_b, se_p, se_q)^2)
  sqrt( as.numeric(g %*% vcov %*% t(g)) )
}
