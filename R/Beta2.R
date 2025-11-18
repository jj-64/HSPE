# Install numDeriv if not already installed
# install.packages("numDeriv")

library(numDeriv)

# CDF_B2 as you defined (helper)
CDF_B2 <- function(y, b, p , q) {
  # params: named vector or list with p, q, b
  z <- y / (b + y)
  # regularized incomplete beta
  pbeta(z, p, q)
}

# Delta-method SE for Beta-2 headcount using numerical gradient
HC_se_B2 <- function(y, b, p, q, se_b, se_p, se_q,
                     eps = 1e-6) {
  # y: numeric vector of poverty lines
  #
  # Returns numeric vector of SEs (same length as y), on same scale as CDF (0..1).

  # --- Input checks ---
  if (!is.numeric(y)) stop("y must be numeric.")

  # Build variance-covariance matrix
    # assume independence if only SEs provided
    Sigma <- diag(as.numeric(c(se_p, se_q, se_b))^2)

  # ensure order is p,q,b in Sigma rows/cols
  # If Sigma has dimnames, try to reorder; otherwise assume order matches.
  if (!is.null(rownames(Sigma))) {
    desired <- c("p","q","b")
    if (!all(desired %in% rownames(Sigma))) {
      warning("vcov provided but does not have rownames p,q,b; assuming column order is (p,q,b).")
    } else {
      Sigma <- Sigma[desired, desired]
    }
  }

  # function to compute scalar CDF for one y and parameter vector (numeric)
  scalar_cdf <- function(theta, y_scalar) {
    names(theta) <- c("p","q","b")
    pbeta( y_scalar / (theta["b"] + y_scalar), theta["p"], theta["q"] )
  }

  # compute SE for each y by numerical gradient
  se_out <- numeric(length(y))
  theta0 <- c(p = as.numeric(p), q = as.numeric(q), b = as.numeric(b))

  for (i in seq_along(y)) {
    yi <- y[i]
    # gradient vector g = d CDF / d theta at theta0
    g <- tryCatch({
      grad(func = function(th) scalar_cdf(th, yi), x = theta0, method = "Richardson")
    }, error = function(e) {
      stop("numerical gradient failed: ", e$message)
    })

    # delta variance
    var_h <- as.numeric( t(g) %*% Sigma %*% g )
    # numerical safety
    if (var_h < 0 && var_h > -1e-12) var_h <- 0
    if (var_h < 0) var_h <- NA_real_
    se_out[i] <- sqrt(var_h)
  }

  se_out
}
