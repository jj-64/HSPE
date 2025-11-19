#install.packages("hypergeo")
library(hypergeo)
library(nloptr)
# packages
if (!requireNamespace("numDeriv", quietly = TRUE)) {
  install.packages("numDeriv")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(numDeriv)
library(dplyr)

## alpha is shape
## beta is scale
## Estimate shape via optimization ---------------
#' Estimate New Pareto Shape Parameter
#'
#' Estimates the shape parameter \eqn{\alpha} of the New Pareto distribution
#' by solving:
#' \deqn{G(\alpha) = Gini}
#'
#' using derivative-free nonlinear optimization (COBYLA).
#'
#' @param Gini The target Gini coefficient in (0,1).
#'
#' @return Estimated shape parameter \eqn{\alpha}.
#'
#' @details
#' The optimizer minimizes:
#' \deqn{|G_NP(\alpha) - Gini|}
#'
#' subject to:
#' \deqn{\alpha \ge 1.1.}
#'
#' @import nloptr
#' @export
#' @examples
#' shape_NP(0.36)
#' #[1] 2.128961
shape_NP <- function(Gini) {

  if (Gini <= 0 || Gini >= 1)
    stop("Gini must be in (0,1).")

    ## Objective function
    eval_f0 <- function( shape,a,b){
      return(abs(Gini_NP(shape[1]) - Gini))
    }

    # constraint function g(x) <=0
    eval_g0 <- function(shape, a,b) {
      return(b- a*shape[1] )
    }

    # define parameters
    a <- c(1)
    b <- c(1)

    # Solve using NLOPT_LN_COBYLA without gradient information
    res <- nloptr( x0=c(1.1),
                    eval_f=eval_f0,
                    lb = c(1.1),
                    ub = c(20),
                    eval_g_ineq = eval_g0,
                    opts = list("algorithm"="NLOPT_LN_COBYLA",
                                "xtol_rel"=1.0e-8),
                    a = a,
                    b=b)


  return(as.numeric(res$solution))

}

# Wrapper for shape as function of G (so derivative uses numeric differentiation)
#' Safe Wrapper for shape_NP()
#'
#' @param Gini Gini coefficient.
#'
#' @return scalar numeric shape or NA.
shape_NP_func_wrapper <- function(Gini) {
  # call your shape_NP implementation; return scalar numeric
  # add tryCatch because shape_NP uses optimizer and may fail
  res <- tryCatch({
    as.numeric(shape_NP(Gini))
  }, error = function(e) NA_real_)
  res
}

## Numerical derivative for α(G) ------------------
#' Numerical Derivative of Shape wrt Gini
#'
#' @param Gini Gini coefficient.
#' @param eps Step size.
#'
#' @return Numeric derivative d(alpha)/d(G).
shape_NP_deriv <- function(Gini){
  eps <- 1e-6
  (shape_NP(Gini + eps) - shape_NP(Gini - eps)) / (2 * eps)
}

## Estimate scale from shape and mean by Optimization -------------
#' Estimate New Pareto Scale Parameter
#'
#' Solves for \eqn{\theta} from:
#' \deqn{E[Y] = \frac{2\alpha\theta \, \mathrm{Re}\{\,_2F_1(2,2-(1+\alpha)/\alpha;3-(1+\alpha)/\alpha;-1)\}}{\alpha - 1}}
#'
#' @param mean_y Mean income.
#' @param shape Previously estimated shape parameter.
#'
#' @return Estimated scale parameter.
#'
#' @import nloptr
#' @importFrom hypergeo hypergeo
#' @export
#' @examples
#' scale_NP(mean_y= 8000, shape= 0.4)
#' 515.4424
scale_NP = function(mean_y , shape){
  if (mean_y <= 0) stop("mean_y must be > 0")
  if (shape <= 1)  stop("shape must be > 1 for finite mean")
  ## Objective function
  eval_f0 <- function( scale,a,b){
    return( abs((2*shape *scale *Re(hypergeo(2,2-((1+shape)/shape),3-((1+shape)/shape),-1))/(shape - 1)) - mean_y ) )
  }

  # constraint function g(x) <=0
  eval_g0 <- function(scale, a,b) {
    return(b- a*scale[1] )
  }

  # define parameters
  a <- c(1)
  b <- c(0)

  # Solve using NLOPT_LN_COBYLA without gradient information
  res1 <- nloptr( x0=c(0.1),
                  eval_f=eval_f0,
                  lb = c(0.1),
                  ub = c(Inf),
                  eval_g_ineq = eval_g0,
                  opts = list("algorithm"="NLOPT_LN_COBYLA",
                              "xtol_rel"=1.0e-8),
                  a = a,
                  b=b)

  return(res1$solution)
}

# Wrapper for scale as function of (mean_y, shape)
#' Safe Wrapper for scale_NP()
#'
#' @param x Numeric vector: c(mean, shape)
#'
#' @return scalar scale or NA
scale_NP_func_wrapper <- function(x) {
  # x is numeric vector: c(mean_y, shape)
  mean_y <- x[1]; shape <- x[2]
  res <- tryCatch({
    as.numeric(scale_NP(mean_y, shape))
  }, error = function(e) NA_real_)
  res
}

## Standard error for shape --------------------------
# Delta-method function
# mean_y: sample mean
# Gini: sample Gini
# se_mean: SE(sample mean)
# se_Gini: SE(sample Gini)
# cov_mean_gini: optional covariance between mean and gini (default 0)
# eps: finite-diff step
#' Standard Error for New Pareto Shape
#'
#' @param Gini Gini coefficient.
#' @param se_Gini SE(Gini).
#' @param shape Optional precomputed shape.
#' @param eps Finite difference step.
#'
#' @return Standard error of shape.
#'
#' @importFrom numDeriv grad
#' @export
#' @examples
#' se_shape_NP(Gini =0.36, se_Gini = 0.1, shape= 0.3)
#' #[1] 0.4810901
se_shape_NP <- function(Gini, se_Gini, shape =NA, eps = 1e-6) {

  # 1) shape (alpha) and its derivative w.r.t Gini
  if(is.na(shape)) shape <- shape_NP_func_wrapper(Gini)
  if (is.na(shape)) stop("shape_NP returned NA or failed at provided Gini")

  d_shape_dG <- num_deriv(function(G) shape_NP_func_wrapper(G), Gini, eps = eps)

  var_shape <- (d_shape_dG^2) * (se_Gini^2)
  se_shape  <- sqrt(var_shape)

  return(se_shape)

}

## Standard error for scale --------------------------
# Delta-method function
# mean_y: sample mean
# Gini: sample Gini
# se_mean: SE(sample mean)
# se_Gini: SE(sample Gini)
# cov_mean_gini: optional covariance between mean and gini (default 0)
# eps: finite-diff step
#' Standard Error for New Pareto Scale Parameter
#'
#' @param mean_y Mean income.
#' @param Gini Gini coefficient.
#' @param se_mean SE(mean).
#' @param se_Gini SE(Gini).
#' @param shape Optional shape (cached).
#' @param se_shape Optional SE(shape).
#' @param scale Optional scale (cached).
#' @param cov_mean_gini Cov(mean, Gini).
#' @param eps finite-difference step.
#'
#' @return Standard error of scale.
#'
#' @importFrom numDeriv grad
#' @export
#' @examples
#' se_scale_NP(mean_y = 6000, se_mean = 250, Gini =0.36, se_Gini = 0.1, shape= 0.3, se_shape = 0.48, scale= 150)
#' #[1] 0.
se_scale_NP <- function(mean_y, Gini, se_mean, se_Gini, shape=NA, se_shape = NA, scale= NA,
                       cov_mean_gini = 0, eps = 1e-6) {
  # 1) shape and its derivative w.r.t Gini
  if(is.na(shape)) shape <- shape_NP_func_wrapper(Gini)
  if (is.na(shape)) stop("shape_NP returned NA or failed at provided Gini")
  if(is.na(se_shape)) se_shape_NP(Gini, se_Gini, shape, eps= eps)

  d_shape_dG <- num_deriv(function(G) shape_NP_func_wrapper(G), Gini, eps = eps)

  # 2) scale and partial derivatives wrt mean_y and shape
  if(is.na(scale)) scale <- scale_NP_func_wrapper(c(mean_y, shape))
  if (is.na(scale)) stop("scale_NP returned NA or failed at provided mean_y & shape")

  # partial derivative wrt mean_y
  d_scale_dMean  <- num_deriv(function(x) scale_NP_func_wrapper(c(x, shape)), mean_y, eps = eps)
  # partial derivative wrt shape (hold mean_y fixed)
  d_scale_dAlpha <- num_deriv(function(a) scale_NP_func_wrapper(c(mean_y, a)), shape, eps = eps)

  # variance/covariance matrix of (mean_y, shape)
  var_mean  <- se_mean^2
  var_shape <- se_shape^2
  cov_ma    <- cov_mean_gini * d_shape_dG  # if you provided Cov(mean, Gini), transform to cov(mean,shape)
  # Explanation: cov(mean, shape) ≈ (d shape / dG) * cov(mean, G)
  # If user provided cov_mean_gini already as Cov(mean,shape), then they should pass that value.
  # We'll accept both: if cov_mean_gini is small relative to var_mean & var_shape, it's effectively 0.

  # If user probably passed Cov(mean,Gini), cov_ma computed above; otherwise assume they passed cov(mean,shape)
  # We'll check a heuristic: if cov_mean_gini magnitude is > 1e-8 and seems like Cov(mean,Gini) (rare),
  # still use transformation. If user passed cov(mean,shape) directly, they can pass cov_mean_gini already transformed.

  # Build gradient and covariance
  grad <- c(d_scale_dMean, d_scale_dAlpha)
  Sigma <- matrix(c(var_mean, cov_ma, cov_ma, var_shape), nrow = 2)

  var_scale <- as.numeric(t(grad) %*% Sigma %*% grad)
  se_scale  <- sqrt(max(0, var_scale))

  return(se_scale)
}

## Lorenz Curve for the country -------------
CDF_NP_integral = function(y, shape, scale) { ## the income group y, shape and scale of the NPareto, mean income

  ## vector for the intergrals

  integrand <- function(z) {z * pdf_NP(z=z, shape=shape, scale=scale) }
  int= c()
  for(i in seq_along(y)) int[i] = integrate(integrand, lower = scale, upper = y[i])$value

  ## Formula
  Lorenz = 1/mean(y) * int
  return(Lorenz)
}

Lorenz_NP_Exact <- function(p = seq(0.1, 1, by =0.1), shape){

  if (shape <= 0) stop("shape must be positive")

  doubleInt <- function(p,shape){
    f <- function(t,shape){((1+t)/(1-t))^(1/shape)}
    fvalue <- integrate(f,lower=0,upper=p,shape=shape)$value
    return(fvalue)
  }

  integrationTrap <- function (func, lower, upper, shape){
    x = seq(lower,upper,by=0.001)  ## divide the x-axis into small intervals
    n = length(x)
    x <- as.matrix(x,n,1)
    Y <- matrix(0,n,1)
    Y <- as.matrix(Y,n,1)
    tra <- matrix(0,n,1)

    for(i in 1:n)
    {
      Y[i] <- func(x[i],shape)  ## the integral over evey portion/interval
    }

    tra[1] <- 0
    for (i in 2:n) {
      tra[i] <-((0.5)*(x[i]-x[i-1])*(Y[i]+Y[i-1]))  }
    return(sum(tra))
  }

  ## Gini integral
  denominator = Re(shape * hypergeo(1,-1/shape,2-1/shape,-1) / (shape-1) )

  L_0 = sapply(p, FUN = function(p) integrationTrap(func=doubleInt,lower=0, upper = p, shape=shape))

  L = L_0*100/denominator
  L[length(L)] = 100- sum(L[-length(L)])
  return(cumsum(L))
  }


# 1) Fast Lorenz (vectorized grid + interpolation)
#    Uses a non-uniform grid with clustering near 1 to capture integrand behaviour.
#    Notes / tuning
#   - n_grid: increase for more accuracy (default 3000). If you need more speed, reduce.
#   - cluster_power: greater -> more grid density near 1; adjust 2..5 as needed.
#   - Hessian: numeric Hessian inversion provides an approximate covariance; for more
#     accurate SEs you can bootstrap (resample grouped residuals) – more costly.
#' Fast Lorenz Curve for New Pareto Distribution
#'
#' Approximates:
#'   L(p) = ∫₀ᵖ ∫₀ᵗ ((1+u)/(1-u))^(1/shape) du dt / denom
#'
#' @param p Probabilities in [0,1].
#' @param shape Shape parameter.
#' @param n_grid Grid resolution.
#' @param cluster_power Grid clustering exponent.
#'
#' @return Lorenz curve values.
#' @export
#' @examples
#' Lorenz_NP(shape=3)
#' # [1] 0.05919209 0.12251535 0.19040316 0.26346712 0.34259232 0.42912127 0.52525375 0.63509513
#' # [9] 0.76851956 1.00000000
Lorenz_NP <- function(p = seq(0.1, 1, by = 0.1),
                      shape,
                      n_grid = 3000,
                      cluster_power = 3,
                      eps = 1e-12) {

  if (shape <= 0)
    stop("shape must be positive")


  # create non-uniform grid on [0, 1) clustering near 1
  u <- seq(0, 1, length.out = n_grid)
  t_grid <- 1 - (1 - u)^cluster_power

  # ensure last grid point slightly less than 1 to avoid blow-ups
  t_grid[length(t_grid)] <- 1 - 1e-10

  # inner integrand = ((1+t)/(1-t))^(1/shape)
  inner_f <- function(t)
    ((1 + t)/(1 - t))^(1/shape)


  # cumulative integral by trapezoid: cumulative integral from 0 to t_i
  # compute trapezoid increments
  fvals <- inner_f(t_grid)
  dx <- diff(t_grid)
  inc <- 0.5 * dx * (fvals[-1] + fvals[-n_grid])
  cumint <- c(0, cumsum(inc))

  ## denominator is cumulative integral at p = 1 (but that's limited from <1)
  # denom <- Re(shape * hypergeo::hypergeo(1,-1/shape,2-1/shape,-1) / (shape-1) ) ## exact
  denom <- cumint[n_grid]
  if (!is.finite(denom) || denom <= 0)
    stop("denominator not finite")

  # for requested p, simple interpolation of cumint vs t_grid
  # ensure p within [0, 1-eps]
  p_req <- p
  p_req[p_req >= 1] <- 1 - 1e-10
  p_req[p_req < 0] <- 0
  Lvals <- approx(x = t_grid, y = cumint, xout = p_req, rule = 2)$y / denom

  #p[p >= 1] <- 1 - 1e-10
  #approx(x = t_grid, y = cumint, xout = p)$y / denom
  return(Lvals)
}

## pdf of NewPareto --------------------
#' New Pareto PDF
#'
#' @param shape Shape parameter.
#' @param scale Scale (minimum income).
#' @param p Income values.
#'
#' @return Density values.
#' @export
#' @examples
#' pdf_NP(2000, shape= 2, scale=100)
#' #[1]0.000004975093
pdf_NP <- function(y, shape, scale) {
  ifelse(
    y >= scale,
    (2 * shape * scale^shape * y^(shape - 1)) /
      (y^shape + scale^shape)^2,
    0
  )
}

## CDF of NewPareto ----------------
#' New Pareto CDF
#'
#' @param y Quantile/income.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#'
#' @return CDF value.
#' @export
#' @examples
#' CDF_NP(2000, shape= 2, scale=100)
#' #[1] 0.9950125
CDF_NP <- function(y, shape, scale) {
  y <- as.numeric(y)
  ifelse(
    y < scale,
    0,
    1 - 2 / (1 + (y/scale)^shape)
  )
}

## Quantile Function from percentiles ------------------
#' New Pareto Quantile function
#'
#' @param p numerical vector, probability between 0 and  1.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#'
#' @return quantile value.
#' @export
#' @examples
#' Quantile_NP(0.2, shape= 2, scale=100)
#' #[1] 122.4745
Quantile_NP <- function(p, shape, scale) {
  if(p <=0 | p>=1) stop("p shoule be between 0 and 1")
  # p in (0,1). Solve 1 - 2/(1+(q/scale)^shape) = p => q = scale * ((2/(1-p)-1)^(1/shape))
  if (any(p >= 1 | p <= 0)) stop("p must be in (0,1)")
  scale * ((2 / (1 - p) - 1)^(1 / shape))
}

## Random Generated NP ----------------------
#NP distributed x: Given a uniform variate U drawn from U(0, 1) distribution, we obtain X, which is NP distributed and is given by
#' New Pareto Random generating function
#'
#' @param u numerical vector, uniform number between 0 and  1.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#'
#' @return NP value.
#' @export
#' @examples
#' random_NP(0.2, shape= 2, scale=100)
#' #[1] 300
random_NP <- function (u,shape, scale){(((2*(scale^shape))/u)-(scale^shape))^(1/shape)}

## Gini from NP distribution ------------------
Gini_NP_Exact <- function(shape){

  integ = Lorenz_NP(p=1, shape= shape)
  gini <- 1- 2*integ

  ## Return
  return(gini)
}

# 2) Gini from Lorenz (numeric integration)

Gini_NP <- function(shape, n_grid = 3000) {
  # compute Lorenz on grid and integrate
  p_grid <- seq(0, 1, length.out = n_grid)
  L_grid <- Lorenz_NP(p_grid, shape = shape, n_grid = n_grid)
  # integrate L(p) via trapezoid
  dp <- diff(p_grid)[1]
  integral_L <- sum(0.5 * (L_grid[-1] + L_grid[-length(L_grid)]) * dp)
  G <- 1 - 2 * integral_L
  return(G)
}

## Mean from NP distribution -----------------
Mean_NP_Exact <- function(shape, scale){
  if (shape <= 1) {
    return(Inf) # mean diverges for shape <= 1
  }

  2*shape *scale *Re(hypergeo::hypergeo(2,2-((1+shape)/shape),3-((1+shape)/shape),-1))/(shape - 1)
}

# 3) Mean via survival integral (numeric)
#    E[X] = scale * (1 + 2 * ∫_0^1 t^(shape-2) / (1 + t^shape) dt)
#    derived by substitution u = 1/x. Works for shape > 1.
Mean_NP <- function(shape, scale, reltol = 1e-8) {
  if (shape <= 1) {
    return(Inf) # mean diverges for shape <= 1
  }
  integrand <- function(t) {
    t^(shape - 2) / (1 + t^shape)
  }
  val <- tryCatch({
    stats::integrate(integrand, lower = 0, upper = 1, rel.tol = reltol)$value
  }, error = function(e) {
    # fallback to a robust trapezoid integration on log-space if integrate fails
    grid <- seq(0, 1, length.out = 2000)
    grid[1] <- 1e-12
    v <- integrand(grid)
    sum(0.5 * (v[-1] + v[-length(v)]) * diff(grid))
  })
  mean_val <- scale * (1 + 2 * val)
  return(mean_val)
}

## MLE estimate -----------------------------------
#Code from paper for maximum likelihood, returns shape and scale estimates
functionMLE_NP <- function(z){                 ## Code from paper for maximum likelihood, returns shape and scale estimates
  z <- sort(z)
  x <- z[2:length(z)]
  b <- min(z)
  shapeguess <- length(z)/sum(log(z)-log(b))

  g<- function(a,b,x=x){
    res <- -2*sum(((b/x)^a*log(b/x))/((b/x)^a + 1 )) + sum(log(b/x)) + length(x)/a
    return(res)
  }

  gd <- function(a,b,x=x){
    b1 <- -length(x)/a^2-2*sum(((b/x)^a *
                                  (log(b/x))^2)/((b/x)^a+1) -
                                 ((b/x)^(2*a)*(log(b/x))^2)/((b/x)^a+1)^2)
    return(b1)
  }

  mm  <- shapeguess
  tol <- 1e-8
  val <- 5
  xn  <- mm
  res <- 1
  cont <- 0

  while (val >tol){
    res <- (xn-(g(a=xn,b=b,x=x)/gd(a=xn,b=b,x=x)))
    val <- (abs(res-xn)/xn)
    xn  <- res
    cont <- cont +1
  }

  minimum <- min(z)
  converge <-c("FALSE")
  shapeResult <- res
  scaleResult <- b
  results <- list(shapeEstimate = shapeResult,
                  scaleEstimate = scaleResult)
  return(results)
}

## Method 1 of estimating micro data ####
# estimated_params <- functionMLE(fyy2)
#
# min_value <- min(fyy2)
# max_value <- max(fyy2)
#
# num_new_values <- length(fyy2)
# new_values <- numeric(num_new_values)

# Inverse transform sampling
# for (i in 1:num_new_values) {
#   u <- runif(1)  ## Generate a uniform random number
#   y <- min_value + u * (max_value - min_value)
#   while (runif(1) > Fy_NPareto(estimated_params$shapeEstimate, estimated_params$scaleEstimate, y)) {
#     u <- runif(1)
#     y <- min_value + u * (max_value - min_value)
#   }
#   new_values[i] <- y
# }


## Method 2 of estimating micro data ####
#
# new_values1 <- random_NP(runif(fyy2), estimated_params$shapeEstimate, estimated_params$scaleEstimate)
#
# sorted_new_values <- sort(new_values)
# sorted_new_values1 <- sort(new_values1)
#
# df <- data.frame(DataVector = sorted_new_values, NewValues = sorted_new_values1)


# Requires: nloptr (for your existing shape_NP / scale_NP), stats


## Bootstrap function (recommended if you have raw income vector) ---------------
# incomes: numeric vector of person-level incomes (or unit-level)
# R: number of bootstrap resamples (e.g., 1000)
# parallel: "no" or use multicore? (not implemented here)
se_newpareto_boot <- function(incomes, R = 1000, seed = 1234, shape_init = NULL) {
  set.seed(seed)
  n <- length(incomes)

  # storage
  shape_bs <- numeric(R)
  scale_bs  <- numeric(R)
  mean_bs  <- numeric(R)
  gini_bs  <- numeric(R)

  # function to compute Gini; uses ineq::Gini if available, otherwise custom
  compute_gini <- function(x) {
    if (requireNamespace("ineq", quietly = TRUE)) {
      return(ineq::Gini(x))
    } else {
      # simple Gini
      x <- sort(x)
      n <- length(x)
      index <- seq_len(n)
      (2 * sum(index * x) / (n * sum(x))) - (n + 1) / n
    }
  }

  pb <- txtProgressBar(min = 0, max = R, style = 3)
  for (r in seq_len(R)) {
    idx <- sample.int(n, n, replace = TRUE)
    samp <- incomes[idx]
    m_hat <- mean(samp)
    g_hat <- compute_gini(samp)

    # try shape and scale; wrap in tryCatch
    al <- tryCatch(shape_NP(g_hat), error = function(e) NA_real_)
    be <- NA_real_
    if (!is.na(al)) {
      be <- tryCatch(scale_NP(m_hat, al), error = function(e) NA_real_)
    }

    mean_bs[r] <- m_hat
    gini_bs[r] <- g_hat
    shape_bs[r] <- al
    scale_bs[r]  <- be

    setTxtProgressBar(pb, r)
  }
  close(pb)

  # Remove failed iterations
  ok <- is.finite(shape_bs) & is.finite(scale_bs)
  if (sum(ok) < max(10, 0.1 * R)) {
    warning("Many bootstrap iterations failed. Check shape_NP / scale_NP stability.")
  }

  shape_bs <- shape_bs[ok]
  scale_bs  <- scale_bs[ok]
  mean_bs  <- mean_bs[ok]
  gini_bs  <- gini_bs[ok]

  # Estimates
  shape <- mean(shape_bs)
  scale_hat  <- mean(scale_bs)

  se_shape <- sd(shape_bs)
  se_scale  <- sd(scale_bs)

  cov_mean_gini <- cov(mean_bs, gini_bs)
  cov_mean_shape <- cov(mean_bs, shape_bs)
  cov_shape_scale <- cov(shape_bs, scale_bs)

  list(
    shape = shape,
    se_shape = se_shape,
    scale_hat = scale_hat,
    se_scale = se_scale,
    cov_mean_gini = cov_mean_gini,
    cov_mean_shape = cov_mean_shape,
    cov_shape_scale = cov_shape_scale,
    raw = list(mean_bs = mean_bs, gini_bs = gini_bs, shape_bs = shape_bs, scale_bs = scale_bs)
  )
}

## Headcount SE ------------------
HC_se_NP <- function(y, shape, scale, se_shape, se_scale, cov_shape_scale = 0) {

  # compute u = (p/scale)^shape
  u <- (y / scale)^shape

  # HC derivative components
  dH_da <- (1 + u)^(-2) * 2 * u * log(y / scale)  ## derivative w.r.t shape
  dH_db <- (1 + u)^(-2) * 2* u * (-shape / scale) ## derivative w.r.t scale

  # delta-method variance
  var_H <- dH_db^2 * se_scale^2 +
    dH_da^2 * se_shape^2 +
    2 * dH_da * dH_db * cov_shape_scale

  sqrt(var_H)  # return standard error
}


## fitgroup.np fitgroup-like wrapper + SEs/CIs  -------------------------------

# a) EWMD-like objective (fit to grouped Lorenz points)
#    We'll match the observed grouped Lorenz (pc.inc = vector of p: .1 ... 1)
#    Optionally include Gini and mean in the objective via weights.
objective_np <- function(par, L_obs,
                         pc.inc = NULL, gini.e = NULL,
                         mean_weight = 1, gini_weight = 1, lorenz_weight = 10,
                         n_grid = 3000) {
  shape <- par[1]
  scale <- par[2]

  p_obs = seq(0.1, 1, by = (1/(length(L_obs))))  ##L_obs is cumulative Lorenz points

  # invalid penalties
  if (shape <= 1 || scale <= 0) return(1e12)

  # theoretical Lorenz at observed p
  L_th <- Lorenz_NP(p_obs, shape = shape, n_grid = n_grid)

  lorenz_err <- sum((L_th - L_obs)^2)

  mean_err <- 0
  if (!is.null(pc.inc)) {
    mean_th <- Mean_NP(shape, scale)
    mean_err <- (mean_th - pc.inc)^2
  }

  gini_err <- 0
  if (!is.null(gini.e)) {
    gini_th <- Gini_NP(shape)
    gini_err <- (gini_th - gini.e)^2
  }

  # weighted sum of squared errors
  loss <- lorenz_weight * lorenz_err + mean_weight * mean_err + gini_weight * gini_err
  return(loss)
}


# b) fitgroup-like wrapper: fitgroup.np()
#    Returns a list with ewmd.estimation ("Coef.", "se") and gini.estimation
# y <- income_shares_nonCum (Lorenz_NP(shape= 3))
# pc.inc <- Mean_NP(shape=3, scale= 1000)
# gini.e <- Gini_NP(shape = 3)
# fit_out <- fitgroup.np(y, pc.inc, gini.e)
# print(fit_out)
# fit_out$ewmd.estimation
# Coefs = fit_out$ewmd.estimation["Coef.",]
# ses =  fit_out$ewmd.estimation["se",]
# conf_bound( Coefs, ses )
fitgroup.np <- function(y,    # observed Non cumulative income shares at y (same length)
                               pc.inc = NULL,  # observed mean (optional, can be used in obj)
                               gini.e = NULL,  # observed gini (optional)
                               gini = TRUE,
                               N = NULL,
                               se.ewmd = TRUE,
                               se.scale = TRUE,
                               nrep = 20,
                               control = list() ) {


  # sanitize inputs
  L_obs = as.numeric(cumsum(y))
  #if (length(y) != length(L_obs)) stop("y and L_obs must have same length")
  #p_obs <- as.numeric(seq(0.1, 1, 1))   # vector of cumulative population shares (p, e.g., 0.1..1)

  # initial guesses:
  # - scale: use pc.inc/2 if mean provided, otherwise median quantile heuristic
  start_scale <- if (!is.null(pc.inc) && is.finite(pc.inc) && pc.inc > 0) {
    pc.inc / 2
  } else {
    1
  }
  start_shape <- 2

  # optimization
  opt <- tryCatch({
    optim(par = c(start_shape, start_scale),
          fn = objective_np,
          L_obs = L_obs,
          pc.inc = pc.inc,
          gini.e = gini.e,
          method = "L-BFGS-B",
          lower = c(1 + 1e-6, 1e-8),
          upper = c(500, 1e12),
          control = modifyList(list(maxit = 1000), control)
    )
  }, error = function(e) {
    list(par = c(NA, NA), value = Inf, converged = FALSE, message = e$message)
  })

  if (!isTRUE(opt$convergence == 0)) {
    warning("fitgroup.np: optimization may not have converged: ", opt$convergence)
  }

  est <- opt$par
  names(est) <- c("shape", "scale")

  # compute standard errors using numeric Hessian on the objective function
  se <- c(NA, NA)
  if (se.ewmd && is.finite(opt$value) && all(is.finite(est))) {
    # numeric Hessian of objective at optimum
    hess <- tryCatch({
      numDeriv::hessian(func = function(p) objective_np(p, L_obs = L_obs,
                                                        pc.inc = pc.inc, gini.e = gini.e),
                        x = est)
    }, error = function(e) NULL)
    ## Ensure symmetric hessian
    for(k in 1:nrow(hess)){
      for( l in 1:ncol(hess)) {
           if(l != k) hess[l,k] = hess[k,l]
    }}

    if (!is.null(hess) && all(is.finite(hess))) {
      # covariance approx = 2 * (H^{-1}) * s2 ? For least-squares-like objective, treat Hessian as 2 * J'J
      # We'll use cov = solve(hess), then se = sqrt(diag(cov))
      # But ensure hess invertible
      inv_h <- tryCatch(solve(hess), error = function(e) NULL)
      if(is.null(inv_h)) inv_h <- tryCatch(MASS::ginv(hess), error = function(e) NULL)
      if (!is.null(inv_h)) {
        se_val <- sqrt(pmax(0, diag(inv_h)))
        # For scaling we don't have residual variance estimate; we present these "raw" se estimates
        se <- se_val
        names(se) <- c("shape", "scale")
      }
    }
  }

  # compute gini and mean estimates for the fitted params
  gini_est <- tryCatch(Gini_NP(shape = est[1]), error = function(e) NA)
  mean_est <- tryCatch(Mean_NP(shape = est[1], scale = est[2]), error = function(e) NA)

  # produce results similar to fitgroup.* structures
  out <- list(
    ewmd.estimation = rbind("Coef." = est, "se" = se),
    gini.estimation = c(mean = as.numeric(mean_est), gini = gini_est),
    objective_value = opt$value,
    converged = (opt$convergence == 0),
    optim_output = opt
  )

  return(out)
}


## fitgroup.np_tidy and return tidy output with CIs ---------------------------------------------------------------
# y <- income_shares_nonCum (Lorenz_NP(shape= 3))
# pc.inc <- Mean_NP(shape=3, scale= 1000)
# gini.e <- Gini_NP(shape = 3)
# fit_out <- fitgroup.np_tidy("CountryX", y, pc.inc, gini.e)
# print(fit_out)
# attr(fit_out, "fit_object")$ewmd.estimation
fitgroup.np_tidy <- function(country_name,
                                y,
                               pc.inc = NULL, gini.e = NULL, N = NULL,
                               se.ewmd = TRUE, control = list()) {

  fit <- fitgroup.np(y = L_vec,
                            pc.inc = pc.inc, gini.e = gini.e, N = N,
                            se.ewmd = se.ewmd, control = control)

  coefs <- as.numeric(fit$ewmd.estimation["Coef.", ])
  names(coefs) <- colnames(fit$ewmd.estimation)
  ses <- as.numeric(fit$ewmd.estimation["se", ])

  ##
  ci_from_fit <- function(coefs, ses, level = 0.95) {
      z <- qnorm((1 + level) / 2)

      data.frame(
        param = names(coefs),
        est = coefs,
        se = ses,
        lower = coefs - z * ses,
        upper = coefs + z * ses,
        row.names = NULL
      )
    }

  ci_tab <- ci_from_fit(coefs, ses )

  tib <- tibble::tibble(
    Country = country_name,
    shape = coefs["shape"],
    scale = coefs["scale"],
    mean_est = fit$gini.estimation["mean"],
    gini_est = fit$gini.estimation["gini"]
  )

  tib_ci <- ci_tab %>%
    tidyr::pivot_wider(names_from = param,
                       values_from = c(est, se, lower, upper),
                       names_sep = "_")

  # merge and return
  out <- dplyr::bind_cols(tib, tib_ci)
  attr(out, "fit_object") <- fit
  return(out)
}



## fit NP----------------

#' Fit the New Pareto (NP) Income Distribution
#'
#' @param y Numeric vector of non-cumulative income shares.
#' @param pc.inc Observed mean income.
#' @param gini.e Observed Gini coefficient.
#' @param N Optional sample size (reserved for future SE estimation).
#' @param se.ewmd Logical; placeholder for compatibility.
#' @param se.scale Logical; placeholder for compatibility.
#' @param nrep Integer; placeholder for bootstrapping future SEs.
#' @param control List of optimizer controls (passed to \code{optim}).
#'
#' @return A list returned by \code{optim}.
#'
#' @details
#' The function minimizes \code{loss_NP} over parameters \code{a} and \code{b}
#' using L-BFGS-B with bounds ensuring valid NP parameters.
#'
#' @examples
#' \dontrun{
#' fit <- fit_np(y = deciles,
#'               pc.inc = 20000,
#'               gini.e = 0.32)
#' fit$par
#' }
fit_np <- function(y, pc.inc, gini.e, N = NULL,
                   se.ewmd = TRUE, se.scale = TRUE,
                   nrep = 1000,
                   control = list(maxit = 500))
{
  p_dec <- seq(0.1, 0.9, length.out = length(y) + 1)

  out <- optim(
    par    = c(a = 2, b = pc.inc / 2),
    fn     = loss_NP,
    pc.inc = pc.inc,
    gini.e = gini.e,
    y      = y,
    p_dec  = p_dec,
    method = "L-BFGS-B",
    lower  = c(1.01, 1e-6),
    upper  = c(50, 1e6),
    control = control
  )

  out
}

## fit_NP_param function ------------------
#' Fast Estimation of NewPareto Parameters Using Mean and Gini
#'
#' Computes:
#' shape as in \code{\link{shape_NP}}
#' scale  as in \code{\link{scale_NP}}
#'
#' Optionally computes delta-method standard errors as in
#'  \code{\link{se_shape_NP}} and  \code{\link{se_scale_NP}}
#'
#' @param mean_y Mean of the data (must be > 0)
#' @param Gini Gini coefficient in (0,1)
#' @param se_mean Standard error of mean (optional)
#' @param se_Gini Standard error of Gini (optional)
#'
#' @return A list with:
#'   \describe{
#'     \item{par}{c(shape, scale)}
#'     \item{se}{c(se_shape, se_scale)}
#'   }
#' @export
#' @examples
#' fit_NP_param(6000, 0.36, se_mean = 250, se_Gini = 0.01)
#' $par
#' # shape       scale
#' #  2.128961  2523.457703
#' #
#' # se_shape     se_scale
#' # 0.04810901 125.06197149
fit_NP_param <- function(mean_y, Gini, se_mean = NA, se_Gini = NA) {

  # --- domain checks ---
  if (mean_y <= 0) stop("mean_y must be > 0.")
  if (Gini <= 0 || Gini >= 1) stop("Gini must be in (0,1).")

  # parameters
  shape <- shape_NP(Gini)
  scale    <- scale_NP(mean_y, shape)

  # SEs
  if (!is.na(se_mean) && !is.na(se_Gini)) {
    se_shape <- se_shape_NP(Gini, se_Gini, shape)
    se_scale <- se_scale_NP(mean_y, Gini, se_mean, se_Gini, shape, se_shape)

  } else {
    se_shape <- NA
    se_scale    <- NA
  }

  list(
    par = c(shape = shape, scale = scale),
    se  = c(se_shape = se_shape , se_scale = se_scale)
  )
}
