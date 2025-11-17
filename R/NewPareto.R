#install.packages("hypergeo")
library(hypergeo)
library(nloptr)

## alpha is shape
## beta is scale
## Estimate shape via optimization ---------------
shape_NP = function(Gini){
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
  res1 <- nloptr( x0=c(1.1),
                  eval_f=eval_f0,
                  lb = c(1.1),
                  ub = c(5),
                  eval_g_ineq = eval_g0,
                  opts = list("algorithm"="NLOPT_LN_COBYLA",
                              "xtol_rel"=1.0e-8),
                  a = a,
                  b=b)

  return(res1$solution)
}

# Wrapper for shape as function of G (so derivative uses numeric differentiation)
shape_NP_func_wrapper <- function(Gini) {
  # call your shape_NP implementation; return scalar numeric
  # add tryCatch because shape_NP uses optimizer and may fail
  res <- tryCatch({
    as.numeric(shape_NP(Gini))
  }, error = function(e) NA_real_)
  res
}

## Numerical derivative for α(G) ------------------
shape_NP_deriv <- function(Gini){
  eps <- 1e-6
  (shape_NP(Gini + eps) - shape_NP(Gini - eps)) / (2 * eps)
}

## Estimate scale from shape and mean by Optimization -------------
scale_NP = function(mean_y , shape){
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

## Function for the Lorenz Curve for the country -------------
Lorenz_NP_fromY = function(y, shape, scale) { ## the income group y, shape and scale of the NPareto, mean income

  ## vector for the intergrals

  integrand <- function(z) {z * pdf_NP(shape, scale, z) }
  int= c()
  for(i in seq_along(y)) int[i] = integrate(integrand, lower = scale, upper = y[i])$value

  ## Formula
  Lorenz = 1/mean(y) * int
  return(Lorenz)
}

Lorenz_NP <- function(p = seq(0.1, 1, by =0.1), shape){

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

  L = sapply(p, FUN = function(p) integrationTrap(func=doubleInt,lower=0, upper = p, shape=shape))

  return(L/denominator)
  }

## Quantile function to build Lorenz curve, returns the percentile ------------
Percentile_NP = function(shape, q) {

  doubleInt <- function(q,shape){
    f <- function(t,shape){((1+t)/(1-t))^(1/shape)}
    fvalue = c()
    for (i in seq_along(q)) fvalue[i] <- integrate(f,lower=0,upper=q[i],shape=shape)$value
    return(fvalue)
  }

  return ( doubleInt(q,shape)/  Re(shape * hypergeo::hypergeo(1,-1/shape,2-1/shape,-1) / (shape-1) ) )

}

## pdf of NewPareto --------------------
pdf_NP = function(shape,scale,p) {
  return(ifelse(p>= scale, (2 * shape * scale^shape *p^(shape-1) ) / (p^shape + scale^shape)^2 ,0))
}

## CDF of NewPareto ----------------
CDF_NP = function(shape,scale,q) {
  #return(ifelse(q>= scale, (q^shape - scale^shape) / (q^shape + scale^shape),0))
  return(ifelse(q>= scale, 1 -2 / ( 1+(q/scale)^shape),0))
}

## Quantile Function from percentiles ------------------
Quantile_NP = function(shape,scale,p) {
  #return(ifelse(p<1, scale * ((-p-1)/(p-1))^(1/shape),100000000))
  scale * ((2/(1-p) - 1)^(1/shape))
}

## Radom Generated NP ----------------------
#NP distributed x: Given a uniform variate U drawn from U(0, 1) distribution, we obtain X, which is NP distributed and is given by
random_NP <- function (u,a,b){(((2*(b^a))/u)-(b^a))^(1/a)}

## Gini from NP distribution ------------------
Gini_NP <- function(shape){

  # doubleInt <- function(p,shape){
  #   f <- function(t,shape){((1+t)/(1-t))^(1/shape)}
  #   fvalue <- integrate(f,lower=0,upper=p,shape=shape)$value
  #   return(fvalue)
  # }
  #
  # integrationTrap <- function (func, lower, upper, shape){
  #   x = seq(lower,upper,by=0.001)  ## divide the x-axis into small intervals
  #   n = length(x)
  #   x <- as.matrix(x,n,1)
  #   Y <- matrix(0,n,1)
  #   Y <- as.matrix(Y,n,1)
  #   tra <- matrix(0,n,1)
  #
  #   for(i in 1:n)
  #   {
  #     Y[i] <- func(x[i],shape)  ## the integral over evey portion/interval
  #   }
  #
  #   tra[1] <- 0
  #   for (i in 2:n) {
  #     tra[i] <-((0.5)*(x[i]-x[i-1])*(Y[i]+Y[i-1]))  }
  #   return(sum(tra))
  # }
  #
  # ## Gini integral
  # integ = integrationTrap(func=doubleInt,lower=0, upper = 1, shape=shape)/Re(shape * hypergeo(1,-1/shape,2-1/shape,-1) / (shape-1) )

  integ = Lorenz_NP(p=1, shape= shape)
  gini <- 1- 2*integ

  ## Return
  return(gini)
}

## Mean from NP distribution -----------------
mean_NP <- function(shape, scale){
  2*shape *scale *Re(hypergeo(2,2-((1+shape)/shape),3-((1+shape)/shape),-1))/(shape - 1)
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

### Method 1 of estimating micro data ####
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


### Method 2 of estimating micro data ####
#
# new_values1 <- random_NP(runif(fyy2), estimated_params$shapeEstimate, estimated_params$scaleEstimate)
#
# sorted_new_values <- sort(new_values)
# sorted_new_values1 <- sort(new_values1)
#
# df <- data.frame(DataVector = sorted_new_values, NewValues = sorted_new_values1)


# Requires: nloptr (for your existing shape_NP / scale_NP), stats


# ---------------------------
# Bootstrap function (recommended if you have raw income vector)
# ---------------------------
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

### Headcount SE ------------------
HC_se_NP <- function(x, shape, scale, se_shape, se_scale, cov_shape_scale = 0) {

  # compute u = (p/scale)^shape
  u <- (x / scale)^shape

  # HC derivative components
  dH_da <- (1 + u)^(-2) * 2 * u * log(x / scale)  ## derivative w.r.t shape
  dH_db <- (1 + u)^(-2) * 2* u * (-shape / scale) ## derivative w.r.t scale

  # delta-method variance
  var_H <- dH_db^2 * se_scale^2 +
    dH_da^2 * se_shape^2 +
    2 * dH_da * dH_db * cov_shape_scale

  sqrt(var_H)  # return standard error
}

