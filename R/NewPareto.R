#install.packages("hypergeo")
library(hypergeo)
library(nloptr)
# Estimate alpha via optimization ---------------
alpha_NP = function(Gini){
  ## Objective function
  eval_f0 <- function( alpha,a,b){
    return(abs(Gini_NP(alpha[1]) - Gini))
  }

  # constraint function g(x) <=0
  eval_g0 <- function(alpha, a,b) {
    return(b- a*alpha[1] )
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

# Wrapper for alpha as function of G (so derivative uses numeric differentiation)
alpha_func_wrapper <- function(Gini) {
  # call your alpha_NP implementation; return scalar numeric
  # add tryCatch because alpha_NP uses optimizer and may fail
  res <- tryCatch({
    as.numeric(alpha_NP(Gini))
  }, error = function(e) NA_real_)
  res
}

# Numerical derivative for α(G) ------------------
alpha_NP_deriv <- function(Gini){
  eps <- 1e-6
  (alpha_NP(Gini + eps) - alpha_NP(Gini - eps)) / (2 * eps)
}

#Estimating beta from alpha and mean by Optimization -------------
beta_NP = function(mean_y , alpha){
  ## Objective function
  eval_f0 <- function( beta,a,b){
    return( abs((2*alpha *beta *Re(hypergeo(2,2-((1+alpha)/alpha),3-((1+alpha)/alpha),-1))/(alpha - 1)) - mean_y ) )
  }

  # constraint function g(x) <=0
  eval_g0 <- function(beta, a,b) {
    return(b- a*beta[1] )
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

# Wrapper for beta as function of (mean_y, alpha)
beta_func_wrapper <- function(x) {
  # x is numeric vector: c(mean_y, alpha)
  mean_y <- x[1]; alpha <- x[2]
  res <- tryCatch({
    as.numeric(beta_NP(mean_y, alpha))
  }, error = function(e) NA_real_)
  res
}

## Standard error for alpha --------------------------
# Delta-method function
# mean_y: sample mean
# Gini: sample Gini
# se_mean: SE(sample mean)
# se_Gini: SE(sample Gini)
# cov_mean_gini: optional covariance between mean and gini (default 0)
# eps: finite-diff step
se_alpha_NP <- function(mean_y, Gini, se_mean, se_Gini, alpha =NA,
                        cov_mean_gini = 0, eps = 1e-6) {
  # 1) alpha and its derivative w.r.t Gini
  if(is.na(alpha)) alpha <- alpha_func_wrapper(Gini)
  if (is.na(alpha)) stop("alpha_NP returned NA or failed at provided Gini")

  d_alpha_dG <- num_deriv(function(G) alpha_func_wrapper(G), Gini, eps = eps)

  var_alpha <- (d_alpha_dG^2) * (se_Gini^2)
  se_alpha  <- sqrt(var_alpha)

  return(se_alpha)

}

## Standard error for beta --------------------------
# Delta-method function
# mean_y: sample mean
# Gini: sample Gini
# se_mean: SE(sample mean)
# se_Gini: SE(sample Gini)
# cov_mean_gini: optional covariance between mean and gini (default 0)
# eps: finite-diff step
se_beta_NP <- function(mean_y, Gini, se_mean, se_Gini, alpha=NA, beta= NA,
                       cov_mean_gini = 0, eps = 1e-6) {
  # 1) alpha and its derivative w.r.t Gini
  if(is.na(alpha)) alpha <- alpha_func_wrapper(Gini)
  if (is.na(alpha)) stop("alpha_NP returned NA or failed at provided Gini")

  d_alpha_dG <- num_deriv(function(G) alpha_func_wrapper(G), Gini, eps = eps)

  var_alpha <- (d_alpha_dG^2) * (se_Gini^2)
  se_alpha  <- sqrt(var_alpha)

  # 2) beta and partial derivatives wrt mean_y and alpha
  if(is.na(beta)) beta <- beta_func_wrapper(c(mean_y, alpha))
  if (is.na(beta)) stop("beta_NP returned NA or failed at provided mean_y & alpha")

  # partial derivative wrt mean_y
  d_beta_dMean  <- num_deriv(function(x) beta_func_wrapper(c(x, alpha)), mean_y, eps = eps)
  # partial derivative wrt alpha (hold mean_y fixed)
  d_beta_dAlpha <- num_deriv(function(a) beta_func_wrapper(c(mean_y, a)), alpha, eps = eps)

  # variance/covariance matrix of (mean_y, alpha)
  var_mean  <- se_mean^2
  var_alpha <- var_alpha
  cov_ma    <- cov_mean_gini * d_alpha_dG  # if you provided Cov(mean, Gini), transform to cov(mean,alpha)
  # Explanation: cov(mean, alpha) ≈ (d alpha / dG) * cov(mean, G)
  # If user provided cov_mean_gini already as Cov(mean,alpha), then they should pass that value.
  # We'll accept both: if cov_mean_gini is small relative to var_mean & var_alpha, it's effectively 0.

  # If user probably passed Cov(mean,Gini), cov_ma computed above; otherwise assume they passed cov(mean,alpha)
  # We'll check a heuristic: if cov_mean_gini magnitude is > 1e-8 and seems like Cov(mean,Gini) (rare),
  # still use transformation. If user passed cov(mean,alpha) directly, they can pass cov_mean_gini already transformed.

  # Build gradient and covariance
  grad <- c(d_beta_dMean, d_beta_dAlpha)
  Sigma <- matrix(c(var_mean, cov_ma, cov_ma, var_alpha), nrow = 2)

  var_beta <- as.numeric(t(grad) %*% Sigma %*% grad)
  se_beta  <- sqrt(max(0, var_beta))

  return(se_beta)
}
## Function for the Lorenz Curve for the country -------------
Lorenz_NP = function(y, alpha, beta) { ## the income group y, alpha and beta of the NPareto, mean income

  ## vector for the intergrals

  integrand <- function(z) {z * pdf_NP(alpha, beta, z) }
  int= c()
  for(i in seq_along(y)) int[i] = integrate(integrand, lower = beta, upper = y[i])$value

  ## Formula
  Lorenz = 1/mean(y) * int
  return(Lorenz)
}

## Quantile function to build Lorenz curve, returns the percentile ------------
Percentile_NP = function(alpha, q) {

  doubleInt <- function(q,alpha){
    f <- function(t,alpha){((1+t)/(1-t))^(1/alpha)}
    fvalue = c()
    for (i in seq_along(q)) fvalue[i] <- integrate(f,lower=0,upper=q[i],alpha=alpha)$value
    return(fvalue)
  }

  return ( doubleInt(q,alpha)/  Re(alpha * hypergeo::hypergeo(1,-1/alpha,2-1/alpha,-1) / (alpha-1) ) )

}

## pdf of NewPareto --------------------
pdf_NP = function(alpha,beta,p) {
  return(ifelse(p>= beta, (2 * alpha * beta^alpha *p^(alpha-1) ) / (p^alpha + beta^alpha)^2 ,0))
}

## Function 10:  CDF of NewPareto ----------------
CDF_NP = function(alpha,beta,q) {
  return(ifelse(q>= beta, (q^alpha - beta^alpha) / (q^alpha + beta^alpha),0))
}

## Quantile Function from percentiles ------------------
Quantile_NP = function(alpha,beta,p) {
  return(ifelse(p<1,beta* ((-p-1)/(p-1))^(1/alpha),100000000))
}

## Radom Generated NP ----------------------
#NP distributed x: Given a uniform variate U drawn from U(0, 1) distribution, we obtain X, which is NP distributed and is given by
random_NP <- function (u,a,b){(((2*(b^a))/u)-(b^a))^(1/a)}

## Compute the Gini from NP distribution ------------------
Gini_NP <- function(alpha){

  doubleInt <- function(p,alpha){
    f <- function(t,alpha){((1+t)/(1-t))^(1/alpha)}
    fvalue <- integrate(f,lower=0,upper=p,alpha=alpha)$value
    return(fvalue)
  }

  integrationTrap <- function (func, lower, upper, alpha){
    x = seq(lower,upper,by=0.001)  ## divide the x-axis into small intervals
    n = length(x)
    x <- as.matrix(x,n,1)
    Y <- matrix(0,n,1)
    Y <- as.matrix(Y,n,1)
    tra <- matrix(0,n,1)

    for(i in 1:n)
    {
      Y[i] <- func(x[i],alpha)  ## the integral over evey portion/interval
    }

    tra[1] <- 0
    for (i in 2:n) {
      tra[i] <-((0.5)*(x[i]-x[i-1])*(Y[i]+Y[i-1]))  }
    return(sum(tra))
  }

  ## Gini integral
  integ = integrationTrap(func=doubleInt,lower=0, upper = 1, alpha=alpha)/Re(alpha * hypergeo(1,-1/alpha,2-1/alpha,-1) / (alpha-1) )
  gini <- 1- 2*integ

  ## Return
  return(gini)
}

##MLE estimate -----------------------------------
#Code from paper for maximum likelihood, returns alpha and beta estimates
functionMLE_NP <- function(z){                 ## Code from paper for maximum likelihood, returns alpha and beta estimates
  z <- sort(z)
  x <- z[2:length(z)]
  b <- min(z)
  alphaguess <- length(z)/sum(log(z)-log(b))

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

  mm  <- alphaguess
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
  alphaResult <- res
  betaResult <- b
  results <- list(alphaEstimate = alphaResult,
                  betaEstimate = betaResult)
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
#   while (runif(1) > Fy_NPareto(estimated_params$alphaEstimate, estimated_params$betaEstimate, y)) {
#     u <- runif(1)
#     y <- min_value + u * (max_value - min_value)
#   }
#   new_values[i] <- y
# }


### Method 2 of estimating micro data ####
#
# new_values1 <- random_NP(runif(fyy2), estimated_params$alphaEstimate, estimated_params$betaEstimate)
#
# sorted_new_values <- sort(new_values)
# sorted_new_values1 <- sort(new_values1)
#
# df <- data.frame(DataVector = sorted_new_values, NewValues = sorted_new_values1)


# Requires: nloptr (for your existing alpha_NP / beta_NP), stats

# Numeric derivative helper (centered difference)
num_deriv <- function(f, x, idx = 1, eps = 1e-6, ...) {
  # f: function taking numeric vector x (or scalar), returns numeric scalar
  # idx: index of x to perturb (for scalar x idx=1)
  x <- as.numeric(x)
  e <- eps
  x_plus  <- x
  x_minus <- x
  x_plus[idx]  <- x_plus[idx]  + e
  x_minus[idx] <- x_minus[idx] - e
  (f(x_plus, ...) - f(x_minus, ...)) / (2*e)
}


# ---------------------------
# Bootstrap function (recommended if you have raw income vector)
# ---------------------------
# incomes: numeric vector of person-level incomes (or unit-level)
# R: number of bootstrap resamples (e.g., 1000)
# parallel: "no" or use multicore? (not implemented here)
se_newpareto_boot <- function(incomes, R = 1000, seed = 1234, alpha_init = NULL) {
  set.seed(seed)
  n <- length(incomes)

  # storage
  alpha_bs <- numeric(R)
  beta_bs  <- numeric(R)
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

    # try alpha and beta; wrap in tryCatch
    al <- tryCatch(alpha_NP(g_hat), error = function(e) NA_real_)
    be <- NA_real_
    if (!is.na(al)) {
      be <- tryCatch(beta_NP(m_hat, al), error = function(e) NA_real_)
    }

    mean_bs[r] <- m_hat
    gini_bs[r] <- g_hat
    alpha_bs[r] <- al
    beta_bs[r]  <- be

    setTxtProgressBar(pb, r)
  }
  close(pb)

  # Remove failed iterations
  ok <- is.finite(alpha_bs) & is.finite(beta_bs)
  if (sum(ok) < max(10, 0.1 * R)) {
    warning("Many bootstrap iterations failed. Check alpha_NP / beta_NP stability.")
  }

  alpha_bs <- alpha_bs[ok]
  beta_bs  <- beta_bs[ok]
  mean_bs  <- mean_bs[ok]
  gini_bs  <- gini_bs[ok]

  # Estimates
  alpha <- mean(alpha_bs)
  beta_hat  <- mean(beta_bs)

  se_alpha <- sd(alpha_bs)
  se_beta  <- sd(beta_bs)

  cov_mean_gini <- cov(mean_bs, gini_bs)
  cov_mean_alpha <- cov(mean_bs, alpha_bs)
  cov_alpha_beta <- cov(alpha_bs, beta_bs)

  list(
    alpha = alpha,
    se_alpha = se_alpha,
    beta_hat = beta_hat,
    se_beta = se_beta,
    cov_mean_gini = cov_mean_gini,
    cov_mean_alpha = cov_mean_alpha,
    cov_alpha_beta = cov_alpha_beta,
    raw = list(mean_bs = mean_bs, gini_bs = gini_bs, alpha_bs = alpha_bs, beta_bs = beta_bs)
  )
}

