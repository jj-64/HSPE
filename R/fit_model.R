# Helper: Fit a grouped model -----------------------------
fit_model <- function(model =c("DA","SM","B2","GB2","FISK","LN","NP"), L_nonCum, mean_y, Gini, N) {

  if (!all(model %in% c("DA","SM","B2","GB2","FISK","LN","NP"))) {
    stop("models must be one of:  DA, SM, Beta2, GB2, FISK ,LN","NP")
  }

  info <- CDF_registry[[model]]

  fitfun <- get(info$fitfun)

  tryCatch({
    fit <- fitfun(
      y = L_nonCum,
      pc.inc = mean_y,
      gini.e = Gini,
      gini = TRUE,
      N = N,
      se.ewmd = TRUE,
      se.scale = TRUE,
      nrep = 100
    )

    list(
      ok = TRUE,
      par = fit$ewmd.estimation["Coef.", info$params],
      se  = fit$ewmd.estimation["se",   info$params],
      gini_est = fit$gini.estimation[3]
    )

  }, error = function(e){
    message(model, " FAILED: ", e$message)
    return(list(ok = FALSE))
  })
}

## fit NP----------------

loss_NP <- function(params, pc.inc, gini.e, y, p_dec = seq(0.1,0.9,0.1)){
  a = params[1]
  b = params[2]

  # Rule out invalid parameters
  if (a <= 1 || b <= 0) return(1e10)

  mean_err  = (mean_NP(a,b) - pc.inc)^2
  gini_err  = (Gini_NP(a) - gini.e)^2

  # theoretical deciles
  dec_th = diff(Lorenz_NP(shape=a)) # Quantile_NP(p_dec, a, b)
  dec_err = sum((dec_th - y)^2)

  # weights (tune if needed)
  w1 = 1; w2 = 1; w3 = 10

  w1*mean_err + w2*gini_err + w3*dec_err
}

fit_np <- function(y, pc.inc, gini.e, N, se.ewmd = TRUE, se.scale = TRUE, nrep=10^3){

  p_dec = seq(0.1,0.9,(1/(length(y)+1)))

  out = optim(
    par = c(a=2, b=pc.inc/2),  # starting values
    fn = loss_NP,
    pc.inc = pc.inc,
    gini.e = gini.e,
    y = y,
    p_dec = p_dec,
    method = "L-BFGS-B",
    lower = c(1.01,1e-6),
    upper = c(50,1e6),
    control = list(maxit=500)
  )

  return(out)
}

#fit = fit_npy=y, pc.inc, gini.e)
#fit$par
