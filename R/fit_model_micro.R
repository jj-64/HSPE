#library("fitdistrplus")

fitdagum <-function(y){
ddagum  <-function(x,a, b, c) VGAM::ddagum(x, scale= b, shape1.a = a, shape2.p = c)
pdagum<-  function(q,a, b, c) VGAM::pdagum(q, scale= b, shape1.a = a, shape2.p = c)
qdagum <-  function(p,a, b, c) VGAM::qdagum(p, scale= b, shape1.a = a, shape2.p = c)

 return(fitdistrplus::fitdist(y, "dagum", start=list(a=2, b=1000, c=2)))
}

# Helper: Fit a grouped model -----------------------------
fit_model_micro <- function(model =c("DA","SM","B2","GB2","FISK","LN","NP"), L_nonCum, mean_y, Gini, N) {

  if (!all(model %in% c("DA","SM","B2","GB2","FISK","LN","NP"))) {
    stop("models must be one of:  DA, SM, Beta2, GB2, FISK ,LN","NP")
  }

  info <- CDF_registry[[model]]

  fitfun <- get(info$fitmicro)

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
