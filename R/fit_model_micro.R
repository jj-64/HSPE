
## Create a wrapper that safely fits each model---------------
## incomplete
safe_fit <- function(y, distr, start = NULL, lower = -Inf, upper = Inf){

  out <- try(
    fitdistrplus::fitdist(y, distr = distr, method = "mle",
                          start = start, lower = lower, upper = upper),
    silent = TRUE
)
  if(inherits(out, "try-error")) return(NULL)
  out
}

## Function that fits ALL models to one dataset--------------
fit_models_micro <- function(y){

  Average <- mean(y)

  y <- y[y>0]

  fits <- list(
    LN   = safe_fit(y, distr = "lnorm")
    ,
    FISK = safe_fit(y, distr = "llogis",
                    start = list(shape1.a = 1, scale = Average))
                    ,   # VGAM name for Fisk
    B2   = safe_fit(y, "b2",
                    start = list(shape1 = 1, shape2=1, scale = Average),
                    lower = c(0.01, 0.01, 0.01))
                    ,
    NP   = safe_fit(y, distr = "newpareto",
                    start = list(shape = 0.1, scale = Average),
                    lower = c(0.01, min(y))  )
    ,
    GB2  = safe_fit(y, "gb2",
                    start = list(shape1 = 2, shape2 = 2, shape3 = 1.5, scale = Average),
                    lower = c(0.01, 0.01, 0.01, 0.01))
    ,
    DA = safe_fit(y, distr= "dagum",
                     start = list(shape1.a = 1.5, shape2.p = 2.0, scale = Average),
                     lower = c(0.01, 0.01, 0.01))
    ,
    SM = safe_fit(y, "sinmad",     # Singh-Maddala in VGAM
                  start = list(shape1.a = 1, shape3.q = 1, scale = Average))
  )
  fits
}

## Extract parameter info from the fit_model_micro -------------
extract_info <- function(fit, model, file){
  if(is.null(fit)) {
    return(tibble(
      file, model,
      status = "FAILED",
      logLik = NA, #AIC = NA, BIC = NA,
      parameter = NA,
      estimate = NA,
      sd = NA
    ))
  }

  params <- names(fit$estimate)
  tibble(
    file,
    model,
    status = "OK",
    logLik = as.numeric(logLik(fit)),
    #AIC = AIC(fit),
    #BIC = BIC(fit),
    parameter = params,
    estimate = fit$estimate,
    sd = sqrt(diag(fit$vcov))
  )
}
