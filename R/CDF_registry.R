library(dplyr)
library(readxl)
library(writexl)
library(GB2group)
library(GB2)
library(laeken)
options(scipen = 999)

# Helper: Compute CDF functions -----------------------------
CDF_DA <- function(y, param){ (1 + (y / param$b)^(-param$a))^(-param$p) }
CDF_SM <- function(y, param){ 1 - (1 + (y / param$b)^param$a)^(-param$Q) }
CDF_B2 <- function(y, param){ pbeta(y/(param$b+y), param$p, param$Q) }
CDF_GB2 <- function(y, param){ pbeta((y/param$b)^param$a/(1+(y/param$b)^param$a), param$p, param$Q) }
CDF_NP_param <- function(y, param) {ifelse(y >= param$scale, 1 - 2/(1 + (y/param$scale)^param$shape), 0) }

# Helper: CDF Registry -----------------------------

CDF_registry <- list(
  DA = list(
    fitfun = "fitgroup.da",
    fitmicro = "llogis",
    params  = c("a","b","p"),
    cdffun  = "CDF_DA"
  ),
  SM = list(
    fitfun = "fitgroup.sm",
    fitmicro = "",
    params  = c("a","b","Q"),
    cdffun  = "CDF_SM"
  ),
  B2 = list(
    fitfun = "fitgroup.b2",
    fitmicro = "beta2",
    params  = c("p","Q","b"),
    cdffun  = "CDF_B2"
  ),
  GB2 = list(
    fitfun = "fitgroup.gb2",
    fitmicro = "GB2",
    params  = c("a","b","p","Q"),
    cdffun  = "CDF_GB2"
  ),
  FISK = list(
    fitfun = "fitgroup.f",
    fitmicro = "llogis",
    params  = c("a","b"),
    cdffun  = "CDF_FISK"
  ),
  LN = list(
    fitfun = "fitgroup.ln",
    fitmicro = "lognormal",
    params  = c("a","b"),
    cdffun  = "CDF_LN"
  ),
  NP = list(
    fitfun = "fitgroup.np",
    fitmicro = "",
    params = c("shape", "scale"),
    cdffun = "CDF_NP_param"
  ),

  Fisk2 <- list(
    fitfun = "fit_Fisk",
    params = c("shape", "scale")
  )

)
