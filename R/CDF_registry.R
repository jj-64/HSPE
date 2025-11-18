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
CDF_FISK_param <- function(y, param){CDF_FISK(y=y, scale=param$b, shape= param$a)}
CDF_LN_param <- function(y, param){CDF_LN(y=y, mean_y= exp(param$mu), s= param$s)}

# Helper: Compute HC SE functions -----------------------------
HC_SE_DA_param <- function(z, param, param_se){HC_se_DA(pov_line=z, b=param$b, a= param$a, p= param$p, se_a= param_se$a, se_b= param_se$b, se_p = param_se$p)}
HC_SE_SM_param <- function(z, param, param_se){HC_se_SM(pov_line=z, b=param$b, a= param$a, q= param$q, se_a= param_se$a, se_b= param_se$b, se_q = param_se$q)}
HC_se_FISK_param <- function(z, param, param_se){HC_se_FISK(pov_line=z, scale=param$b, shape= param$a, se_shape = param_se$a, se_scale = param_se$b)}
HC_se_LN_param <- function(z, param, param_se){HC_se_LN(pov_line=z, mean_y=exp(param$mu), sigma= param$s, se_mean = param_se$a, se_sigma = param_se$b)}
HC_se_NP_param <- function(z, param, param_se){HC_se_NP(pov_line=z, scale=param$b, shape= param$a, se_shape = param_se$a, se_scale = param_se$b)}

# Helper: CDF Registry -----------------------------

CDF_registry <- list(
  DA = list(
    fitfun = "fitgroup.da",
    fitmicro = "llogis",
    params  = c("a","b","p"),
    cdffun  = "CDF_DA",
    HCsefun = "HC_SE_DA_param"
  ),
  SM = list(
    fitfun = "fitgroup.sm",
    fitmicro = "",
    params  = c("a","b","q"),
    cdffun  = "CDF_SM"
  ),
  B2 = list(
    fitfun = "fitgroup.b2",
    fitmicro = "beta2",
    params  = c("p","q","b"),
    cdffun  = "CDF_B2"
  ),
  GB2 = list(
    fitfun = "fitgroup.gb2",
    fitmicro = "GB2",
    params  = c("a","b","p","q"),
    cdffun  = "CDF_GB2"
  ),
  FISK = list(
    fitfun = "fitgroup.f",
    fitmicro = "llogis",
    params  = c("a","b"),
    cdffun  = "CDF_FISK_param",
    HCsefun = "HC_se_FISK_param"
  ),
  LN = list(
    fitfun = "fitgroup.ln",
    fitmicro = "lognormal",
    params  = c("s","mu"),
    cdffun  = "CDF_LN_param",
    HCsefun = "HC_se_LN_param"
  ),
  NP = list(
    fitfun = "fitgroup.np",
    fitmicro = "",
    params = c("shape", "scale"),
    cdffun = "CDF_NP_param",
    HCsefun = "HC_se_NP_param"
  )
)
