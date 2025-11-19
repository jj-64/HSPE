library(dplyr)
library(readxl)
library(writexl)
library(GB2group)
library(GB2)
library(laeken)
options(scipen = 999)

# Helper: Compute CDF functions -----------------------------
CDF_DA_param <- function(y, param){CDF_DA(y=y, a = param$a, b= param$b, p = param$p) }
CDF_SM_param <- function(y, param){CDF_SM(y = y, a = param$a, b= param$b, q = param$q )}
CDF_B2_param <- function(y, param){ CDF_B2 (y=y, b= param$b, p= param$p, q = param$q) }
CDF_GB2_param <- function(y, param){ CDF_GB2(y =y, a = param$a, b= param$b, p= param$p, q = param$q) }
CDF_NP_param <- function(y, param) {ifelse(y >= param$scale, 1 - 2/(1 + (y/param$scale)^param$shape), 0) }
CDF_FISK_param <- function(y, param){CDF_FISK(y=y, scale=param$b, shape= param$a)}
CDF_LN_param <- function(y, param){CDF_LN(y=y, mean_y= exp(param$mu), s= param$s)}

# Helper: Compute HC SE functions -----------------------------
HC_SE_DA_param <- function(z, param, param_se){HC_se_DA(pov_line=z, b=param$b, a= param$a, p= param$p, se_a= param_se$a, se_b= param_se$b, se_p = param_se$p)}
HC_SE_SM_param <- function(z, param, param_se){HC_se_SM(pov_line=z, b=param$b, a= param$a, q= param$q, se_a= param_se$a, se_b= param_se$b, se_q = param_se$q)}
HC_SE_B2_param <- function(z, param, param_se){HC_se_B2(y=z, b=param$b, p= param$p, q= param$q, se_p= param_se$p, se_b= param_se$b, se_q = param_se$q)}
HC_SE_GB2_param <- function(z, param, param_se){HC_se_GB2(y=z, a=param$a, b=param$b, p= param$p, q= param$q, se_a = param_se$a, se_p= param_se$p, se_b= param_se$b, se_q = param_se$q)}

HC_se_FISK_param <- function(z, param, param_se){HC_se_FISK(pov_line=z, scale=param$b, shape= param$a, se_shape = param_se$a, se_scale = param_se$b)}
HC_se_LN_param <- function(z, param, param_se){HC_se_LN(pov_line=z, mean_y=exp(param$mu), sigma= param$s, se_mean = param_se$a, se_sigma = param_se$b)}
HC_se_NP_param <- function(z, param, param_se){HC_se_NP(pov_line=z, scale=param$scale, shape= param$shape, se_shape = param_se$scale, se_scale = param_se$shape)}

# Helper: CDF Registry -----------------------------
#
# CDF_registry <- list(
#   DA = list(
#     fitfun = "fitgroup.da",
#     fitmicro = "llogis",
#     params  = c("a","b","p"),  ## b is scale
#     cdffun  = "CDF_DA_param",
#     HCsefun = "HC_SE_DA_param"
#   ),
#   SM = list(
#     fitfun = "fitgroup.sm",
#     fitmicro = "",
#     params  = c("a","b","q"),
#     cdffun  = "CDF_SM_param",
#     HCsefun = "HC_SE_SM_param"
#   ),
#   B2 = list(
#     fitfun = "fitgroup.b2",
#     fitmicro = "beta2",
#     params  = c("p","q","b"),
#     cdffun  = "CDF_B2_param",
#     HCsefun = "HC_SE_B2_param"
#   ),
#   GB2 = list(
#     fitfun = "fitgroup.gb2",
#     fitmicro = "GB2",
#     params  = c("a","b","p","q"),
#     cdffun  = "CDF_GB2_param",
#     HCsefun = "HC_SE_GB2_param"
#   ),
#   FISK = list(
#     fitfun = "fitgroup.f",
#     fitmicro = "llogis",
#     params  = c("a","b"),
#     cdffun  = "CDF_FISK_param",
#     HCsefun = "HC_se_FISK_param"
#   ),
#   LN = list(
#     fitfun = "fitgroup.ln",
#     fitmicro = "lognormal",
#     params  = c("s","mu"),
#     cdffun  = "CDF_LN_param",
#     HCsefun = "HC_se_LN_param"
#   ),
#   NP = list(
#     fitfun = "fitgroup.np",
#     fitmicro = "",
#     params = c("shape", "scale"),
#     cdffun = "CDF_NP_param",
#     HCsefun = "HC_se_NP_param"
#   )
# )


CDF_registry <- list(

  # ---------------------------------------------------------
  # 1. LOGNORMAL
  # ---------------------------------------------------------
  LN = list(
    params = c("mu", "s"),

    fitfun = "fitgroup.ln",

    cdffun = function(y, pars) {
      CDF_LN(y, mean_y = exp(pars$mu), sigma = pars$s)
    },

    pdffun = function(y, pars) {
      pdf_LN(y, mean_y = exp(pars$mu), sigma = pars$s)
    },

    lorenzfun = function(p, pars) {
      pnorm(qnorm(p) - pars$s)
    },

    HCsefun = function(z, pars, se) {
      HC_se_LN(z, mean_y = exp(pars$mu), sigma = pars$s, se_mean = se$mu, se_sigma = se$s)
    }
  ),


  # ---------------------------------------------------------
  # 2. FISK (Log-logistic)
  # ---------------------------------------------------------
  FISK = list(
    params = c("a", "b"),

    fitfun = "fitgroup.f",

    cdffun = function(y, pars) {
      CDF_FISK(y, scale = pars$b, shape = pars$a )
    },

    pdffun = function(y, pars) {
      pdf_FISK(y, scale = pars$b, shape = pars$a )
    },

    lorenzfun = function(p, pars) {
      Lorenz_FISK(p, pars$a)
    },

    HCsefun = function(z, pars, se) {
      HC_se_FISK(z, scale = pars$b, shape = pars$a, se_scale = se$b, se_shape = se$a)
    }
  ),

  # ---------------------------------------------------------
  # 3. New Pareto
  # ---------------------------------------------------------
  NP = list(
    params = c("shape", "scale"),

    fitfun = "fitgroup.np",

    cdffun = function(y, pars) {
      CDF_NP(y, scale = pars$scale, shape = pars$shape )
    },

    pdffun = function(y, pars) {
      pdf_NP(y, scale = pars$scale, shape = pars$shape )
    },

    lorenzfun = function(p, pars) {
      Lorenz_NP(p, shape= pars$shape)
    },

    HCsefun = function(z, pars, se) {
      HC_se_NP(z, scale = pars$scale, shape = pars$shape, se_scale = se$scale, se_shape = se$shape)
    }
  ),

  # ---------------------------------------------------------
  # 4. SINGH–MADDALA (Burr XII)
  # ---------------------------------------------------------
  SM = list(
    params = c("a","b","q"),

    fitfun = "fitgroup.sm",

    cdffun = function(y, pars) {
      1 - (1 + (y/pars$b)^pars$a)^(-pars$q)
    },

    pdffun = function(y, pars) {
      a <- pars$a; b <- pars$b; q <- pars$q
      (a*q/b) * (y/b)^(a-1) * (1 + (y/b)^a)^(-q-1)
    },

    lorenzfun = function(p, pars) {
      Lorenz_SM(p, a=pars$a, b = pars$b, q =  pars$q)
    },

    HCsefun = function(z, pars, se) {
      HC_se_SM(y=z, a=pars$a, b = pars$b, q =  pars$q, se_a=se$a, se_b=se$b, se_q=se$q)
    }
  ),


  # ---------------------------------------------------------
  # 5. DAGUM
  # ---------------------------------------------------------
  DA = list(
    params = c("a","b","p"),

    fitfun = "fitgroup.da",

    cdffun = function(y, pars) {
      (1 + (pars$b / y)^pars$a)^(-pars$p)
    },

    pdffun = function(y, pars) {
      a <- pars$a; b <- pars$b; p <- pars$p
      a*p*b^a * y^(-a*p - 1) / (b^a + y^a)^(p+1)
    },

    lorenzfun = function(u, pars) {
      Lorenz_DA(u, a = pars$a, b = pars$b,  p =pars$p )
    },

    HCsefun = function(z, pars, se) {
      HC_se_DA(y=z, a=pars$a, b = pars$b, p =  pars$p, se_a=se$a, se_b=se$b, se_p=se$p)
    }
  ),


  # ---------------------------------------------------------
  # 6. GB2
  # ---------------------------------------------------------
  GB2 = list(
    params = c("a","b","p","q"),

    fitfun = "fitgroup.gb2",

    cdffun = function(y, pars) {
      z <- (y/pars$b)^pars$a
      u <- z / (1 + z)
      pbeta(u, pars$p, pars$q)
    },

    pdffun = function(y, pars) {
      a <- pars$a; b <- pars$b; p <- pars$p; q <- pars$q
      z <- (y/b)^a
      (a * y^(a*p - 1)) / (b^(a*p) * beta(p,q) * (1 + z)^(p+q))
    },

    lorenzfun = function(u, pars) {
      Lorenz_GB2(u, a = pars$a, b = pars$b,  p =pars$p ,q=pars$q)
    },

    HCsefun = function(z, pars, se) {
      HC_se_GB2(y=z, a=pars$a, b = pars$b, p =  pars$p, q=pars$q, se_a=se$a, se_b=se$b, se_p=se$p, se_q = se$q)
    }
  ),


  # ---------------------------------------------------------
  # 7. BETA PRIME (Beta-2)
  # ---------------------------------------------------------
  B2 = list(
    params = c("b","p","q"),

    fitfun = "fitgroup.b2",

    cdffun = function(y, pars) {
      CDF_B2(y, b= pars$b, p=pars$p, q=pars$q)
    },

    pdffun = function(y, pars) {
      p <- pars$p; q <- pars$q; b= pars$b
      (y/b)^(p-1) / (1+y/b)^(p+q) / (b*beta(p,q))
    },

    lorenzfun = function(u, pars) {
      Lorenz_B2(u, b = pars$b,  p =pars$p ,q=pars$q)
    },

    HCsefun = function(z, pars, se) {
      HC_se_B2(y=z, b = pars$b, p =  pars$p, q=pars$q, se_b=se$b, se_p=se$p, se_q = se$q)
    }
  )

)
