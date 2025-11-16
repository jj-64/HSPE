# -------------------------
# Helper function: poverty headcounts
# -------------------------
compute_headcounts <- function(PL_vals, Average, a, b, sigma, alpha, beta) {
  data.frame(
    Observed = PL_vals$obs,
    FISK_H   = sapply(PL_vals$pl,  function(p) CDF_Fisk(y=p, alpha = a, beta = b) * 100),
    LN_H     = sapply(PL_vals$pl,  function(p) pnorm((log(p / Average)/sigma) + sigma/2) * 100),
    NP_H     = sapply(PL_vals$pl,  function(p) CDF_NP(alpha, beta, q=p) * 100),
  )
}

# -------------------------
# Helper function: parameters summary
# -------------------------
compute_param_summary <- function(a_FISK, b_FISK, mu_LN, sigma_LN, alpha_NP, beta_NP,
                                  se_a = NA, se_b = NA, se_mu = NA, se_sigma = NA, se_alpha = NA, se_beta = NA) {
  data.frame(
    row = c("Shape/μ", "Scale/σ"),
    Observed_FISK = c(a_FISK, b_FISK),
    SE_FISK = c(se_a, se_b),
    Observed_LN = c(mu_LN, sigma_LN),
    SE_LN = c(se_mu, se_sigma),
    Observed_NP = c(alpha_NP, beta_NP),
    SE_NP = c(se_alpha, se_beta)
  )
}

# -------------------------
# Helper function: Headcount Ratio Variance
# -------------------------
headcounts_SE = function(H, N){
  if (H>1 | H <0) stop("H should be a proportion between 0 and 1")
  sqrt(H*(1-H)/N)
}

compute_headcounts_SE = function(H,N){
  data.frame(
FISK_SE = sapply(H, FUN = function(p) headcounts_SE(H=p, N=N)),
LN_SE   = NA,
NP_SE   = NA
)
}

library(dplyr)
##Function to expand data according to HH weight vector -----------
replicate_rows <- function(df) {
  df = df %>%
    slice(rep(1:n(), hwgt))  %>%
    ungroup()
  return(as.data.frame(df))
}

##Function to expand data according to HH weight vector -------------
replicate_rows_p <- function(df) {
  df = df %>%
    slice(rep(1:n(), pwgt))  %>%
    ungroup()
  return(as.data.frame(df))
}
