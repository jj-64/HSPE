load("data/SumData.rda")

#library(dplyr)
#library(stringr)

# -------------------------
# MAIN LOOP
# -------------------------
combined_HC = list()
combined_Param = list()
combined_CI = list()

for (i in 1:nrow(data)  ) {

  # Extract inputs
  Average <- data$Mean[i]
  se_Average <- data$Mean_SE[i]
  Gini1   <- data$Gini[i]
  se_Gini1 <- data$Gini_SE[i]
  N = data$N[i]
  Country = data$Country[i]

  ## Extract PL values
  ci_lower_cols <- grep("CI_lower_", names(data), value = TRUE)
  ci_upper_cols <- grep("CI_upper_", names(data), value = TRUE)
  #thresholds <- as.numeric(stringr::str_extract(ci_lower_cols, "\\d+"))/100 ## get thresholds from the suffix

  ## Extract observed Headcount
  #HC_obs_cols <- grep("^HC_", names(data), value = TRUE)
  #HC_obs_cols <- HC_obs_cols[!grepl("SE", HC_obs_cols)]  # remove HC_SE_*
  #HC_obs_cols <- which(colnames(data) %in% paste0("HC_", thresholds*100) )
  #HC_obs = as.numeric(data[i, HC_obs_cols])

  HC_obs = get_observed_HC(data, Country)$observed_HC
  thresholds = get_observed_HC(data, Country)$threshold

  # --- Distribution parameters ---
  fit <- fit_model_limited(dist = c("LN", "NP", "FISK"), mean_y = Average,
                    Gini = Gini1,
                    se_mean = se_Average,
                    se_Gini = se_Gini1)
        ## Fisk
  Fisk_shape  <- fit$FISK$par["shape"]
  Fisk_scale  <- fit$FISK$par["scale"]
        ##LN
  LN_sigma <- fit$LN$par["sigma"]
  LN_mu    <- fit$LN$par["mu"]
        ##NP
  NP_shape <- fit$NP$par["shape"]
  NP_scale <- fit$NP$par["scale"]

  # --- Poverty lines ---
  PL_vals <- list(
    ratio = thresholds,
    pl  = thresholds * Average,
    obs = HC_obs
  )

  # --- Compute headcounts ---
  SUM_H <- compute_headcounts_limited(PL_vals = PL_vals,
                              models = c("all"), # "LN", "FISK", "NP"
                              Average = Average,
                              Fisk_scale = Fisk_scale,
                              Fisk_shape  = Fisk_shape,
                              LN_sigma   = LN_sigma,
                              NP_scale   = NP_scale,
                              NP_shape    = NP_shape)
  SUM_H$Country <- data$Country[i]

  ## Standard Errors for Estimated Parameters
      ## Fisk
  if(!is.na(se_Gini1) & !is.na(se_Average)){
    se_Fisk_shape = fit$FISK$se["se_shape"]
    se_Fisk_scale = fit$FISK$se["se_scale"]
    ## LN
    se_LN_sigma = fit$LN$se["se_sigma"]
    se_LN_mu = fit$LN$se["se_mu"]
        ## NP
    se_NP_shape = fit$NP$se["se_shape"]
    se_NP_scale = fit$NP$se["se_scale"]

  } else{
    se_Fisk_scale =NA; se_Fisk_shape=NA; se_LN_mu=NA; se_LN_sigma=NA; se_NP_shape=NA; se_NP_scale=NA
  }

  # --- Compute parameter summary ---

  SUM_Param <- compute_param_summary(Fisk_scale = Fisk_scale, Fisk_shape = Fisk_shape,
                                         LN_mu = LN_mu, LN_sigma = LN_sigma,
                                         NP_shape = NP_shape, NP_scale = NP_scale,
                                         se_Fisk_shape = se_Fisk_shape, se_Fisk_scale = se_Fisk_scale,
                                         se_LN_mu = se_LN_mu, se_LN_sigma = se_LN_sigma,
                                         se_NP_shape = se_NP_shape, se_NP_scale = se_NP_scale)

  SUM_Param$Country <- data$Country[i]

  # --- Compute Hedacount standard error ----
  if(!is.na(se_Gini1) & !is.na(se_Average)){
  FISK_SE <- sapply(
    PL_vals$pl,
    function(x) HC_se_FISK (
      y = x,
      scale = Fisk_scale,
      shape = Fisk_shape,
      se_scale = se_Fisk_scale,
      se_shape = se_Fisk_shape
    )
  )
  SUM_H$"FISK_H_SE" = as.numeric(FISK_SE)
  SUM_CI_Fisk = conf_bound(SUM_H[,"FISK_H"],FISK_SE)
  colnames(SUM_CI_Fisk) = paste0("FISK_H_",c("lower", "upper"))

  LN_SE <- sapply(
    PL_vals$pl,
    function(x)
      HC_se_LN (x, mean_y = Average,
                sigma = LN_sigma,
                se_mean = se_Average,
                se_sigma = se_LN_sigma,
                cov_mu_sigma = 0)
  )
  SUM_H$"LN_H_SE" = as.numeric(LN_SE)
  SUM_CI_LN = conf_bound(SUM_H[,"LN_H"],LN_SE)
  colnames(SUM_CI_LN) = paste0("LN_H_",c("lower", "upper"))

  NP_SE <- sapply(
    PL_vals$pl,
    function(x) HC_se_NP (
      y = x,
      scale = NP_scale,
      shape = NP_shape,
      se_scale = se_NP_scale,
      se_shape = se_NP_shape
    )
  )
  SUM_H$"NP_H_SE" = as.numeric(NP_SE)
  SUM_CI_NP = conf_bound(SUM_H[,"NP_H"],NP_SE)
  colnames(SUM_CI_NP) = paste0("NP_H_",c("lower", "upper"))

  OBS_CI <- data.frame(
    Country = data$Country[i],
    threshold = thresholds,
    Obs_lower = as.numeric(data[i, ci_lower_cols]),
    Obs_upper = as.numeric(data[i, ci_upper_cols])
  )
  }

  # Store results
  combined_HC[[i]]  <- SUM_H
  combined_Param[[i]] <- SUM_Param
  combined_CI[[i]]  <- dplyr::bind_cols(
    OBS_CI,
    SUM_CI_Fisk,
    SUM_CI_LN,
    SUM_CI_NP
  )

  message(paste("Done Country", data$Country[i],"..."))
}

# --- Combine all results ---
HC_limited_data  <- dplyr::bind_rows(combined_HC)
Param_limited_data <- dplyr::bind_rows(combined_Param)
CI_limited_data <- dplyr::bind_rows(combined_CI)

## --- Output file ---

writexl::write_xlsx(
  list(
    "H"          = HC_limited_data,
    "Parameters" = Param_limited_data,
    "CI" = CI_limited_data
  ),
  path = paste0(here::here("DATA"),"/Limited data.xlsx")
)

save(HC_limited_data, file = "data/HC_limited_data.rda")
save(Param_limited_data, file = "data/Param_limited_data.rda")
save(CI_limited_data, file = "data/CI_limited_data.rda")
