load("DataProcessed/SumData.rda")

library(dplyr)
#library(stringr)

# -------------------------
# MAIN LOOP
# -------------------------
combined_HC = list()
combined_Param = list()
combined_CI = list()

for (i in seq_len(nrow(data))) {

  # Extract inputs
  Average <- data$Mean[i]
  se_Average <- data$Mean_SE[i]
  Gini1   <- data$Gini[i]
  se_Gini1 <- data$Gini_SE[i]
  N = data$N[i]
  Country = data$Country[i]

  ## Extract PL values
  #ci_lower_cols <- grep("CI_lower_", names(data), value = TRUE)
  #ci_upper_cols <- grep("CI_upper_", names(data), value = TRUE)
  #thresholds <- as.numeric(stringr::str_extract(ci_lower_cols, "\\d+"))/100 ## get thresholds from the suffix

  ## Extract observed Headcount
  #HC_obs_cols <- grep("^HC_", names(data), value = TRUE)
  #HC_obs_cols <- HC_obs_cols[!grepl("SE", HC_obs_cols)]  # remove HC_SE_*
  #HC_obs_cols <- which(colnames(data) %in% paste0("HC_", thresholds*100) )
  #HC_obs = as.numeric(data[i, HC_obs_cols])

  HC_obs = get_observed_HC(data, Country)$observed_HC
  thresholds = get_observed_HC(data, Country)$threshold

  # --- Distribution parameters ---
        ## Fisk
  Fisk_shape  <- shape_Fisk(Gini1)
  Fisk_scale  <- scale_Fisk(shape = Fisk_shape, mean_y = Average)
        ##LN
  LN_sigma <- Sigma_LN(Gini1)
  LN_mu    <- mu_LN(mean_y = Average, LN_sigma)
        ##NP
  NP_shape <- shape_NP(Gini1)
  NP_scale <- scale_NP(mean_y = Average, shape = NP_shape)

  # --- Poverty lines ---
  PL_vals <- list(
    ratio = thresholds,
    pl  = thresholds * Average,
    obs = HC_obs * 100
  )

  # --- Compute headcounts ---
  SUM_H <- compute_headcounts(PL_vals = PL_vals,
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
    se_Fisk_shape = se_shape_Fisk(Gini1, se_Gini1)
    se_Fisk_scale = se_scale_Fisk(mean_y = Average, Gini = Gini1, se_mean = se_Average, se_Gini = se_Gini1, shape_Fisk = Fisk_shape)
        ## LN
    se_LN_sigma = se_sigma_LN(Average, Gini1,  se_Gini1)
    se_LN_mu = se_mu_LN(Average, Gini1, se_Average, se_Gini1)
        ## NP
    se_NP_shape = se_shape_NP(Gini1, se_Gini1, shape = NP_shape)
    se_NP_scale = se_scale_NP(mean_y = Average, Gini = Gini1, se_mean = se_Average, se_Gini1, shape=NP_shape, se_shape = se_NP_shape, scale= NP_scale)

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
  FISK_SE <- sapply(
    PL_vals$pl,
    function(x) HC_se_Fisk (
      x = x,
      scale = Fisk_scale,
      shape = Fisk_shape,
      se_scale = se_Fisk_scale,
      se_shape = se_Fisk_shape
    ) * 100     # since HC is scaled to percent
  )
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
  SUM_CI_LN = conf_bound(SUM_H[,"LN_H"],LN_SE)
  colnames(SUM_CI_LN) = paste0("LN_H_",c("lower", "upper"))

  NP_SE <- sapply(
    PL_vals$pl,
    function(x) HC_se_Fisk (
      x = x,
      scale = NP_scale,
      shape = NP_shape,
      se_scale = se_NP_scale,
      se_shape = se_NP_shape
    ) * 100     # since HC is scaled to percent
  )
  SUM_CI_NP = conf_bound(SUM_H[,"NP_H"],NP_SE)
  colnames(SUM_CI_NP) = paste0("NP_H_",c("lower", "upper"))

  OBS_CI <- data.frame(
    threshold = thresholds,
    CI_lower = as.numeric(data[i, ci_lower_cols])*100,
    CI_upper = as.numeric(data[i, ci_upper_cols])*100
  )

  # Store results
  combined_HC[[i]]  <- SUM_H
  combined_Param[[i]] <- SUM_Param
  combined_CI[[i]]  <- bind_cols(
    OBS_CI,
    SUM_CI_Fisk,
    SUM_CI_LN,
    SUM_CI_NP
  )

  message(paste("Done Country", data$Country[i],"..."))
}

# --- Combine all results ---
combined_HC_df  <- bind_rows(combined_HC)
combined_Param_df <- bind_rows(combined_Param)
combined_CI_df <- bind_rows(combined_CI)

## --- Output file ---

writexl::write_xlsx(
  list(
    "H"          = combined_HC_df,
    "Parameters" = combined_Param_df,
    "CI" = combined_CI_df
  ),
  path = paste0(here::here("DataProcessed"),"/2Param.xlsx")
)

save(combined_HC_df, file = "DataProcessed/HC_TwoParam.rda")
save(combined_Param_df, file = "DataProcessed/Param_TwoParam.rda")
save(combined_CI_df, file = "DataProcessed/CI_TwoParam.rda")
