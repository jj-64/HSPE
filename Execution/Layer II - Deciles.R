load("DataProcessed/SumData.rda")

library(GB2group)
#library(stringr)

# -------------------------
# MAIN LOOP
# -------------------------
combined_HC = list()
combined_Param = list()
combined_HC_CI = list()

for(i in seq_along(data)){

  # Extract inputs
  Average <- data$Mean[i]
  se_Average <- data$Mean_SE[i]
  Gini1   <- data$Gini[i]
  se_Gini1 <- data$Gini_SE[i]
  N = data$N[i]
  Country = data$Country[i]

  ## Extract observed Headcount
  observed_HC <- get_observed_HC(data, Country)  # threshold + observed_HC
  thresholds <- observed_HC$threshold

  # Lorenz shares
  L_values = as.numeric(get_observed_deciles(data, Country))
  L_nonCum <- income_shares_nonCum(L_values)

  ## Poverty Lines
  PL_vals <-  thresholds * Average

  # Storage
  # H_row <- data.frame(
  #   Dist = dist_names,
  #   PL20 = NA, PL50 = NA, PL80 = NA
  # )

  # PAR_row <- data.frame(matrix(NA, 6, 4)); rownames(PAR_row)=dist_names
  # SE_row  <- data.frame(matrix(NA, 6, 4)); rownames(SE_row)=dist_names
  # KS_row  <- rep(NA, 6)

  # Fit each distribution
    models <- names(CDF_registry)

    SUM_Param <- data.frame(Parameter = unique(unlist(lapply(CDF_registry, `[[`, "params"))))

    results = list()

    ## Loop over all models
    for(model in models) {

      fit <- fit_model_grouped(model, L_nonCum, mean_y = Average,
                               Gini = Gini1, N = N, nrep=1000)

      if(!fit$ok){
        SUM_Param <- compute_param_summary2(SUM_Param, model, NA, NA)
        HC <- NA
               } else {

      # Add parameter estimates to table
      SUM_Param <- compute_param_summary2(SUM_Param, model, fit$par, fit$se)

      # Compute headcounts
      cdf_fun <- get(CDF_registry[[model]]$cdffun)
      HC <- 100 * cdf_fun(PL_vals, as.list(fit$par))

       ## Compute headcounts SE
      HCSE_fun <- get(CDF_registry[[model]]$HCsefun)
      HC_SE <- HCSE_fun(PL_vals, as.list(fit$par), as.list(fit$se))
               }

      results[[model]] <- tibble::tibble(
        Country = Country,
        threshold = thresholds,
        model = model,
        model_H = HC
      )

      SUM_Param$Country = Country

      }

    # combine models output
    model_df <- dplyr::bind_rows(results)

    # merge observed HC
    SUM_H <- model_df %>%
        dplyr::left_join(observed_HC, by = "threshold") %>%
        dplyr::select(Country, threshold, observed_HC, model, model_H)

    ## Save Results
    combined_HC[[i]]  <- SUM_H
    combined_Param[[i]] <- SUM_Param
  }

combined_HC_df  <- bind_rows(combined_HC)
combined_Param_df <- bind_rows(combined_Param)

######### Fit statistics -----------
    # KS statistic (synthetic sample)
    if(model %in% c("DA","SM")){
      sample <- switch(name,
                       Dagum = rdagum(N, par$b, par$a, par$p),
                       SM    = rsm(N, par$b, par$a, par$p))
      emp <- ecdf(sample)(sample)
      theo <- cdf_fun(sample, param)
      KS_row[j] <- max(abs(emp - theo))
    }




  # Record outputs

  SE_row$Country <- df0$r1[i]
  KS_df <- data.frame(Dist=dist_names, KS=KS_row, Country=df0$r1[i])

  results_SE[[i]]  <- SE_row
  results_KS[[i]]  <- KS_df


# -----------------------------
# Combine ALL outputs
# -----------------------------
df_H   <- bind_rows(combined_HC)
df_PAR <- bind_rows(results_PAR)
df_SE  <- bind_rows(results_SE)
df_KS  <- bind_rows(results_KS)

# -----------------------------
# Export
# -----------------------------
write_xlsx(
  list(
    H = df_H,
    Parameters = df_PAR,
    SE = df_SE,
    KS = df_KS
  ),
  "resultsFISKadjusted999.xlsx"
)

