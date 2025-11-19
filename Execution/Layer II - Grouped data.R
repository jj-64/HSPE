load("DataProcessed/SumData.rda")

# final storage as clean long-format tibbles
PARAM_ROWS <- list()
HC_ROWS    <- list()
DIAG_ROWS  <- list()

models <- names(CDF_registry)

## Loop over all countries #
for (i in seq_len(nrow(data))) {

  Country    <- data$Country[i]
  Average    <- data$Mean[i]
  se_Average <- data$Mean_SE[i]
  Gini1      <- data$Gini[i]
  se_Gini1   <- data$Gini_SE[i]
  N          <- data$N[i]

  ## observed headcount threshold data
  observed_HC <- get_observed_HC(data, Country)
  thresholds  <- observed_HC$threshold
  PL_vals     <- thresholds * Average

  ## income shares
  L_values  <- get_observed_deciles(data, Country)
  L_nonCum  <- income_shares_nonCum(as.numeric(L_values))

  # ----------------------------
  # Loop over all distribution models
  # ----------------------------
  for (model in models) {

    fit <- fit_model_grouped(
      model,
      L_nonCum,
      mean_y = Average,
      Gini   = Gini1,
      N      = N,
      nrep   = 10
    )

    ## --- failed model
    if (!fit$ok) {

      # (A) parameter rows
      PARAM_ROWS[[length(PARAM_ROWS) + 1]] <- tibble::tibble(
        Country = Country,
        model   = model,
        param   = CDF_registry[[model]]$params,
        value   = NA_real_,
        se      = NA_real_
      )

      # (B) HC rows
      HC_ROWS[[length(HC_ROWS) + 1]] <- tibble::tibble(
        Country    = Country,
        model      = model,
        threshold  = thresholds,
        PL         = PL_vals,
        HC         = NA_real_,
        HC_se      = NA_real_,
        observed_HC = observed_HC$observed_HC
      )

      # (C) diagnostics
      DIAG_ROWS[[length(DIAG_ROWS)+1]] <- tibble::tibble(
        Country = Country,
        model   = model,
        ok      = FALSE
      )
      message("Failed Country ",i, "-" ,Country, " for model ", model)
      next
    }

    ## --- successful model fit ----
    parlist <- as.list(fit$par)
    selist  <- as.list(fit$se)

    # --- A) PARAMETERS long format --------------
    PARAM_ROWS[[length(PARAM_ROWS) + 1]] <- tibble::tibble(
      Country = Country,
      model   = model,
      param   = names(parlist),
      value   = unlist(parlist),
      se      = unlist(selist)
    )

    # --- B) HEADCOUNTS --------------------------
    cdf_fun  <- CDF_registry[[model]]$cdffun #get()
    HC       <-  cdf_fun(PL_vals, parlist)

    HCse_fun <- CDF_registry[[model]]$HCsefun#get()
    HC_se    <- sapply(PL_vals, function(z) HCse_fun(z, parlist, selist))

    HC_ROWS[[length(HC_ROWS) + 1]] <- tibble::tibble(
      Country    = Country,
      model      = model,
      threshold  = thresholds,
      PL         = PL_vals,
      HC         = HC,
      HC_se      = HC_se,
      observed_HC = observed_HC$observed_HC
    )

    # --- C) DIAGNOSTICS -------------------------
    DIAG <-diagnostics_grouped(L_obs = L_values, Lorenz_fun = CDF_registry[[model]]$lorenzfun, parlist, N = N)

    DIAG_ROWS[[length(DIAG_ROWS) + 1]] <- tibble::tibble(
      Country = Country,
      model   = model,
      ok      = TRUE,
      logLik  = DIAG$logLik,
      KS      = DIAG$KS,
      MSE     = DIAG$MSE,
      RMSE    = DIAG$RMSE,
      AIC     = DIAG$AIC,
      BIC     = DIAG$BIC,
    )

    message("Done Country ",i, "-" ,Country, " for model ", model)
  }
}

## bind final outputs -----------
ALL_PARAM <- dplyr::bind_rows(PARAM_ROWS)
HC_ROWS <- lapply(HC_ROWS, function(x) {
  x$HC_se <- as.numeric(x$HC_se)
  x
})
ALL_HC    <- dplyr::bind_rows(HC_ROWS)
ALL_DIAG  <- dplyr::bind_rows(DIAG_ROWS)


## Output file -----

writexl::write_xlsx(
  list(
    "H"          = ALL_HC,
    "Parameters" = ALL_PARAM,
    "Diagonistic" = ALL_DIAG
  ),
  path = paste0(here::here("DataProcessed"),"/Grouped data.xlsx")
)

save(ALL_HC, file = "DataProcessed/HC_Grouped data.rda")
save(ALL_PARAM, file = "DataProcessed/Param_Grouped data.rda")
save(ALL_DIAG, file = "DataProcessed/Diag_Grouped data.rda")
