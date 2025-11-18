summarize_data <- function(indices = NA,
                           variable = "dhi",
                           raw = FALSE,
                           country_filter = NULL,
                           percentiles = FALSE,
                           pov_lines = c(0.2, 0.5, 0.8)) {

  # ---- 1. Load microdata using your previous function ----
  data_list <- read_data(indices = indices,
                         variable = variable,
                         raw = raw,
                         country_filter = country_filter)

  out <- list()


  for (nm in names(data_list)) {
    df <- data_list[[nm]]

    # Skip empty files
    if (is.null(df) || nrow(df) == 0) next

    y <- df[[variable]]
    y=y[y>0]

    for(models in modelss){
      mlfit.gb2(fyy2)
      fit_FISK <- fitdistrplus::fitdist(as.numeric(fyy2), "llogis", method = "mle", control = list(trace = 1))

      fit_LN <- fitdistr(y, densfun = "lognormal")
      fit_LN <- fitdistrplus::fitdist(y, distr = "lnorm")
      fit_GB2 <- fitdistrplus::fitdist(y, distr = "gb2")
      fit_GB2 <- fitdistrplus::fitdist(y, distr = "dagum")

      #fit_NP <- fitdist(as.numeric(fyy2), "NPareto",start=list(alpha=0.1,Beta=min(fyy2)),
                        #method="mge",control = list(trace = 1))
    }
  }
}
