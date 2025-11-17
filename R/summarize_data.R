#' Summarize LIS Microdata Files
#'
#' Reads LIS microdata (via read_LIS_data()) and computes summary statistics for each file:
#' mean, standard error, Gini index, Gini SE, distribution quantiles, and poverty headcounts.
#'
#' @param indices File indices passed to read_LIS_data().
#' @param variable Income variable name (default "dhi").
#' @param raw Whether to expand by household weights.
#' @param country_filter Optional filter to select specific files.
#' @param percentiles Logical.
#'        If TRUE: compute p1–p100.
#'        If FALSE: compute d1–d10 (deciles).
#' @param pov_lines Poverty lines expressed as multiples of the mean.
#'        Default = c(0.2, 0.5, 0.8).
#'
#' @return A data.frame of summary results, one row per dataset.
#' @example
#' sum_data = summarize_data(country_filter = "at97")
#' View(sum_data)
#' @export
#'
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
    N <- length(y)

    # ---- 2. Mean & SE ----
    avg <- mean(y)
    se_avg <- sd(y) / sqrt(N)

    # ---- 3. Gini index & SE ----
    Gini <- ineq::ineq(y, type = "Gini")
    Gini_se <- sqrt(giniVarCI::igini(y, interval = "zjackknife")$Variance)

    # ---- 4. Quantiles ----
    if (percentiles) {
      qs <- income_decile_shares(y,probs = seq(0.01, 1, 0.01), names = TRUE) #quantile(y, probs = seq(0.01, 1, 0.01), names = TRUE)
      #quant_names <- names(qs)
    } else {
      qs <- income_decile_shares(y,probs = seq(0.01, 1, 0.1), names = TRUE)
      #quant_names <- paste0("d", 1:10)
      #names(qs) <- quant_names
    }

    # ---- 5. Poverty lines & headcounts ----
    PL_values <- avg * pov_lines
    names(PL_values) <- paste0("PL_", pov_lines * 100)

    HC <- sapply(PL_values, function(pl) mean(y < pl))
    SE_HC <- sqrt(HC * (1 - HC) / N)

    ## Confidence iterval constrained between 0 and 1
    CI_HC <- conf_bound(HC, SE_HC, significance = 0.95)
    CI_HC$lower <- pmax(CI_HC$lower, 0)
    CI_HC$upper <- pmin(CI_HC$upper, 1)

    names(HC) <- paste0("HC_", pov_lines * 100)
    names(SE_HC) <- paste0("HC_SE_", pov_lines * 100)
    CI_vec <- c(CI_HC$lower, CI_HC$upper)
    names(CI_vec) <- c(
      paste0("CI_lower_", pov_lines * 100),
      paste0("CI_upper_", pov_lines * 100)
    )

    # ---- 6. Country name extracted from file name ----
    country <- sub("\\.dta$", "", nm)

    # ---- 7. Build a clean row ----
    out[[nm]] <- c(
      Country = country,
      N = N,
      Mean = avg,
      Mean_SE = se_avg,
      Gini = Gini,
      Gini_SE = Gini_se,
      qs,
      HC,
      SE_HC,
      CI_vec
    )
  }

  # ---- 8. Return as a dataframe ----
  result <- do.call(rbind, out)
  result <- as.data.frame(result)

  # Convert numeric columns correctly
  numeric_cols <- setdiff(names(result), "Country")
  result[numeric_cols] <- lapply(result[numeric_cols], as.numeric)

  return(result)
}
