# HSPE Plot Functions, Vignette, and Documentation

This document contains **all materials requested**:
- Full `plot_functions.R` file (ready to copy into `R/`)
- Full vignette `Analysis_with_HSPE.Rmd` (ready for `vignettes/`)
- Full **roxygen2 documentation** for all functions

---

# ==============================
# 1. FILE: R/plot_functions.R
# ==============================

```r
#' Lorenz Curve Plot from Decile Shares
#'
#' @description Creates a Lorenz curve using cumulative decile income shares.
#' @param data_row A single-row data frame from SumData.
#' @return A ggplot2 object.
#' @export
plot_lorenz <- function(data_row) {
  library(ggplot2)
  shares <- as.numeric(data_row[paste0("d", 1:10)])
  df <- data.frame(
    p = seq(0.1, 1, by = 0.1),
    L = cumsum(shares) / 100
  )
  ggplot(df, aes(x = p, y = L)) +
    geom_line(size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      title = paste("Lorenz Curve –", data_row$Country),
      x = "Cumulative population share",
      y = "Cumulative income share"
    ) +
    theme_minimal(base_size = 14)
}

#' Decile Income Share Barplot
#'
#' @description Plots the % income received by each decile.
#' @param data_row A single-row data frame from SumData.
#' @return ggplot2 object.
#' @export
plot_decile_shares <- function(data_row) {
  library(ggplot2)
  shares <- as.numeric(data_row[paste0("d",1:10)])
  df <- data.frame(decile = paste0("D", 1:10), share = shares)
  ggplot(df, aes(x = decile, y = share)) +
    geom_col() +
    labs(
      title = paste("Income Distribution by Decile –", data_row$Country),
      y = "% of total income"
    ) + theme_minimal(base_size = 14)
}

#' Observed vs Parametric Headcount
#'
#' @description Compares observed and estimated poverty headcounts.
#' @param HC_row Row from HC_limited_data.
#' @param CI_row Corresponding row(s) from CI_limited_data.
#' @return ggplot2 object.
#' @export
plot_headcount_models <- function(HC_row, CI_row) {
  library(ggplot2)
  df <- data.frame(
    Threshold = CI_row$threshold,
    Observed = HC_row$Observed,
    Fisk = HC_row$FISK_H,
    Lognormal = HC_row$LN_H,
    NP = HC_row$NP_H
  )
  df_long <- tidyr::pivot_longer(df, -Threshold, names_to="Model", values_to="Headcount")

  ggplot(df_long, aes(x = Threshold, y = Headcount, color = Model)) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    scale_x_continuous(labels = scales::percent) +
    labs(
      title = paste("Observed vs Parametric Headcounts –", CI_row$Country[1]),
      y = "Headcount ratio"
    ) + theme_minimal(base_size = 14)
}

#' Parametric Distribution Fit Plot
#'
#' @description Generates density curves using Fisk, Lognormal, and NP parameters.
#' @param param_row Row from Param_limited_data.
#' @return ggplot2 object.
#' @export
plot_parametric_fit <- function(param_row) {
  library(ggplot2)
  x <- seq(0.01, 5 * exp(param_row$LN), length.out = 400)
  df <- data.frame(
    x = x,
    Fisk = VGAM::dfisk(x, shape = param_row$Fisk, scale = param_row$Fisk_SE),
    LN = dlnorm(x, meanlog = param_row$LN, sdlog = param_row$LN_SE),
    NP = VGAM::dpareto(x, shape = param_row$NP, scale = param_row$NP_SE)
  )
  df_long <- tidyr::pivot_longer(df, -x, names_to="Model", values_to="Density")

  ggplot(df_long, aes(x=x, y=Density, color=Model)) +
    geom_line(size = 1.1) +
    labs(
      title = paste("Parametric Fit Comparison –", param_row$Country),
      x = "Income scale (normalized)"
    ) + theme_minimal(base_size = 14)
}

#' Complete Summary Panel for a Country
#'
#' @description Combines Lorenz curve, decile shares, headcount comparison, and parametric fit.
#' @param country_code Country ID (e.g., "at94").
#' @return patchwork layout of 4 ggplots.
#' @export
plot_country_summary <- function(country_code) {
  library(patchwork)
  p1 <- plot_lorenz(SumData[SumData$Country==country_code, ])
  p2 <- plot_decile_shares(SumData[SumData$Country==country_code, ])
  p3 <- plot_headcount_models(
    HC_limited_data[HC_limited_data$Country==country_code, ],
    CI_limited_data[CI_limited_data$Country==country_code, ])
  p4 <- plot_parametric_fit(
    Param_limited_data[Param_limited_data$Country==country_code, ])

  (p1 + p2) / (p3 + p4)
}
```

---

# ===================================
# 2. FILE: vignettes/Analysis_with_HSPE.Rmd
# ===================================

```r
---
title: "Analysis with the HSPE Package"
author: "HSPE Team"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{HSPE Analysis}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Introduction
The **HSPE package** provides a set of harmonized inequality indicators, poverty headcount measures, parametric model fits, and distribution summaries across multiple countries.

This vignette demonstrates:
- Loading HSPE data
- Plotting inequality
- Comparing poverty headcounts
- Visualizing parametric fits
- Generating complete country summaries

# Load the package
```{r}
library(HSPE)
```

# Available datasets
```{r}
ls("package:HSPE")
```

# Example: Inspect summary data
```{r}
head(SumData)
```

# Lorenz curve
```{r}
plot_lorenz(SumData[1, ])
```

# Decile income shares
```{r}
plot_decile_shares(SumData[1, ])
```

# Observed vs parametric headcounts
```{r}
country <- "at94"
plot_headcount_models(
  HC_limited_data[HC_limited_data$Country==country, ],
  CI_limited_data[CI_limited_data$Country==country, ]
)
```

# Parametric density fit
```{r}
plot_parametric_fit(Param_limited_data[Param_limited_data$Country=="at94", ])
```

# Complete Summary
```{r}
plot_country_summary("at94")
```

# Conclusion
This vignette illustrates how HSPE can be used to explore income distribution, inequality, poverty, and parametric modeling.
```

---

# =================================
# 3. Roxygen2 Documentation (all functions)
# =================================

All documentation was inserted directly into the R file.
No additional work required here.

---

If you'd like, I can also generate:
- **pkgdown website** configuration
- **unit tests** (testthat)
- **DESCRIPTION imports cleanup**
- **automated examples** for each dataset
- **README.md** with badges and installation instructions

Just tell me!

