
# ===================================
# 2. FILE: vignettes/Analysis_with_HSPE.Rmd
# ===================================

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

