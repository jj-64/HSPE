
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

# Available functions
```{r}
ls("package:HSPE")
```

# Example: Inspect summary data
```{r}
head(SumData)
```

# Example: Inspect Computations from limited data
```{r}
head(CI_limited_data)
head(HC_limited_data)
head(Param_limited_data)
```

# Example: Inspect Computations from grouped data
```{r}
head(CI_grouped_data)
head(HC_grouped_data)
head(Param_grouped_data)
```

# Plot: Lorenz curve
```{r}
plot_lorenz(SumData[1, ])
```

# Plot: Decile income shares
```{r}
plot_decile_shares(SumData[1, ])
```

# Plot: Observed vs parametric headcounts
```{r}
country <- "at94"
plot_headcount_models(
  HC_limited_data[HC_limited_data$Country==country, ],
  CI_limited_data[CI_limited_data$Country==country, ]
)
```

# Plot: Parametric density fit
```{r}
plot_parametric_fit(Param_limited_data[Param_limited_data$Country=="at94", ])
```

# Plot: Complete Summary
```{r}
plot_country_summary("at94")
```

# Conclusion
This vignette illustrates how HSPE can be used to explore income distribution, inequality, poverty, and parametric modeling.
```



