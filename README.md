# HSPE
Household Survey Parameter Estimation


````markdown
# HSPE

HSPE is an R package for analyzing income distributions, poverty, and record-driven time series methods.

[![CRAN status](https://www.r-pkg.org/badges/version/HSPE)](#)
[![Build Status](https://github.com/jj-64/HSPE/actions/workflows/R-CMD-check.yaml/badge.svg)](#)


## Installation

You can install the development version from GitHub using `remotes` or `devtools`:

```r
install.packages("remotes")
remotes::install_github("jj-64/HSPE")
````
```r
install.packages("devtools")
devtools::install_github("jj-64/HSPE", upgrade = "ask")
````

If you prefer to build locally:

```bash
git clone https://github.com/jj-64/HSPE.git
cd HSPE
R CMD build .
R CMD check HSPE_0.0.0.tar.gz
```

## Quick start

```r
library(HSPE)
# List data included
data()

# Load example dataset
data(SumData)
head(SumData)

# Quick plots
plot_lorenz(SumData[1, ])
plot_country_summary("at94")
```

## Citation

Add citation instructions here.

## Contributing

Please open issues or PRs on GitHub.

````


# 1 - Create summarized data
In the "Execution"" folder, the code "Create Summarized data.R" will read each LIS.dat file from the "DATA" folder.
If "summarize_data()" is run, all files will be summarized. 
The code will save the generated dataframe in a database in "DataProcessed" folder to be used later on easily.
we can call this data using load("DataProcessed/SumData.rda") command.

If only certain files are needed, we can specify an index (as the file number) or a country code.
For more details, go to the function documentation.
