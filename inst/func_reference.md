# HSPE Function Reference


## Plotting


- `plot_lorenz(data_row)` — produce Lorenz curve from decile shares.
- `plot_decile_shares(data_row)` — bar plot of decile income shares.
- `plot_headcount_models(HC_row, CI_row)` — compare observed and model headcounts.
- `plot_parametric_fit(param_row)` — show parametric model densities (Fisk, LN, NP).
- `plot_country_summary(country_code)` — combined 4-panel summary.


## Data


- `SumData` — deciles and summary stats.
- `HC_limited_data` — headcount estimates and SEs.
- `CI_limited_data` — confidence intervals for headcounts.
- `Param_limited_data` — fitted parameters and SEs.