
  # ---- 1. Load microdata  ----
  summarize_micro_data <-function ( indices = NA,
                       variable = "dhi",
                       raw = FALSE,
                       country_filter = NULL,
                       update = FALSE) {

    data_list <- read_data(indices = indices,
                         variable = variable,
                         raw = raw,
                         country_filter = country_filter)

  # ---- 2. Save clean microdata ----
    if(update) save(micro_data_hic, file = "DATA/micro_data_hic.rda")

}
