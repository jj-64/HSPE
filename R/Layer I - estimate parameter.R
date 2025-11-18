estimate_param <- function(dist = c("LN", "NP", "FISK"), mean_y, Gini, se_mean = NA, se_Gini = NA){
  dist <- match.arg(dist)

  if(dist == "FISK"){
    fit_Fisk_param (mean_y, Gini, se_mean = NA, se_Gini = NA)
  }

  if(dist == "LN"){
    fit_LN_param(mean_y, Gini, se_mean = NA, se_Gini = NA)
  }
}


