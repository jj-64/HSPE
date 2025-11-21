# library(fitdistrplus)
# library(VGAM)        # for Fisk, Dagum, Singh-Maddala, GB2
# library(tidyverse)

load("DATA/DataProcessed/micro_data_hic.rda")

dllogis = VGAM::dfisk; pllogis = VGAM::pfisk ; qllogis = VGAM::qfisk
dnewpareto = pdf_NP; pnewpareto = CDF_NP; qnewpareto = Quantile_NP
ddagum <- function(x, shape1.a, shape2.p, scale, log = FALSE) {
  if (length(x) == 0) return(numeric(0))
  VGAM::ddagum(x, shape1.a = shape1.a, shape2.p = shape2.p, scale=scale, log)
}

pdagum <- function(q, shape1.a, shape2.p, scale) {
  if ( length(q) == 0) return(numeric(0))
  VGAM::pdagum(q, shape1.a = shape1.a, shape2.p = shape2.p, scale=scale, lower.tail = TRUE, log.p = FALSE)
}

qdagum = function(p, shape1.a, shape2.p, scale, lower.tail = TRUE, log.p = FALSE) {
  if (length(p) == 0) return(numeric(0))
  VGAM::qdagum(p, shape1.a = shape1.a, shape2.p = shape2.p, scale=scale, lower.tail, log.p)
}

db2 = extraDistr::dbetapr; pb2 = extraDistr::pbetapr; qb2 = extraDistr::qbetapr

ddagum <- function (x, scale = 1, shape1.a, shape2.p)
    { if (length(x) ==0) return(0)
      LLL <- max(length(x), length(shape1.a), length(scale), length(shape2.p))
      if (length(x) != LLL)
        x <- rep_len(x, LLL)
      if (length(shape1.a) != LLL)
        shape1.a <- rep_len(shape1.a, LLL)
      if (length(scale) != LLL)
        scale <- rep_len(scale, LLL)
      if (length(shape2.p) != LLL)
        shape2.p <- rep_len(shape2.p, LLL)
      Loglik <- rep_len(log(0), LLL)
      xok <- (x > 0) & !is.na(x)
      Loglik[xok] <- log(shape1.a[xok]) + log(shape2.p[xok]) +
        (shape1.a[xok] * shape2.p[xok] - 1) * log(x[xok]) - shape1.a[xok] *
        shape2.p[xok] * log(scale[xok]) - (1 + shape2.p[xok]) *
        log1p((x[xok]/scale[xok])^shape1.a[xok])
      Loglik[shape2.p <= 0] <- NaN
      x.eq.0 <- (x == 0) & !is.na(x)
      Loglik[x.eq.0] <- log(shape1.a[x.eq.0]) + log(shape2.p[x.eq.0]) -
        shape1.a[x.eq.0] * shape2.p[x.eq.0] * log(scale[x.eq.0])
      Loglik[is.na(x)] <- NA
      Loglik[is.nan(x)] <- NaN
      Loglik[x == Inf] <- log(0)
      return(exp(Loglik))
}

## Loop over all dat files

process_all_files <- function(path = "."){

  #files <- list.files(path, pattern = "\\.dat$", full.names = TRUE)
  files <- names(micro_data_hic)

  results <- purrr::map_dfr(files, function(f){

    #y <- scan(f, quiet = TRUE)
    y <- micro_data_hic[[i]]

    fits <- fit_model_micro(y)

    purrr::map2_dfr(fits, names(fits),
             ~ extract_info(.x, model = .y, file = basename(f)))
  })

  results
}

## Run
all_results <- process_all_files("path/to/dat/files")

write_csv(all_results, "fit_results.csv")


