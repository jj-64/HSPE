

## Compute Gini index from raw data without any expansion but using HHWeight ---------
gini_index_weighted <- function(income_vector, weight_vector) {
  ## Combine income and weight into a data frame and sort by income
  df <- data.frame(income = income_vector, weight = weight_vector)
  df <- df[order(df$income), ]

  # Extract sorted incomes and weights
  sorted_incomes <- df$income
  sorted_weights <- df$weight
  n <- length(sorted_incomes)

  # Compute the weighted cumulative sums
  cum_weight <- cumsum(sorted_weights)
  cum_income <- cumsum(sorted_incomes * sorted_weights)
  total_weight <- sum(sorted_weights)
  total_income <- sum(sorted_incomes * sorted_weights)

  # Compute the weighted Gini index
  weighted_sum <- sum((cum_income - (sorted_incomes * sorted_weights)) * sorted_weights)
  gini <- 1 - (2 / (total_weight * total_income)) * weighted_sum

  return(gini)
}

## Lorenz curve from microdata ignoring negative values
Lorenz_curve = function(income_vector){
  income_vector = income_vector[income_vector>=0]
  income_vector <- sort(income_vector,FALSE)
  income <- as.data.frame(income_vector)
  income$FFy = cumsum(income$income_vector)
  income$L = income$FFy/sum(income$income_vector)

  ## Percentile P
  income$rank <- rank(income$income_vector)
  income$P <- income$rank/max(income$rank)

  return(income[,c("P","L")])
}

# Reduce P&L into deciles
Lorenz_curve_reduced <- function(x, P, L) {
  closest_index <- sapply(x, function(p) which.min(abs(P - p)))
  closest_values <- P[closest_index]
  closest_incomes <- L[closest_index]
  return(data.frame(P = x, P_exact = closest_values, L = closest_incomes))
}

## Descriptive Stat of fyy2

Descriptive = function(fyy2){
  ## Make sure to remove negative (& zero) values from the vector before computing the Gini automatically
  fyy2 <- as.numeric(data2[,2])
  fyy2 <- fyy2[fyy2 > 0]

  ## Average
  Average <- mean(fyy2)

  ## Median
  median_value <- median(fyy2)

  ## Standard deviation
  std_dev <- sd(fyy2)

  ## Quantiles
  q10 <- quantile(fyy2, probs = 0.1)
  q20 <- quantile(fyy2, probs = 0.2)
  q25 <- quantile(fyy2, probs = 0.25)
  q30 <- quantile(fyy2, probs = 0.3)
  q40 <- quantile(fyy2, probs = 0.4)
  q75 <- quantile(fyy2, probs = 0.75)

  ## Gini
  Gini1 <- Gini(fyy2, unbiased=FALSE)

  return(list("Mean"= Average, "Median"=median_value, "STD"=std_dev,"Gini"=Gini1,
              "Quantile"=c(q10,q20,q25,q30,q40,q75)))
}

## Outliers
Outlier_range = function(fyy2,r=1.5){
  r1 = mean(fyy2) - r*IQR(fyy2)
  r2 = mean(fyy2) + r*IQR(fyy2)
  return(c(r1,r2))
}


## Income deciles from observed micro data
income_decile_shares <- function(y, probs = seq(0.1, 1, by = 0.1) , names = TRUE) {
  y <- sort(y)
  n <- length(y)

  # Decile cutpoints
  decile_indices <- floor(probs * n)
  decile_indices[decile_indices == 0] <- 1

  # Compute cumulative income shares
  cum_income <- cumsum(y)
  total_income <- sum(y)

  # Income share at each decile
  shares <- cum_income[decile_indices] / total_income * 100

  if(names){
  names(shares) <- paste0("d", 1:length(probs)) }
  return(shares)
}

