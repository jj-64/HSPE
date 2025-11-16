

## CDF of Fisk -----------
CDF_Fisk = function(y,alpha,beta){
  1-( 1+(y/beta)^alpha )^(-1)
}

## Pdf of Fisk ------------
pdf_Fisk = function(y,alpha,beta){
  alpha*y^(alpha-1)/(beta^alpha*(1+(y/beta)^alpha )^2)
}

## function to compute the mean --------------
mean_Fisk = function(alpha,beta){

  if(alpha>1){
    beta*gamma(1+1/alpha)*gamma(1-1/alpha)
  }
}

## Gini --------------
Gini_Fisk = function(alpha){
  1/alpha
}


#### Fisk mean-matching scale parameter --------------
b_Fisk <- function(alpha, mean_y){
  mean_y / (gamma(1 + 1/alpha) * gamma(1 - 1/alpha))
}

#### Fisk mean-matching scale parameter --------------
a_Fisk <- function(Gini){
  1/Gini
}

#######  Standard Error of a point estimate--------------
se_a_Fisk = function(Gini, se_Gini){
  # Derivative da/dG
  d_a_dG <- -1 / (Gini^2)

  se_a_fisk <- abs(d_a_dG) * se_Gini
  return(se_a_fisk)
}

#######  Standard Error of b point estimate --------------
se_b_Fisk = function(mean_y, Gini, se_mean, se_Gini, a_fisk = NA){

   if(is.na(a_fisk)) {a_fisk = 1/Gini}

  d_a_dG <- -1 / (Gini^2)

  # Derivative db/dmean
  h_a <- gamma(1 + 1/a_fisk) * gamma(1 - 1/a_fisk)
  d_b_dmean <- 1 / h_a

  # Derivative db/da (numerical)
  eps <- 1e-6
  b_plus <- b_Fisk(a_fisk + eps, mean_y)
  b_minus <- b_Fisk(a_fisk - eps, mean_y)
  d_b_da <- (b_plus - b_minus) / (2 * eps)

  # Chain rule db/dG = db/da * da/dG
  d_b_dG <- d_b_da * d_a_dG

  se_b_fisk <- sqrt( (d_b_dmean * se_mean)^2 + (d_b_dG * se_Gini)^2 )
  return(se_b_fisk)
}

## Compute the Lorenz Curve - FISK ---------------
Lorenz_Fisk = function(alpha, P) {

  Lorenz=0
  for(i in 1:length(P)) {
    Lorenz[i] = pbeta(P[i],(1+1/alpha),(1-1/alpha))
  }
  return(Lorenz)
}
