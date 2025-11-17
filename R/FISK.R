
## α > 0 scale
## β >  0  shape

## CDF of Fisk -----------
CDF_Fisk = function(y,scale,shape){
  #1-( 1+(y/shape)^scale )^(-1)
  1/(1+(y/scale)^(-shape))
}

## Pdf of Fisk ------------
pdf_Fisk = function(y,shape,scale){
  shape*y^(shape-1)/(scale^shape*(1+(y/scale)^shape )^2)
}

## function to compute the mean --------------
mean_Fisk = function(shape,scale){

  if(shape>1){
    scale*gamma(1+1/shape)*gamma(1-1/shape)
  }
}

## Gini --------------
Gini_Fisk = function(shape){
  1/shape
}


#### Fisk mean-matching scale (alpha, a) parameter --------------
scale_Fisk <- function(shape, mean_y){
  mean_y / (gamma(1 + 1/shape) * gamma(1 - 1/shape))
}

#### Fisk mean-matching shape (beta, b) parameter --------------
shape_Fisk <- function(Gini){
  1/Gini
}

#######  Standard Error of a point estimate--------------
se_shape_Fisk = function(Gini, se_Gini){
  # Derivative dbeta/dG
  d_b_dG <- -1 / (Gini^2)

  se_shape_Fisk <- abs(d_b_dG) * se_Gini
  return(se_shape_Fisk)
}

#######  Standard Error of a point estimate --------------
se_scale_Fisk = function(mean_y, Gini, se_mean, se_Gini, shape_Fisk = NA){

   if(is.na(shape_Fisk)) {shape_Fisk = 1/Gini}

  d_b_dG <- -1 / (Gini^2)

  # Derivative db/dmean
  h_b <- gamma(1 + 1/shape_Fisk) * gamma(1 - 1/shape_Fisk)
  d_a_dmean <- 1 / h_b

  # Derivative db/da (numerical)
  eps <- 1e-6
  a_plus <- scale_Fisk(shape_Fisk + eps, mean_y)
  a_minus <- scale_Fisk(shape_Fisk - eps, mean_y)
  d_a_db <- (a_plus - a_minus) / (2 * eps)

  # Chain rule db/dG = db/da * da/dG
  d_a_dG <- d_a_db * d_b_dG

  se_shape_Fisk <- sqrt( (d_a_dmean * se_mean)^2 + (d_a_dG * se_Gini)^2 )
  return(se_shape_Fisk)
}

## Compute the Lorenz Curve - FISK ---------------
Lorenz_Fisk = function(shape, P) {

  Lorenz=0
  for(i in 1:length(P)) {
    Lorenz[i] = pbeta(P[i],(1+1/shape),(1-1/shape))
  }
  return(Lorenz)
}

### Headcount SE ------------------
HC_se_Fisk <- function(x, shape, scale, se_shape, se_scale, cov_shape_scale = 0) {

  # compute u = (p/scale)^shape
  u <- (x / scale)^shape

  # HC derivative components
  dH_db <- (1 + u)^(-2) * u * log(x / scale)  ## derivative w.r.t shape
  dH_da <- (1 + u)^(-2) * u * (-shape / scale) ## derivative w.r.t scale

  # delta-method variance
  var_H <- dH_da^2 * se_scale^2 +
    dH_db^2 * se_shape^2 +
    2 * dH_da * dH_db * cov_shape_scale

  sqrt(var_H)  # return standard error
}

