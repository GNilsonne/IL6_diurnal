## Function to calculate the acrophase of from a cosinor model where sin & cos are 
## the estimated coefficients from covariates describing the 24h period of the 
## sine & cosine functions

acroPhase <- function(sin, cos) {
  if (is.na(sin) | is.na(cos)) { return(NA) }
  if (cos >= 0 & sin  > 0) {K=0 }
  if (cos <  0 & sin >= 0) {K=12} 
  if (cos <= 0 & sin  < 0) {K=12} 
  if (cos >  0 & sin <= 0) {K=24} 
  return(atan(sin/cos)/(2*pi)*24 + K)
}


## Demonstrating acrophase calculations
for (p in seq(0,24, .25)){
  x=seq(-0, 24,.1)
  y = cos((x-p)/24*pi*2) + rnorm(length(x))*.001 ## generate data with phase=p + small error

  sin = sin(x/24*pi*2)
  cos = cos(x/24*pi*2)
  
  m=lm(y~sin+cos)

  ap=acroPhase(coef(m)[["sin"]], coef(m)[["cos"]])
  
  print(paste0("True phase / Estimated phase:  ",  format(p, nsmall=6), " / ", format(ap, nsmall=6)))
}