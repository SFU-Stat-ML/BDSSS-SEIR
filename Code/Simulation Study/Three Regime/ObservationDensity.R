

# Observation density 

obs.density <- function(yt, It, lambda, p){
  
  return(dbeta(yt, shape1 = lambda*p*It, shape2 = lambda*(1-p*It)))
  
}

log.obs.density <- function(yt, It, lambda, p){
  
  return(dbeta(yt, shape1 = lambda*p*It, shape2 = lambda*(1-p*It), log=TRUE))
  
}
