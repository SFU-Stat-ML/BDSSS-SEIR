

# Observation density f(yt|xt=[St,Et,Rt,It]')

obs.density <- function(yt, It, lambda, p){
  
  return(dbeta(yt, shape1 = lambda*p*It, shape2 = lambda*(1-p*It)))
  
}

log.obs.density <- function(yt, It, lambda, p){
  
  return(dbeta(yt, shape1 = lambda*p*It, shape2 = lambda*(1-p*It), log=TRUE))
  
}
