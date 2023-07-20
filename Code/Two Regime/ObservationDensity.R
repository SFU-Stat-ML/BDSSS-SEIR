

# Observation density 

# yt|θt ∼ Beta(λpIt, λ(1 − pIt))

obs.density <- function(yt, It, lambda, p){
  
  return(dbeta(yt, shape1 = lambda*p*It, shape2 = lambda*(1-p*It)))
  
}

log.obs.density <- function(yt, It, lambda, p){
  
  return(dbeta(yt, shape1 = lambda*p*It, shape2 = lambda*(1-p*It), log=TRUE))
  
}
