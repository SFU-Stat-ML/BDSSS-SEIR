
# Takes a vector of log weights, and return the vector of normalized exponentiated values.

normalize_weight <- function(logweights){
  mlw <- max(logweights)
  w <- exp(logweights - mlw)
  nw <- w / sum(w) # normalized weights
  return(nw)
}