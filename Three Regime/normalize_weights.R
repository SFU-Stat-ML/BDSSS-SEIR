#############################################################
# Author: Jingxue (Grace) Feng
#         Simon Fraser University, Burnaby, BC, Canada
#         Email: jingxuef@sfu.ca
#         Date: January 17, 2022
#############################################################

# Takes a vector of log weights, and return the vector of normalized exponentiated values.

normalize_weight <- function(logweights){
  mlw <- max(logweights)
  w <- exp(logweights - mlw)
  nw <- w / sum(w) # normalized weights
  return(nw)
}