
# Transition Density

# θt|θt−1, xt ∼ Dirichlet(κ r(θt−1; α, β, γ, fxt ))

trans.density <- function(St, Et, It, Rt, 
                          kappa, 
                          eta.S, eta.E, eta.I, eta.R){
  
  return(DirichletReg::ddirichlet(matrix(c(St, Et, It, Rt),nrow=1),
                                  alpha = kappa * c(eta.S, eta.E, eta.I, eta.R)))
  
}


log.trans.density <- function(St, Et, It, Rt, 
                          kappa, 
                          eta.S, eta.E, eta.I, eta.R){
  
  return(DirichletReg::ddirichlet(matrix(c(St, Et, It, Rt),nrow=1),
                                  alpha = kappa * c(eta.S, eta.E, eta.I, eta.R),
                                  log = TRUE))
  
}