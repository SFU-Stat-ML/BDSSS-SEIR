#############################################################
# Author: Jingxue (Grace) Feng
#         Simon Fraser University, Burnaby, BC, Canada
#         Email: jingxuef@sfu.ca
#############################################################

## BDSSSM-SEIR
## Update the precision parameter kappa in PG-CSMC-AS

library("DirichletReg")

update.kappa <- function(y, x,             # y_1:T x_0:T
                         S, E, I, R,
                         alpha, m.alpha, sigma.alpha,
                         beta, m.beta, sigma.beta,
                         gamma, m.gamma, sigma.gamma,
                         kappa, a.kappa, b.kappa,
                         lambda, a.lambda, b.lambda,
                         p, a.p, b.p, m.p, sigma.p,
                         Px, delta.mat,
                         f, a.f, b.f,
                         pop.size,
                         step.size){
  
  # Proposal
  (new.kappa  <- rtruncnorm(1, a=0, b=Inf, mean = kappa, sd = step.size)) # proposal

  # Acceptance probability
  log.r <- min(c(log(dtruncnorm(kappa, a=0, b=Inf, mean = new.kappa, sd = step.size)) +
                   log.full.conditional(y, x,             # y_1:T x_0:T
                                        S, E, I, R, 
                                        alpha, m.alpha, sigma.alpha,
                                        beta, m.beta, sigma.beta,
                                        gamma, m.gamma, sigma.gamma,
                                        new.kappa, a.kappa, b.kappa,
                                        lambda, a.lambda, b.lambda,
                                        p, a.p, b.p, m.p, sigma.p,
                                        Px, delta.mat,
                                        f, a.f, b.f,
                                        pop.size) -
                   log(dtruncnorm(new.kappa, a=0, b=Inf, mean = kappa, sd = step.size)) -
                   log.full.conditional(y, x,             # y_1:T x_0:T
                                        S, E, I, R, 
                                        alpha, m.alpha, sigma.alpha,
                                        beta, m.beta, sigma.beta,
                                        gamma, m.gamma, sigma.gamma,
                                        kappa, a.kappa, b.kappa,
                                        lambda, a.lambda, b.lambda,
                                        p, a.p, b.p, m.p, sigma.p,
                                        Px, delta.mat,
                                        f, a.f, b.f,
                                        pop.size), log(1)))

  if(log(runif(1)) < log.r){
    new.kappa = new.kappa       # accept move with probability min(1,r)
    indicator = 1                     # indicator of acceptance
  } else{
    new.kappa = kappa          # otherwise "reject" move, and stay where we are
    indicator = 0
  }
  return(list("new.kappa" = new.kappa,
              "indicator" = indicator))

}


# update.kappa <- function(y, x,             # y_1:T x_0:T
#                          S, E, I, R, 
#                          etaS, etaE, etaI, etaR,
#                          alpha, m.alpha, sigma.alpha,
#                          beta, m.beta, sigma.beta,
#                          gamma, m.gamma, sigma.gamma,
#                          kappa, a.kappa, b.kappa,
#                          lambda, a.lambda, b.lambda,
#                          p, a.p, b.p,
#                          Px, delta.mat,
#                          transm.modifier,
#                          pop.size,
#                          step.size){
# 
#   (new.kappa  <- rtruncnorm(1, a=step.size, b=Inf, mean = kappa, sd = step.size))
#   
#   # (new.kappa  <- rgamma(1, shape=kappa^2/step.size^2, rate=kappa/step.size^2))
#   
#   ll <- dgamma(kappa, shape = a.kappa, rate = b.kappa, log=TRUE)
#   ll.new <- dgamma(new.kappa, shape = a.kappa, rate = b.kappa, log=TRUE)
#   
#   T <- length(y)
#   for (t in 1:(T-1)) {
# 
#     ll <- ll + ddirichlet(matrix(c(S[t+1], E[t+1], I[t+1], R[t+1]), nrow = 1),
#                         alpha = kappa*c(etaS[t], etaE[t], etaI[t], etaR[t]),
#                         log=TRUE)
# 
#     ll.new <- ll.new + ddirichlet(matrix(c(S[t+1], E[t+1], I[t+1], R[t+1]), nrow = 1),
#                                   alpha = new.kappa*c(etaS[t], etaE[t], etaI[t], etaR[t]),
#                                   log=TRUE)
#   }
# 
#   log.r <- min(c(log(dtruncnorm(kappa, a=step.size, b=Inf, mean = new.kappa, sd = step.size)) +
#                    ll.new -
#                    log(dtruncnorm(new.kappa, a=step.size, b=Inf, mean = kappa, sd = step.size)) -
#                    ll), log(1))
# 
#   # log.r <- min(c(dgamma(kappa, shape=new.kappa^2/step.size^2, rate=new.kappa/step.size^2, log=TRUE) +
#   #                  ll.new -
#   #                  dgamma(new.kappa, shape=kappa^2/step.size^2, rate=kappa/step.size^2, log=TRUE)  -
#   #                  ll , log(1)))
#  
#   u <- log(runif(1))
#   
#   if(u < log.r){
#     new.kappa = new.kappa       # accept move with probability min(1,r)
#     indicator = 1                     # indicator of acceptance
#   } else{
#     new.kappa = kappa          # otherwise "reject" move, and stay where we are
#     indicator = 0
#   }
#   
#   return(list("new.kappa" = new.kappa,
#               "indicator" = indicator))
# 
# }
# 
