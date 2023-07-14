#############################################################
# Author: Jingxue (Grace) Feng
#         Simon Fraser University, Burnaby, BC, Canada
#         Email: jingxuef@sfu.ca
#############################################################

## Update the transmission parameter beta in PG-CSMC-AS
## β|θ0:T,x0:T,mβ,σ2β∼N(mβ,σ2β)Iβ>0
source("log_full_conditional.R")

update.beta <- function(y, x,             # y_1:T x_0:T
                        S, E, I, R, 
                        alpha, m.alpha, sigma.alpha,
                        beta, m.beta, sigma.beta,
                        gamma, m.gamma, sigma.gamma,
                        kappa, a.kappa, b.kappa,
                        lambda, a.lambda, b.lambda,
                        p, a.p, b.p,
                        Px, delta.mat,
                        f, a.f, b.f,
                        pop.size,
                        step.size){

  new.beta <- rtruncnorm(1, a=0, b=Inf, mean = beta, sd = step.size)

    log.r <- min(c(log(dtruncnorm(beta, a=0, b=Inf, mean = new.beta, sd = step.size)) +
                     log.full.conditional(y, x,             # y_1:T x_0:T
                                          S, E, I, R, 
                                          alpha, m.alpha, sigma.alpha,
                                          new.beta, m.beta, sigma.beta,
                                          gamma, m.gamma, sigma.gamma,
                                          kappa, a.kappa, b.kappa,
                                          lambda, a.lambda, b.lambda,
                                          p, a.p, b.p,
                                          Px, delta.mat,
                                          f, a.f, b.f,
                                          pop.size)-
                     log(dtruncnorm(new.beta, a=0, b=Inf, mean = beta, sd = step.size)) -
                     log.full.conditional(y, x,             # y_1:T x_0:T
                                          S, E, I, R,
                                          alpha, m.alpha, sigma.alpha,
                                          beta, m.beta, sigma.beta,
                                          gamma, m.gamma, sigma.gamma,
                                          kappa, a.kappa, b.kappa,
                                          lambda, a.lambda, b.lambda,
                                          p, a.p, b.p,
                                          Px, delta.mat,
                                          f, a.f, b.f,
                                          pop.size), log(1)))

    if(log(runif(1)) < log.r){
      new.beta = new.beta       # accept move with probability min(1,r)
      indicator = 1                     # indicator of acceptance
    } else{
      new.beta = beta         # otherwise "reject" move, and stay where we are
      indicator = 0
    }

  return(list("new.beta" = new.beta,
              "indicator" = indicator))

}

# update.beta <- function(y, x,             # y_1:T x_0:T
#                         S, E, I, R,
#                         alpha, m.alpha, sigma.alpha,
#                         beta, m.beta, sigma.beta,
#                         gamma, m.gamma, sigma.gamma,
#                         kappa, a.kappa, b.kappa,
#                         lambda, a.lambda, b.lambda,
#                         p, a.p, b.p,
#                         Px, delta.mat,
#                         transm.modifier,
#                         step.size){
#
#   new.beta <- rnorm(1, mean = beta, sd = step.size)
#
#   if (new.beta > 0){
#     log.r <- min(c(dnorm(beta, mean = new.beta, sd = step.size, log=TRUE) +
#                      log.full.conditional(y, x,             # y_1:T x_0:T
#                                           S, E, I, R,
#                                           alpha, m.alpha, sigma.alpha,
#                                           new.beta, m.beta, sigma.beta,
#                                           gamma, m.gamma, sigma.gamma,
#                                           kappa, a.kappa, b.kappa,
#                                           lambda, a.lambda, b.lambda,
#                                           p, a.p, b.p,
#                                           Px, delta.mat,
#                                           transm.modifier)-
#                      dnorm(new.beta, mean = beta, sd = step.size, log=TRUE) -
#                      log.full.conditional(y, x,             # y_1:T x_0:T
#                                           S, E, I, R,
#                                           alpha, m.alpha, sigma.alpha,
#                                           beta, m.beta, sigma.beta,
#                                           gamma, m.gamma, sigma.gamma,
#                                           kappa, a.kappa, b.kappa,
#                                           lambda, a.lambda, b.lambda,
#                                           p, a.p, b.p,
#                                           Px, delta.mat,
#                                           transm.modifier), log(1)))
#
#     if(log(runif(1)) < log.r){
#       new.beta = new.beta       # accept move with probability min(1,r)
#       indicator = 1                     # indicator of acceptance
#     } else{
#       new.beta = beta         # otherwise "reject" move, and stay where we are
#       indicator = 0
#     }
#
#   }else{
#     new.beta = beta
#     indicator = 0
#   }
#
#   return(list("new.beta" = new.beta,
#               "indicator" = indicator))
#
# }

