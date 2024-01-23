
## BDSSSM-SEIR
## Update the recovery parameter gamma in PG-CSMC-AS
## γ|θ0:T,x0:T,mγ,γ2β∼N(mγ,σ2γ)Iγ>0
source("log_full_conditional.R")

update.gamma <- function(y, x,             # y_1:T x_0:T
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

  new.gamma <- rtruncnorm(1, a=0, b=Inf, mean = gamma, sd = step.size)

    log.r <- min(c(log(dtruncnorm(gamma, a=0, b=Inf, mean = new.gamma, sd = step.size)) +
                     log.full.conditional(y, x,             # y_1:T x_0:T
                                          S, E, I, R, 
                                          alpha, m.alpha, sigma.alpha,
                                          beta, m.beta, sigma.beta,
                                          new.gamma, m.gamma, sigma.gamma,
                                          kappa, a.kappa, b.kappa,
                                          lambda, a.lambda, b.lambda,
                                          p, a.p, b.p, m.p, sigma.p,
                                          Px, delta.mat,
                                          f, a.f, b.f,
                                          pop.size)-
                     log(dtruncnorm(new.gamma, a=0, b=Inf, mean = gamma, sd = step.size)) -
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
      new.gamma = new.gamma       # accept move with probability min(1,r)
      indicator = 1                     # indicator of acceptance
    } else{
      new.gamma = gamma         # otherwise "reject" move, and stay where we are
      indicator = 0
    }



  return(list("new.gamma" = new.gamma,
              "indicator" = indicator))

}


# 
# update.gamma <- function(y, x,             # y_1:T x_0:T
#                          S, E, I, R, S0, E0, I0, R0,
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
#   new.gamma <- rtruncnorm(1, a=0, b=Inf, mean=gamma, sd=step.size)
# 
#   log.r <- min(c(log(dtruncnorm(gamma, a=0, b=Inf, mean = new.gamma, sd = step.size)) +
#                    log(dtruncnorm(new.gamma, a=0, b=Inf, mean = m.gamma, sd = sigma.gamma))-
#                    log(dtruncnorm(new.gamma, a=0, b=Inf, mean = gamma, sd = step.size)) -
#                    log(dtruncnorm(gamma, a=0, b=Inf, mean = m.gamma, sd = sigma.gamma)), log(1)))
# 
#   if(log(runif(1)) < log.r){
#     new.gamma = new.gamma      # accept move with probability min(1,r)
#     indicator = 1                     # indicator of acceptance
#   } else{
#     new.gamma = gamma         # otherwise "reject" move, and stay where we are
#     indicator = 0
#   }
# 
#   return(list("new.gamma" = new.gamma,
#               "indicator" = indicator))
# 
# }
