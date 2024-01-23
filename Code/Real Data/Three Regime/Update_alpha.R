

## Update the latency parameter alpha in PG-CSMC-AS
## α|θ0:T,x0:T,mα,σ2α∼N(mα,σ2α)Iα>0
source("log_full_conditional.R")

update.alpha <- function(y, x,            # y_1:T x_0:T
                         S, E, I, R, 
                         alpha, m.alpha, sigma.alpha,
                         beta, m.beta, sigma.beta,
                         gamma, m.gamma, sigma.gamma,
                         kappa, a.kappa, b.kappa,
                         lambda, a.lambda, b.lambda,
                         p1, a.p1, b.p1, m.p1, sigma.p1,
                         p2, a.p2, b.p2, m.p2, sigma.p2,
                         Px, delta.mat,
                         f, a.f, b.f,
                         pop.size,
                         step.size){

    new.alpha <- rtruncnorm(1, a=0, b=Inf, mean=alpha, sd=step.size)


    log.r <- min(c(log(dtruncnorm(alpha,  a=0, b=Inf, mean = new.alpha, sd = step.size)) +
                     log.full.conditional(y, x,             # y_1:T x_0:T
                                          S, E, I, R, 
                                          new.alpha, m.alpha, sigma.alpha,
                                          beta, m.beta, sigma.beta,
                                          gamma, m.gamma, sigma.gamma,
                                          kappa, a.kappa, b.kappa,
                                          lambda, a.lambda, b.lambda,
                                          p1, a.p1, b.p1, m.p1, sigma.p1,
                                          p2, a.p2, b.p2, m.p2, sigma.p2,
                                          Px, delta.mat,
                                          f, a.f, b.f,
                                          pop.size) -
                     log(dtruncnorm(new.alpha,  a=0, b=Inf, mean = alpha, sd = step.size)) -
                     log.full.conditional(y, x,             # y_1:T x_0:T
                                          S, E, I, R, 
                                          alpha, m.alpha, sigma.alpha,
                                          beta, m.beta, sigma.beta,
                                          gamma, m.gamma, sigma.gamma,
                                          kappa, a.kappa, b.kappa,
                                          lambda, a.lambda, b.lambda,
                                          p1, a.p1, b.p1, m.p1, sigma.p1,
                                          p2, a.p2, b.p2, m.p2, sigma.p2,
                                          Px, delta.mat,
                                          f, a.f, b.f,
                                          pop.size), log(1)))

    if(log(runif(1)) < log.r){
      new.alpha = new.alpha       # accept move with probability min(1,r)
      indicator = 1                     # indicator of acceptance
    } else{
      new.alpha = alpha         # otherwise "reject" move, and stay where we are
      indicator = 0
    }


  return(list("new.alpha" = new.alpha,
              "indicator" = indicator))

}

# update.alpha <- function(y, x,             # y_1:T x_0:T
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
#   new.alpha <- rtruncnorm(1, a=0, b=Inf, mean=alpha, sd=step.size)
# 
#   log.r <- min(c(log(dtruncnorm(alpha, a=0, b=Inf, mean = new.alpha, sd = step.size)) +
#                    log(dtruncnorm(new.alpha, a=0, b=Inf, mean = m.alpha, sd = sigma.alpha))-
#                    log(dtruncnorm(new.alpha, a=0, b=Inf, mean = alpha, sd = step.size)) -
#                    log(dtruncnorm(alpha, a=0, b=Inf, mean = m.alpha, sd = sigma.alpha)), log(1)))
# 
#   if(log(runif(1)) < log.r){
#     new.alpha = new.alpha       # accept move with probability min(1,r)
#     indicator = 1                     # indicator of acceptance
#   } else{
#     new.alpha = alpha         # otherwise "reject" move, and stay where we are
#     indicator = 0
#   }
# 
#   return(list("new.alpha" = new.alpha,
#               "indicator" = indicator))
# 
# }

