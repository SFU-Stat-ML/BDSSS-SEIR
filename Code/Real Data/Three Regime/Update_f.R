
#############################################################
# Author: Jingxue (Grace) Feng
#         Simon Fraser University, Burnaby, BC, Canada
#         Email: jingxuef@sfu.ca
#############################################################

## Update the transmission modifier in PG-CSMC-AS

# update.f <- function(y, x,             # y_1:T x_0:T
#                       S, E, I, R,
#                       alpha, m.alpha, sigma.alpha,
#                       beta, m.beta, sigma.beta,
#                       gamma, m.gamma, sigma.gamma,
#                       kappa, a.kappa, b.kappa,
#                       lambda, a.lambda, b.lambda,
#                       p, a.p, b.p,
#                       Px, delta.mat,
#                       f, a.f, b.f, # a.f, b.f are vectors of length >=1
#                       pop.size,
#                       step.size){
# 
#   K <- length(f)
# 
#   new.f <- f
#   log.r <- matrix(NA, nrow=K, ncol=1)
#   indicator <- matrix(NA, nrow=K, ncol=1) # no indicator for f1=1
# 
#   # Given f1=1, update f2, f3,...,fK
#   for (k in 2:K){
# 
#     new.f[k] <- rtruncnorm(1, a=a.f[k-1], b=b.f[k-1], mean=f[k], sd=step.size)
#     log.r[k,1] <- min(c(sum(log(dtruncnorm(f[k], a=a.f[k-1], b=b.f[k-1], mean = new.f[k], sd = step.size))) +
#                    log.full.conditional(y, x,             # y_1:T x_0:T
#                                         S, E, I, R,
#                                         alpha, m.alpha, sigma.alpha,
#                                         beta, m.beta, sigma.beta,
#                                         gamma, m.gamma, sigma.gamma,
#                                         kappa, a.kappa, b.kappa,
#                                         lambda, a.lambda, b.lambda,
#                                         p, a.p, b.p,
#                                         Px, delta.mat,
#                                         new.f, a.f, b.f,
#                                         pop.size)-
#                    sum(log(dtruncnorm(new.f[k], a=a.f[k-1], b=b.f[k-1], mean = f[k], sd = step.size))) -
#                    log.full.conditional(y, x,             # y_1:T x_0:T
#                                         S, E, I, R,
#                                         alpha, m.alpha, sigma.alpha,
#                                         beta, m.beta, sigma.beta,
#                                         gamma, m.gamma, sigma.gamma,
#                                         kappa, a.kappa, b.kappa,
#                                         lambda, a.lambda, b.lambda,
#                                         p, a.p, b.p,
#                                         Px, delta.mat,
#                                         f, a.f, b.f,
#                                         pop.size), log(1)))
#     if(log(runif(1)) < log.r[k,1]){
#       new.f = new.f     # accept move with probability min(1,r)
#       indicator[k,1] = 1                     # indicator of acceptance
#     } else{
#       new.f = f         # otherwise "reject" move, and stay where we are
#       indicator[k,1] = 0
#     }
# 
# 
#     }
# 
#   return(list("new.f" = new.f,
#               "indicator" = indicator))
# 
# }



update.f2 <- function(y, x,             # y_1:T x_0:T
                      S, E, I, R,
                      alpha, m.alpha, sigma.alpha,
                      beta, m.beta, sigma.beta,
                      gamma, m.gamma, sigma.gamma,
                      kappa, a.kappa, b.kappa,
                      lambda, a.lambda, b.lambda,
                      p1, a.p1, b.p1, m.p1, sigma.p1,
                      p2, a.p2, b.p2, m.p2, sigma.p2,
                      Px, delta.mat,
                      f, a.f, b.f, # a.f, b.f are vectors of length >=1
                      pop.size,
                      step.size){

  new.f2 <- rtruncnorm(1, a=a.f[1], b=b.f[1], mean=f[2], sd=step.size)
  new.f <- c(f[1],  new.f2 , f[3])
  
  log.r <- min(c(log(dtruncnorm(f[2], a=a.f[1], b=b.f[1], mean = new.f2, sd = step.size)) +
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
                                        new.f, a.f, b.f,
                                        pop.size)-
                   log(dtruncnorm(new.f2, a=a.f[1], b=b.f[1], mean = f[2], sd = step.size)) -
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
      new.f2 = new.f2     # accept move with probability min(1,r)
      indicator = 1                     # indicator of acceptance
    } else{
      new.f2 = f[2]         # otherwise "reject" move, and stay where we are
      indicator = 0
    }

  return(list("new.f2" = new.f2,
              "indicator" = indicator))

}



update.f3 <- function(y, x,             # y_1:T x_0:T
                      S, E, I, R,
                      alpha, m.alpha, sigma.alpha,
                      beta, m.beta, sigma.beta,
                      gamma, m.gamma, sigma.gamma,
                      kappa, a.kappa, b.kappa,
                      lambda, a.lambda, b.lambda,
                      p1, a.p1, b.p1, m.p1, sigma.p1,
                      p2, a.p2, b.p2, m.p2, sigma.p2,
                      Px, delta.mat,
                      f, a.f, b.f, # a.f, b.f are vectors of length >=1
                      pop.size,
                      step.size){
  
  new.f3 <- rtruncnorm(1, a=a.f[2], b=b.f[2], mean=f[3], sd=step.size)
  new.f <- c(f[1],  f[2] , new.f3)
  
  log.r <- min(c(sum(log(dtruncnorm(f[3], a=a.f[2], b=b.f[2], mean = new.f3, sd = step.size))) +
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
                                        new.f, a.f, b.f,
                                        pop.size)-
                   sum(log(dtruncnorm(new.f3, a=a.f[2], b=b.f[2], mean = f[3], sd = step.size))) -
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
    new.f3 = new.f3     # accept move with probability min(1,r)
    indicator = 1                     # indicator of acceptance
  } else{
    new.f3 = f[3]         # otherwise "reject" move, and stay where we are
    indicator = 0
  }
  
  return(list("new.f3" = new.f3,
              "indicator" = indicator))
  
}
