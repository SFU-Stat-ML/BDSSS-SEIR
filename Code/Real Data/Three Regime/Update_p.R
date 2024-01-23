

## Update the detection rate p in PG-CSMC-AS

update.p1 <- function(y, x,             # y_1:T x_0:T
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
  
  new.p1 <- rtruncnorm(1, a=a.p1, b=b.p1, mean=p1, sd=step.size)
  
  log.r <- min(c(log(dtruncnorm(p1, a=a.p1, b=b.p1, mean = new.p1, sd = step.size)) + 
                   log.full.conditional(y, x,             # y_1:T x_0:T
                                        S, E, I, R, 
                                        alpha, m.alpha, sigma.alpha,
                                        beta, m.beta, sigma.beta,
                                        gamma, m.gamma, sigma.gamma,
                                        kappa, a.kappa, b.kappa, 
                                        lambda, a.lambda, b.lambda,
                                        new.p1, a.p1, b.p1, m.p1, sigma.p1,
                                        p2, a.p2, b.p2, m.p2, sigma.p2,
                                        Px, delta.mat,
                                        f, a.f, b.f,
                                        pop.size)-
                   log(dtruncnorm(new.p1, a=a.p1, b=b.p1, mean = p1, sd = step.size)) -
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
    new.p1 = new.p1     # accept move with probability min(1,r)
    indicator = 1                     # indicator of acceptance
  } else{
    new.p1 = p1         # otherwise "reject" move, and stay where we are
    indicator = 0
  }
  return(list("new.p1" = new.p1,
              "indicator" = indicator))
  
}




update.p2 <- function(y, x,             # y_1:T x_0:T
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
  
  new.p2 <- rtruncnorm(1, a=a.p2, b=b.p2, mean=p2, sd=step.size)
  
  log.r <- min(c(log(dtruncnorm(p2, a=a.p2, b=b.p2, mean = new.p2, sd = step.size)) + 
                   log.full.conditional(y, x,             # y_1:T x_0:T
                                        S, E, I, R, 
                                        alpha, m.alpha, sigma.alpha,
                                        beta, m.beta, sigma.beta,
                                        gamma, m.gamma, sigma.gamma,
                                        kappa, a.kappa, b.kappa, 
                                        lambda, a.lambda, b.lambda,
                                        p1, a.p1, b.p1, m.p1, sigma.p1,
                                        new.p2, a.p2, b.p2, m.p2, sigma.p2,
                                        Px, delta.mat,
                                        f, a.f, b.f,
                                        pop.size)-
                   log(dtruncnorm(new.p2, a=a.p2, b=b.p2, mean = p2, sd = step.size)) -
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
    new.p2 = new.p2     # accept move with probability min(1,r)
    indicator = 1                     # indicator of acceptance
  } else{
    new.p2 = p2         # otherwise "reject" move, and stay where we are
    indicator = 0
  }
  return(list("new.p2" = new.p2,
              "indicator" = indicator))
  
}




