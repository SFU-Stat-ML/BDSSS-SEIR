

## Update the transmission modifier in PG-CSMC-AS

update.f <- function(y, x,             # y_1:T x_0:T
                      S, E, I, R, 
                      alpha, m.alpha, sigma.alpha,
                      beta, m.beta, sigma.beta,
                      gamma, m.gamma, sigma.gamma,
                      kappa, a.kappa, b.kappa, 
                      lambda, a.lambda, b.lambda,
                      p, a.p, b.p, m.p, sigma.p,
                      Px, delta.mat, 
                      f, a.f, b.f, # a.f, b.f are vectors of length >=1
                      pop.size,
                      step.size){
  
  K <- length(f)
  
  new.f <- f
  log.r <- matrix(NA, nrow=K, ncol=1)
  indicator <- matrix(NA, nrow=K, ncol=1) # no indicator for f1=1
  
  # Given f1=1, update f2
  for (k in 2:K){
    
    new.f[k] <- rtruncnorm(1, a=a.f, b=b.f, mean=f[k], sd=step.size)
    log.r[k,1] <- min(c(sum(log(dtruncnorm(f[k], a=a.f, b=b.f, mean = new.f[k], sd = step.size))) + 
                   log.full.conditional(y, x,             # y_1:T x_0:T
                                        S, E, I, R, 
                                        alpha, m.alpha, sigma.alpha,
                                        beta, m.beta, sigma.beta,
                                        gamma, m.gamma, sigma.gamma,
                                        kappa, a.kappa, b.kappa, 
                                        lambda, a.lambda, b.lambda,
                                        p, a.p, b.p, m.p, sigma.p,
                                        Px, delta.mat,
                                        new.f, a.f, b.f,
                                        pop.size)-
                   sum(log(dtruncnorm(new.f[k], a=a.f, b=b.f, mean = f[k], sd = step.size))) -
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
    if(log(runif(1)) < log.r[k,1]){
      new.f = new.f     # accept move with probability min(1,r)
      indicator[k,1] = 1                     # indicator of acceptance
    } else{
      new.f = f         # otherwise "reject" move, and stay where we are
      indicator[k,1] = 0
    }
    
    
    }

  return(list("new.f" = new.f,
              "indicator" = indicator))
  
}
