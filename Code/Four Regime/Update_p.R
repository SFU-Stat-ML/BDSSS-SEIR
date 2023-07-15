

## Update the detection rate p in PG-CSMC-AS

update.p <- function(y, x,             
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
  
  new.p <- rtruncnorm(1, a=a.p, b=b.p, mean=p, sd=step.size)
  
  log.r <- min(c(log(dtruncnorm(p, a=a.p, b=b.p, mean = new.p, sd = step.size)) + 
                   log.full.conditional(y, x,             
                                        S, E, I, R, 
                                        alpha, m.alpha, sigma.alpha,
                                        beta, m.beta, sigma.beta,
                                        gamma, m.gamma, sigma.gamma,
                                        kappa, a.kappa, b.kappa, 
                                        lambda, a.lambda, b.lambda,
                                        new.p, a.p, b.p,
                                        Px, delta.mat,
                                        f, a.f, b.f,
                                        pop.size)-
                   log(dtruncnorm(new.p, a=a.p, b=b.p, mean = p, sd = step.size)) -
                   log.full.conditional(y, x,             
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
    new.p = new.p     # accept move with probability min(1,r)
    indicator = 1     # indicator of acceptance
  } else{
    new.p = p         # otherwise "reject" move, and stay where we are
    indicator = 0
  }
  return(list("new.p" = new.p,
              "indicator" = indicator))
  
}
