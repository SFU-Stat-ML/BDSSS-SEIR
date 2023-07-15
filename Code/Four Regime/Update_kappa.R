
## Update the precision parameter kappa in PG-CSMC-AS

library("DirichletReg")

update.kappa <- function(y, x,             
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
                   log.full.conditional(y, x,             
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
                   log.full.conditional(y, x,             
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
    indicator = 1               # indicator of acceptance
  } else{
    new.kappa = kappa          # otherwise "reject" move, and stay where we are
    indicator = 0
  }
  return(list("new.kappa" = new.kappa,
              "indicator" = indicator))

}


