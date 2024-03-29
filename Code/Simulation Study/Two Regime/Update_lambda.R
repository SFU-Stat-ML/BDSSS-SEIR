

## BDSSSM-SEIR
## Update the precision parameter lambda in PG-CSMC-AS

update.lambda <- function(y, x,             # y_1:T x_0:T
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

  (new.lambda  <- rtruncnorm(1, a=0, b=Inf, mean = lambda, sd = step.size))
  # (new.lambda  <- rgamma(1, shape=lambda^2/step.size^2, rate=lambda/step.size^2))


  log.r <- min(c(log(dtruncnorm(lambda, a=0, b=Inf, mean = new.lambda, sd = step.size)) +
    log.full.conditional(y, x,             # y_1:T x_0:T
                         S, E, I, R,
                         alpha, m.alpha, sigma.alpha,
                         beta, m.beta, sigma.beta,
                         gamma, m.gamma, sigma.gamma,
                         kappa, a.kappa, b.kappa,
                         new.lambda, a.lambda, b.lambda,
                         p, a.p, b.p, m.p, sigma.p,
                         Px, delta.mat,
                         f, a.f, b.f,
                         pop.size) -
    log(dtruncnorm(new.lambda, a=0, b=Inf, mean = lambda, sd = step.size)) -
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

  # (log.r <- min(c(dgamma(lambda, shape=new.lambda^2/step.size^2, rate=new.lambda/step.size^2, log=TRUE) +
  #                  log.full.conditional(y, x,             # y_1:T x_0:T
  #                                       S, E, I, R, 
  #                                       alpha, m.alpha, sigma.alpha,
  #                                       beta, m.beta, sigma.beta,
  #                                       gamma, m.gamma, sigma.gamma,
  #                                       kappa, a.kappa, b.kappa,
  #                                       new.lambda, a.lambda, b.lambda,
  #                                       p, a.p, b.p,
  #                                       Px, delta.mat,
  #                                       f, a.f, b.f, 
  #                                       pop.size)-
  #                   dgamma(new.lambda, shape=lambda^2/step.size^2, rate=lambda/step.size^2, log=TRUE) -
  #                  log.full.conditional(y, x,             # y_1:T x_0:T
  #                                       S, E, I, R, 
  #                                       alpha, m.alpha, sigma.alpha,
  #                                       beta, m.beta, sigma.beta,
  #                                       gamma, m.gamma, sigma.gamma,
  #                                       kappa, a.kappa, b.kappa,
  #                                       lambda, a.lambda, b.lambda,
  #                                       p, a.p, b.p,
  #                                       Px, delta.mat,
  #                                       f, a.f, b.f, 
  #                                       pop.size), log(1))))

  if(log(runif(1)) < log.r){
    new.lambda= new.lambda       # accept move with probability min(1,r)
    indicator = 1                     # indicator of acceptance
  } else{
    new.lambda = lambda          # otherwise "reject" move, and stay where we are
    indicator = 0
  }
  return(list("new.lambda" = new.lambda,
              "indicator" = indicator))


}



