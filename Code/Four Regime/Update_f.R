

## Update the transmission modifier in PG-CSMC-AS

## For four regimes

update.f2 <- function(y, x,            
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

  new.f2 <- rtruncnorm(1, a=a.f[1], b=b.f[1], mean=f[2], sd=step.size)
  new.f <- c(f[1],  new.f2 , f[3], f[4])
  
  log.r <- min(c(log(dtruncnorm(f[2], a=a.f[1], b=b.f[1], mean = new.f2, sd = step.size)) +
                   log.full.conditional(y, x,             
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
                   log(dtruncnorm(new.f2, a=a.f[1], b=b.f[1], mean = f[2], sd = step.size)) -
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
      new.f2 = new.f2     # accept move with probability min(1,r)
      indicator = 1       # indicator of acceptance
    } else{
      new.f2 = f[2]       # otherwise "reject" move, and stay where we are
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
                      p, a.p, b.p, m.p, sigma.p,
                      Px, delta.mat,
                      f, a.f, b.f, # a.f, b.f are vectors of length >=1
                      pop.size,
                      step.size){
  
  new.f3 <- rtruncnorm(1, a=a.f[2], b=b.f[2], mean=f[3], sd=step.size)
  new.f <- c(f[1],  f[2] , new.f3, f[4])
  
  log.r <- min(c(sum(log(dtruncnorm(f[3], a=a.f[2], b=b.f[2], mean = new.f3, sd = step.size))) +
                   log.full.conditional(y, x,           
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
                   sum(log(dtruncnorm(new.f3, a=a.f[2], b=b.f[2], mean = f[3], sd = step.size))) -
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
    new.f3 = new.f3     # accept move with probability min(1,r)
    indicator = 1       # indicator of acceptance
  } else{
    new.f3 = f[3]       # otherwise "reject" move, and stay where we are
    indicator = 0
  }
  
  return(list("new.f3" = new.f3,
              "indicator" = indicator))
  
}


update.f4 <- function(y, x,            
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
  
  new.f4 <- rtruncnorm(1, a=a.f[3], b=b.f[3], mean=f[4], sd=step.size)
  new.f <- c(f[1],  f[2] , f[3], new.f4)
  
  log.r <- min(c(sum(log(dtruncnorm(f[4], a=a.f[3], b=b.f[3], mean = new.f4, sd = step.size))) +
                   log.full.conditional(y, x,             
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
                   sum(log(dtruncnorm(new.f4, a=a.f[3], b=b.f[3], mean = f[4], sd = step.size))) -
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
    new.f4 = new.f4     # accept move with probability min(1,r)
    indicator = 1       # indicator of acceptance
  } else{
    new.f4 = f[4]       # otherwise "reject" move, and stay where we are
    indicator = 0
  }
  
  return(list("new.f4" = new.f4,
              "indicator" = indicator))
  
}
