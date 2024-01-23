

## Update the transition probability matrix (row vector) in PG-CSMC-AS

library("DirichletReg")

## Metropolis-hasting step for Px

# For three regimes
update.pi.k <- function(y, x,             # y_1:T x_0:T
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

  k <- ceiling(runif(1)*length(f))
  pi.k <- Px[k,]          # \pi_k

  if(k==1){
    new.pi.k.11 <- rtruncnorm(1, a=0, b=1, mean=pi.k[1], sd=step.size)
    new.pi.k.12 <- rtruncnorm(1, a=0, b=1-new.pi.k.11, mean=pi.k[2], sd=step.size)
    new.pi.k.13 <- 1-new.pi.k.11-new.pi.k.12
    new.pi.k <- c(new.pi.k.11,  new.pi.k.12, new.pi.k.13)

    new.Px <- Px
    new.Px[k,] <- new.pi.k                  # Replace the k-th row with new pi

    # compute log of the acceptance ratio
    log.r <- min(c(log(dtruncnorm(pi.k[1], a=0, b=1, mean = new.pi.k.11, sd = step.size)) +
                     log(dtruncnorm(pi.k[2], a=0, b=1-new.pi.k.11, mean = new.pi.k.12, sd = step.size)) +
                     log.full.conditional(y, x,             # y_1:T x_0:T
                                          S, E, I, R,
                                          alpha, m.alpha, sigma.alpha,
                                          beta, m.beta, sigma.beta,
                                          gamma, m.gamma, sigma.gamma,
                                          kappa, a.kappa, b.kappa,
                                          lambda, a.lambda, b.lambda,
                                          p, a.p, b.p, m.p, sigma.p,
                                          new.Px, delta.mat,
                                          f, a.f, b.f,
                                          pop.size)-
                     log(dtruncnorm(new.pi.k.11, a=0, b=1, mean = pi.k[1], sd = step.size))-
                     log(dtruncnorm(new.pi.k.12, a=0, b=1-new.pi.k.11, mean = pi.k[2], sd = step.size))-
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
      new.pi.k = new.pi.k       # accept move with probabily min(1,r)
      indicator = 1             # indicator of acceptance
    } else{
      new.pi.k = pi.k           # otherwise "reject" move, and stay where we are
      indicator = 0
    }

    new.Px[k,] <- new.pi.k

  }else if(k==2){
    new.pi.k.21 <- rtruncnorm(1, a=0, b=1, mean=pi.k[1], sd=step.size)
    new.pi.k.22 <- rtruncnorm(1, a=0, b=1-new.pi.k.21, mean=pi.k[2], sd=step.size)
    new.pi.k.23 <- 1-new.pi.k.21-new.pi.k.22
    new.pi.k <- c(new.pi.k.21,  new.pi.k.22, new.pi.k.23)

    new.Px <- Px
    new.Px[k,] <- new.pi.k                  # Replace the k-th row with new pi

    # compute log of the acceptance ratio
    log.r <- min(c(log(dtruncnorm(pi.k[1], a=0, b=1, mean = new.pi.k.21, sd = step.size)) +
                     log(dtruncnorm(pi.k[2], a=0, b=1-new.pi.k.21, mean = new.pi.k.22, sd = step.size)) +
                     log.full.conditional(y, x,             
                                          S, E, I, R,
                                          alpha, m.alpha, sigma.alpha,
                                          beta, m.beta, sigma.beta,
                                          gamma, m.gamma, sigma.gamma,
                                          kappa, a.kappa, b.kappa,
                                          lambda, a.lambda, b.lambda,
                                          p, a.p, b.p, m.p, sigma.p,
                                          new.Px, delta.mat,
                                          f, a.f, b.f,
                                          pop.size)-
                     log(dtruncnorm(new.pi.k.21, a=0, b=1, mean = pi.k[1], sd = step.size)) -
                     log(dtruncnorm(new.pi.k.22, a=0, b=1-new.pi.k.21, mean = pi.k[2], sd = step.size)) -
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
      new.pi.k = new.pi.k       # accept move with probabily min(1,r)
      indicator = 1             # indicator of acceptance
    } else{
      new.pi.k = pi.k           # otherwise "reject" move, and stay where we are
      indicator = 0
    }

    new.Px[k,] <- new.pi.k
  }else if(k==3){
    new.pi.k.31 <- rtruncnorm(1, a=0, b=1, mean=pi.k[1], sd=step.size)
    new.pi.k.32 <- rtruncnorm(1, a=0, b=1-new.pi.k.31, mean=pi.k[2], sd=step.size)
    new.pi.k.33 <- 1-new.pi.k.31-new.pi.k.32
    new.pi.k <- c(new.pi.k.31,  new.pi.k.32, new.pi.k.33)

    new.Px <- Px
    new.Px[k,] <- new.pi.k                  # Replace the k-th row with new pi

    # compute log of the acceptance ratio
    log.r <- min(c(log(dtruncnorm(pi.k[1], a=0, b=1, mean = new.pi.k.31, sd = step.size)) +
                     log(dtruncnorm(pi.k[2], a=0, b=1-new.pi.k.31, mean = new.pi.k.32, sd = step.size)) +
                     log.full.conditional(y, x,             
                                          S, E, I, R,
                                          alpha, m.alpha, sigma.alpha,
                                          beta, m.beta, sigma.beta,
                                          gamma, m.gamma, sigma.gamma,
                                          kappa, a.kappa, b.kappa,
                                          lambda, a.lambda, b.lambda,
                                          p, a.p, b.p, m.p, sigma.p,
                                          new.Px, delta.mat,
                                          f, a.f, b.f,
                                          pop.size)-
                     log(dtruncnorm(new.pi.k.31, a=0, b=1, mean = pi.k[1], sd = step.size)) -
                     log(dtruncnorm(new.pi.k.32, a=0, b=1-new.pi.k.31, mean = pi.k[2], sd = step.size)) -
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
      new.pi.k = new.pi.k       # accept move with probabily min(1,r)
      indicator = 1             # indicator of acceptance
    } else{
      new.pi.k = pi.k           # otherwise "reject" move, and stay where we are
      indicator = 0
    }

    new.Px[k,] <- new.pi.k
  }


  return(list("indicator" = indicator,
              "k" = k,
              "newPx" = new.Px))
}


