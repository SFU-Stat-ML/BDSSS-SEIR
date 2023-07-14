

## Update the transition probability matrix (row vector) 
source("log_full_conditional.R")

update.pi.k <- function(y, x,             
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
  pi.k <- Px[k,]          

  if(k==1){
    new.pi.k.11 <- rtruncnorm(1, a=0, b=1, mean=pi.k[1], sd=step.size)
    new.pi.k.12 <- 1-new.pi.k.11
    new.pi.k <- c(new.pi.k.11,  new.pi.k.12)

    new.Px <- Px
    new.Px[k,] <- new.pi.k                  

    # compute log of the acceptance ratio
    log.r <- min(c(log(dtruncnorm(pi.k[1], mean = new.pi.k.11, sd = step.size)) +
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
                     log(dtruncnorm(new.pi.k.11, mean = pi.k[1], sd = step.size))-
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

  }else if(k==2){
    new.pi.k.22 <- rtruncnorm(1, a=0, b=1, mean=pi.k[2], sd=step.size)
    new.pi.k.21 <- 1-new.pi.k.22
    new.pi.k <- c(new.pi.k.21,  new.pi.k.22)

    new.Px <- Px
    new.Px[k,] <- new.pi.k                

    # compute log of the acceptance ratio
    log.r <- min(c(log(dtruncnorm(pi.k[2], mean = new.pi.k.22, sd = step.size)) +
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
                     log(dtruncnorm(new.pi.k.22, mean = pi.k[2], sd = step.size)) -
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

