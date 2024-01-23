

## Update the transition probability matrix (row vector) in PG-CSMC-AS

library("DirichletReg")

## Metropolis-hasting step for Px
# a.k: scale parameter
# k: index of state
# Delta.k: Dirichlet prior parameter of pi_k
# currentPx: transition probability matrix
# RefParticleX: reference trajectory of x (discrete switching states)

# update.pi.k <- function(y, x,             # y_1:T x_0:T
#                         S, E, I, R, 
#                         alpha, m.alpha, sigma.alpha,
#                         beta, m.beta, sigma.beta,
#                         gamma, m.gamma, sigma.gamma,
#                         kappa, a.kappa, b.kappa, 
#                         lambda, a.lambda, b.lambda,
#                         p, a.p, b.p,
#                         Px, delta.mat, 
#                         f, a.f, b.f,
#                         pop.size,
#                         scale.param){
#   
#   k <- ceiling(runif(1)*length(f))
#   pi.k <- Px[k,]          # \pi_k
#   new.pi.k <- DirichletReg::rdirichlet(1, scale.param*pi.k) # proposal distribution
#   new.Px <- Px
#   new.Px[k,] <- new.pi.k                  # Replace the k-th row with new pi
#   
#   log.r <- min(c(DirichletReg::ddirichlet(matrix(pi.k,1,length(pi.k)), scale.param*new.pi.k, log=TRUE) + 
#                    log.full.conditional(y, x,             # y_1:T x_0:T
#                                         S, E, I, R, 
#                                         alpha, m.alpha, sigma.alpha,
#                                         beta, m.beta, sigma.beta,
#                                         gamma, m.gamma, sigma.gamma,
#                                         kappa, a.kappa, b.kappa, 
#                                         lambda, a.lambda, b.lambda,
#                                         p, a.p, b.p,
#                                         new.Px, delta.mat, 
#                                         f, a.f, b.f,
#                                         pop.size)-
#                    DirichletReg::ddirichlet(matrix(new.pi.k,1,length(new.pi.k)), scale.param*pi.k, log=TRUE) -
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
#   
#   if(log(runif(1)) < log.r){
#     new.pi.k = new.pi.k       # accept move with probabily min(1,r)
#     indicator = 1    # indicator of acceptance
#   } else{
#     new.pi.k = pi.k        # otherwise "reject" move, and stay where we are
#     indicator = 0
#   }
#   
#   new.Px[k,] <- new.pi.k
#   
#   return(list("indicator" = indicator,
#               "k" = k,
#               "newPx" = new.Px))
# }


## Metropolis-hasting step for Px
# Gaussian proposal- generates a new proposal for the k-th element of Px by drawing 
# from a normal distribution with mean pi.k and standard deviation of 0.1. 
# We take the absolute value of the result to ensure that the proposal lies between 0 and 1.

# update.pi.k <- function(y, x,             # y_1:T x_0:T
#                         S, E, I, R, 
#                         alpha, m.alpha, sigma.alpha,
#                         beta, m.beta, sigma.beta,
#                         gamma, m.gamma, sigma.gamma,
#                         kappa, a.kappa, b.kappa, 
#                         lambda, a.lambda, b.lambda,
#                         p, a.p, b.p,
#                         Px, delta.mat, 
#                         f, a.f, b.f,
#                         pop.size,
#                         step.size){
#   
#   k <- ceiling(runif(1)*length(f))
#   pi.k <- Px[k,]          # \pi_k
#   new.pi.k <- abs(rnorm(length(pi.k), mean=pi.k, sd=step.size)) # proposal distribution
#   new.pi.k <- new.pi.k / sum(new.pi.k)  # rescale to ensure the sum is 1
#   
#   new.Px <- Px
#   new.Px[k,] <- new.pi.k                  # Replace the k-th row with new pi
#   
#   # compute log of the acceptance ratio
#   log.r <- min(c(sum(dnorm(pi.k, mean = new.pi.k, sd = step.size, log = TRUE)) + 
#                    log.full.conditional(y, x,             # y_1:T x_0:T
#                                         S, E, I, R, 
#                                         alpha, m.alpha, sigma.alpha,
#                                         beta, m.beta, sigma.beta,
#                                         gamma, m.gamma, sigma.gamma,
#                                         kappa, a.kappa, b.kappa, 
#                                         lambda, a.lambda, b.lambda,
#                                         p, a.p, b.p,
#                                         new.Px, delta.mat, 
#                                         f, a.f, b.f,
#                                         pop.size)-
#                    sum(dnorm(new.pi.k, mean = pi.k, sd = step.size, log = TRUE)) -
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
#   
#   if(log(runif(1)) < log.r){
#     new.pi.k = new.pi.k       # accept move with probabily min(1,r)
#     indicator = 1    # indicator of acceptance
#   } else{
#     new.pi.k = pi.k        # otherwise "reject" move, and stay where we are
#     indicator = 0
#   }
#   
#   new.Px[k,] <- new.pi.k
#   
#   return(list("indicator" = indicator,
#               "k" = k,
#               "newPx" = new.Px))
# }



update.pi.k <- function(y, x,             # y_1:T x_0:T
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

  k <- ceiling(runif(1)*length(f))
  pi.k <- Px[k,]          # \pi_k
  # new.pi.k <- abs(rnorm(length(pi.k), mean=pi.k, sd=step.size)) # proposal distribution
  # new.pi.k <- new.pi.k / sum(new.pi.k)  # rescale to ensure the sum is 1

  if(k==1){
    new.pi.k.11 <- rtruncnorm(1, a=0, b=1, mean=pi.k[1], sd=step.size)
    new.pi.k.12 <- 1-new.pi.k.11
    new.pi.k <- c(new.pi.k.11,  new.pi.k.12)

    new.Px <- Px
    new.Px[k,] <- new.pi.k                  # Replace the k-th row with new pi

    # compute log of the acceptance ratio
    log.r <- min(c(log(dtruncnorm(pi.k[1], mean = new.pi.k.11, sd = step.size)) +
                     log.full.conditional(y, x,             # y_1:T x_0:T
                                          S, E, I, R,
                                          alpha, m.alpha, sigma.alpha,
                                          beta, m.beta, sigma.beta,
                                          gamma, m.gamma, sigma.gamma,
                                          kappa, a.kappa, b.kappa,
                                          lambda, a.lambda, b.lambda,
                                          p1, a.p1, b.p1, m.p1, sigma.p1,
                                          p2, a.p2, b.p2, m.p2, sigma.p2,
                                          new.Px, delta.mat,
                                          f, a.f, b.f,
                                          pop.size)-
                     log(dtruncnorm(new.pi.k.11, mean = pi.k[1], sd = step.size))-
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
      new.pi.k = new.pi.k       # accept move with probabily min(1,r)
      indicator = 1    # indicator of acceptance
    } else{
      new.pi.k = pi.k        # otherwise "reject" move, and stay where we are
      indicator = 0
    }

    new.Px[k,] <- new.pi.k

  }else if(k==2){
    new.pi.k.22 <- rtruncnorm(1, a=0, b=1, mean=pi.k[2], sd=step.size)
    new.pi.k.21 <- 1-new.pi.k.22
    new.pi.k <- c(new.pi.k.21,  new.pi.k.22)

    new.Px <- Px
    new.Px[k,] <- new.pi.k                  # Replace the k-th row with new pi

    # compute log of the acceptance ratio
    log.r <- min(c(log(dtruncnorm(pi.k[2], mean = new.pi.k.22, sd = step.size)) +
                     log.full.conditional(y, x,             # y_1:T x_0:T
                                          S, E, I, R,
                                          alpha, m.alpha, sigma.alpha,
                                          beta, m.beta, sigma.beta,
                                          gamma, m.gamma, sigma.gamma,
                                          kappa, a.kappa, b.kappa,
                                          lambda, a.lambda, b.lambda,
                                          p1, a.p1, b.p1, m.p1, sigma.p1,
                                          p2, a.p2, b.p2, m.p2, sigma.p2,
                                          new.Px, delta.mat,
                                          f, a.f, b.f,
                                          pop.size)-
                     log(dtruncnorm(new.pi.k.22, mean = pi.k[2], sd = step.size)) -
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
      new.pi.k = new.pi.k       # accept move with probabily min(1,r)
      indicator = 1    # indicator of acceptance
    } else{
      new.pi.k = pi.k        # otherwise "reject" move, and stay where we are
      indicator = 0
    }

    new.Px[k,] <- new.pi.k
  }


   return(list("indicator" = indicator,
              "k" = k,
              "newPx" = new.Px))
}

