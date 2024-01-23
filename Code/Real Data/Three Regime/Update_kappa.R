

## BDSSSM-SEIR
## Update the precision parameter kappa in PG-CSMC-AS

library("DirichletReg")

update.kappa <- function(y, x,             # y_1:T x_0:T
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

  (new.kappa  <- rtruncnorm(1, a=0, b=Inf, mean = kappa, sd = step.size))

   # Calculate acceptance probability
  log.r <- min(c(log(dtruncnorm(kappa, a=0, b=Inf, mean = new.kappa, sd = step.size)) +
                 log.full.conditional(y, x,             # y_1:T x_0:T
                                      S, E, I, R,
                                      alpha, m.alpha, sigma.alpha,
                                      beta, m.beta, sigma.beta,
                                      gamma, m.gamma, sigma.gamma,
                                      new.kappa, a.kappa, b.kappa,
                                      lambda, a.lambda, b.lambda,
                                      p1, a.p1, b.p1, m.p1, sigma.p1,
                                      p2, a.p2, b.p2, m.p2, sigma.p2,
                                      Px, delta.mat,
                                      f, a.f, b.f,
                                      pop.size) -
                 log(dtruncnorm(new.kappa, a=0, b=Inf, mean = kappa, sd = step.size)) -
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

  # # Generate new kappa using log-normal proposal distribution (mean and standard deviation of the distribution are on the log scale )
  # log.kappa <- log(kappa)
  # new.kappa <- rlnorm(1, meanlog=log.kappa, sdlog=step.size)
  # new.log.kappa <- log(new.kappa)
  # 
  # # Calculate acceptance probability
  # log.r <- min(c(dlnorm(kappa, meanlog = new.log.kappa, sdlog = step.size, log = TRUE) +
  #                  log.full.conditional(y, x,             # y_1:T x_0:T
  #                                       S, E, I, R,
  #                                       alpha, m.alpha, sigma.alpha,
  #                                       beta, m.beta, sigma.beta,
  #                                       gamma, m.gamma, sigma.gamma,
  #                                       new.log.kappa, a.kappa, b.kappa,
  #                                       lambda, a.lambda, b.lambda,
  #                                       p, a.p, b.p,
  #                                       Px, delta.mat,
  #                                       f, a.f, b.f,
  #                                       pop.size) -
  #                  dlnorm(new.kappa, meanlog = log.kappa, sdlog = step.size, log = TRUE) -
  #                  log.full.conditional(y, x,             # y_1:T x_0:T
  #                                       S, E, I, R,
  #                                       alpha, m.alpha, sigma.alpha,
  #                                       beta, m.beta, sigma.beta,
  #                                       gamma, m.gamma, sigma.gamma,
  #                                       log.kappa, a.kappa, b.kappa,
  #                                       lambda, a.lambda, b.lambda,
  #                                       p, a.p, b.p,
  #                                       Px, delta.mat,
  #                                       f, a.f, b.f,
  #                                       pop.size), log(1)))

  if(log(runif(1)) < log.r){
    new.kappa = new.kappa       # accept move with probability min(1,r)
    indicator = 1                     # indicator of acceptance
  } else{
    new.kappa = kappa          # otherwise "reject" move, and stay where we are
    indicator = 0
  }
  return(list("new.kappa" = new.kappa,
              "indicator" = indicator))

}


# update.kappa <- function(y, x,             # y_1:T x_0:T
#                          S, E, I, R, 
#                          etaS, etaE, etaI, etaR,
#                          alpha, m.alpha, sigma.alpha,
#                          beta, m.beta, sigma.beta,
#                          gamma, m.gamma, sigma.gamma,
#                          kappa, a.kappa, b.kappa,
#                          lambda, a.lambda, b.lambda,
#                          p, a.p, b.p,
#                          Px, delta.mat,
#                          transm.modifier,
#                          pop.size,
#                          step.size){
# 
#   (new.kappa  <- rtruncnorm(1, a=step.size, b=Inf, mean = kappa, sd = step.size))
#   
#   # (new.kappa  <- rgamma(1, shape=kappa^2/step.size^2, rate=kappa/step.size^2))
#   
#   ll <- dgamma(kappa, shape = a.kappa, rate = b.kappa, log=TRUE)
#   ll.new <- dgamma(new.kappa, shape = a.kappa, rate = b.kappa, log=TRUE)
#   
#   T <- length(y)
#   for (t in 1:(T-1)) {
# 
#     ll <- ll + ddirichlet(matrix(c(S[t+1], E[t+1], I[t+1], R[t+1]), nrow = 1),
#                         alpha = kappa*c(etaS[t], etaE[t], etaI[t], etaR[t]),
#                         log=TRUE)
# 
#     ll.new <- ll.new + ddirichlet(matrix(c(S[t+1], E[t+1], I[t+1], R[t+1]), nrow = 1),
#                                   alpha = new.kappa*c(etaS[t], etaE[t], etaI[t], etaR[t]),
#                                   log=TRUE)
#   }
# 
#   log.r <- min(c(log(dtruncnorm(kappa, a=step.size, b=Inf, mean = new.kappa, sd = step.size)) +
#                    ll.new -
#                    log(dtruncnorm(new.kappa, a=step.size, b=Inf, mean = kappa, sd = step.size)) -
#                    ll), log(1))
# 
#   # log.r <- min(c(dgamma(kappa, shape=new.kappa^2/step.size^2, rate=new.kappa/step.size^2, log=TRUE) +
#   #                  ll.new -
#   #                  dgamma(new.kappa, shape=kappa^2/step.size^2, rate=kappa/step.size^2, log=TRUE)  -
#   #                  ll , log(1)))
#  
#   u <- log(runif(1))
#   
#   if(u < log.r){
#     new.kappa = new.kappa       # accept move with probability min(1,r)
#     indicator = 1                     # indicator of acceptance
#   } else{
#     new.kappa = kappa          # otherwise "reject" move, and stay where we are
#     indicator = 0
#   }
#   
#   return(list("new.kappa" = new.kappa,
#               "indicator" = indicator))
# 
# }
# 

# 
# # Adaptive MH step
# update.kappa <- function(y, x,             # y_1:T x_0:T
#                          S, E, I, R,
#                          alpha, m.alpha, sigma.alpha,
#                          beta, m.beta, sigma.beta,
#                          gamma, m.gamma, sigma.gamma,
#                          kappa, a.kappa, b.kappa,
#                          lambda, a.lambda, b.lambda,
#                          p, a.p, b.p,
#                          Px, delta.mat,
#                          f, a.f, b.f,
#                          pop.size,
#                          step.size){
#   
# # Initialization
# kappa <- kappa  # Starting value of kappa
# niter <- 10000  # Number of iterations
# accept <- rep(0, niter)  # Vector to keep track of acceptance rates
# tune <- 5  # Initial tuning parameter
# 
# # Adaptive Metropolis algorithm
# for (i in 1:niter) {
#   
#   # Update kappa using a normal random walk with adaptive step size
#   new.kappa <- rtruncnorm(1, a=0, b=Inf, mean = kappa, sd = tune)
#   
#   # Calculate acceptance probability
#   log.r <- min(c(log(dtruncnorm(kappa, a=0, b=Inf, mean = new.kappa, sd = tune)) +
#                    log.full.conditional(y, x,             # y_1:T x_0:T
#                                         S, E, I, R,
#                                         alpha, m.alpha, sigma.alpha,
#                                         beta, m.beta, sigma.beta,
#                                         gamma, m.gamma, sigma.gamma,
#                                         new.kappa, a.kappa, b.kappa,
#                                         lambda, a.lambda, b.lambda,
#                                         p, a.p, b.p,
#                                         Px, delta.mat,
#                                         f, a.f, b.f,
#                                         pop.size) -
#                    log(dtruncnorm(new.kappa, a=0, b=Inf, mean = kappa, sd = tune)) -
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
#   # Accept or reject move based on acceptance probability
#   if (log(runif(1)) < log.r) {
#     kappa <- new.kappa
#     accept[i] <- 1  # Mark move as accepted
#   } else {
#     accept[i] <- 0  # Mark move as rejected
#   }
#   
#   # Update tuning parameter using adaptive Metropolis criterion
#   if (i %% 100 == 0) {
#     acc.rate <- mean(accept[(i-99):i])
#     if (acc.rate < 0.2) {
#       tune <- tune / 2
#     } else if (acc.rate > 0.5) {
#       tune <- tune * 2
#     }
#   }
#   
# }
# 
#    # Return final value of kappa and acceptance rate
#    list("kappa" = kappa, "acceptance_rate" = mean(accept))
# }



