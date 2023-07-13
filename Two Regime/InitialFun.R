

# Initialization on theta1 = [S1, E1, I1, R1] and x1

# Note: All initial values must be > 0

# Dirichlet hyperparameter

initial.theta <- function(N){
  
  xi <- c(100, 1, 1, 1) 
  theta.init <- rdirichlet(N, xi)
  theta.init.density <- ddirichlet(theta.init, xi)
  
  return(list("theta.init" = theta.init,
              "theta.init.density" = theta.init.density))
}

initial.x <- function(regimes, N){
  x.init <- sample(1:length(regimes), N, replace = TRUE)
  x.init.density <- rep(1/length(regimes), N)
  return(list("x.init" = x.init,
              "x.init.density" = x.init.density))
}


# initial.fun <- function(N, y.init, p){
#   
#   # For theta0
#   S.init <- rep(0.99, N)
#   E.init <- rep(0.005, N)
#   
#   # Truncated beta function
#   rtruncbeta <- function(n, a, b, lower, upper) {
#     u <- runif(n, 0, 1)
#     x <- qbeta(u, a, b)
#     x <- pmin(pmax(x, lower), upper)
#     return(x)
#   }
#   
#   I.init <- rtruncbeta(N, 1, 1/y.init, lower=y.init, upper=1-S.init-E.init)
#   # I.init <- rbeta(N, 1, 1/y.init)/p
#   # I.init <- rep(0.01, N)
#   R.init <- 1-S.init-E.init-I.init
#   
#   
#   neg.R <- sum(R.init < 0)
#   while (neg.R > 0){
#     
#     print("There are negative inital values in R compartment.")
#     
#     I.init <- rtruncbeta(N, 1, 1/y.init, lower=y.init, upper=1-S.init-E.init)
#     R.init <- 1-S.init-E.init-I.init
#     
#     neg.R <- sum(R.init < 0)
#   }
#   
#   theta.init <- matrix(c(S.init,  E.init, I.init, R.init), nrow=N, ncol=4)
#   colnames(theta.init) <- c("S", "E", "I", "R")
#   
#   # For switching variable
#   x.init <- sample(1:2, N, replace = TRUE)
#   
#   # # Joint proposal density q(θ_1, x_1)
#   # ## q(θ1)
#   # truncated_beta_density <- function(x, a, b, lower, upper) {
#   #   dbeta(x, a, b) / (pbeta(upper, a, b) - pbeta(lower, a, b))
#   # }
#   # 
#   # theta.init.density <- truncated_beta_density(I.init,  
#   #                                          1, 1/y.init, 
#   #                                          lower=y.init, upper=1-S.init-E.init)
#   # ## q(x0)
#   # x.init.density <- rep(1/2, N)
#   # 
#   # init.proposal.density <- theta.init.density* x.init.density
#   
#   
#   return(list("x.init" = x.init, 
#               "theta.init"= theta.init  ))
#   
# }
