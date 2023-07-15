

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


