##################### Switching state-space SEIR model #########################
#
# yt|θt,ψ∼Beta(λpIt,λ(1−pIt))
# θt|θt−1,xt,ψ∼Dir(κg(θt−1,α,β,γ,fxt))

# ------------ Call libraries -------------
library("DirichletReg") # For log density of Dirichlet distribution

# ------------- Load scripts --------------
source("normalize_weights.R") # Weights normalization
source("Update_SEIR_RK4.R")
source("InitialFun.R")
source("ObservationDensity.R")
# --------------- Algorithm: CSMC-AS with Replicator M ------------------------]
# ----- Input: 
#         y - measurements, length T
#         regimes - values X could take from {1,...,K}
#         N - number of particles 
#         Px - transition probabilities
#         alpha, beta, gamma - SEIR parameters 
#         lambda, kappa - precision parameters
#         p - detection rate
#         f - a 1xK vector of transmission modifier in K regimes
#         pop.size - population size

# SMC
SMC <- function(y, 
                regimes,
                N, 
                Px,
                alpha,
                beta,
                gamma,
                lambda,
                kappa,
                p,
                f,
                pop.size) {
  
  
  T <- length(y) - 1     
  
  
  # Initialize variables
  particlesS <- matrix(0, nrow = N, ncol = T+1)   # Nx(T+1)
  particlesE <- matrix(0, nrow = N, ncol = T+1)   # Nx(T+1)
  particlesI <- matrix(0, nrow = N, ncol = T+1)   # Nx(T+1)
  particlesR <- matrix(0, nrow = N, ncol = T+1)   # Nx(T+1)
  particlesX <- matrix(0, nrow = N, ncol = T+1)   # Nx(T+1)
  AL <- matrix(NA, nrow = N, ncol = T)  # for ancestral lineage
  Weights <- matrix(NA, nrow = N, ncol = T+1)
  normalisedWeights <- matrix(NA, nrow = N, ncol = T+1)
  logLikelihood <- matrix(NA, nrow = 1, ncol = T+1)
  
  # Initialization, t=0
  t <- 1
  
  # Initialize switching variable at t=1
  x1 <- initial.x(regimes, N)
  particlesX[, 1] <- x1$x.init
  
  # Initialize particle set for states S, E, I, R at t=1; 
  theta1 <- initial.theta(N)
  particlesS[, 1] <- theta1$theta.init[,1]
  particlesE[, 1] <- theta1$theta.init[,2]
  particlesI[, 1] <- theta1$theta.init[,3]
  particlesR[, 1] <- theta1$theta.init[,4]

  # Compute initial importance weights as p_\psi(theta_0, x_0)/q(theta_0, x_0)
  logWeights <- log.obs.density(y[1], particlesI[,1], lambda, p)
  max_logWeight <- max(logWeights, na.rm = TRUE)
  Weights[,1] <-  exp(logWeights - max_logWeight)
  normalisedWeights[, 1] <- Weights[,1]/sum(Weights[,1])
  
  # Compute marginal likelihood
  logLikelihood[1,1] <- log(sum(Weights[,1])) + max_logWeight - log(N)
  
  
  while (t <= T){
    
    t <- t+1
    
    # Resample N particles
    ind <- sample(1:N, N, replace = TRUE, prob = normalisedWeights[, t-1])
    
    # Propagation step
    
    for (i in 1:N){
      
      # Draw {x(i)t } from g(x(i)t |x(a(i))0:t−1)
      particlesX[i, t] <- which(rmultinom(1, size=1, Px[particlesX[ind[i], t-1],])==1)
      
      # Generate theta_t according to q(θ(i)t|θ(a(i)t)0:t−1,x(a(i)t)0:t−1,x(i)t)
      runge.kutta <- update.SEIR.rk4(particlesS[ind[i], t-1],
                                     particlesE[ind[i], t-1],
                                     particlesI[ind[i], t-1],
                                     particlesR[ind[i], t-1],
                                     alpha, beta, gamma,
                                     f[particlesX[i, t]],
                                     pop.size)
      
      thetaPred <- rdirichlet(1, kappa * c(runge.kutta$S.new, runge.kutta$E.new, 
                                           runge.kutta$I.new, runge.kutta$R.new))
      
      particlesS[i, t] <- thetaPred[1]
      particlesE[i, t] <- thetaPred[2]
      particlesI[i, t] <- thetaPred[3]
      particlesR[i, t] <- thetaPred[4]
    }
    
    # Store indices of ancestral particles at time t-1
    AL[ ,t-1] <- ind  
    
    # Weighting step (avoid numerical instability)
    
    ## Get target transition probability
    logWeights <- log.obs.density(y[t], particlesI[,t], lambda, p)
    max_logWeight <- max(logWeights, na.rm = TRUE)
    Weights[,t] <- exp(logWeights - max_logWeight)
    normalisedWeights[, t] <- Weights[, t] /sum(Weights[, t])
    
    # Compute marginal log-likelihood
    logLikelihood[1,t] <- log(sum(Weights[,t])) + max_logWeight - log(N)
    

  }
  
  results <- list("particlesS" = particlesS,
                  "particlesE" = particlesE,
                  "particlesI" = particlesI,
                  "particlesR" = particlesR,
                  "particlesX" = particlesX,
                  "Weights" = Weights,
                  "normalisedWeights" = normalisedWeights,
                  "AncestorLineage" = AL,
                  "logLikY" = logLikelihood)
  return(results)
}
