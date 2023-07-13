
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
source("ReplaceZeroInTrajectory.R") 
source("ObservationDensity.R")

# --------------- Algorithm: CSMC-AS with Replicator M ------------------------]
# ----- Input: 
#         y - measurements, length T
#         regimes - values X could take from {1,...,K}
#         S_ref, E_ref, I_ref, R_ref - reference particle on S, E, I, R, length T
#         X_ref - reference particle on latent variable X (discrete), length T
#         M - number of particles of Xt produced for each value of Xt (replicator)
#         Px - transition probabilities
#         alpha, beta, gamma - SEIR parameters 
#         lambda, kappa - precision parameters
#         p - detection rate
#         f - a 1xK vector of transmission rate modifier in K regimes

# Conditional SMC
CSMC.AS.repM <- function(y,
                         regimes,
                         S_ref,
                         E_ref,
                         I_ref,
                         R_ref,
                         X_ref,
                         M,
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
  K <- length(regimes) # Number of states
  N <- K*M             # Number of total particles
  
  
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
  
  
  # Initialization, t=1
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
  
  
  # Replace m-th particle with the element in reference trajectory
  m <- X_ref[1]*M
  
  particlesS[m, 1] <- as.vector(S_ref)[1]
  particlesE[m, 1] <- as.vector(E_ref)[1]
  particlesI[m, 1] <- as.vector(I_ref)[1]
  particlesR[m, 1] <- as.vector(R_ref)[1]
  particlesX[m, 1] <- as.vector(X_ref)[1]
  
  
  # Compute initial importance weights as p_\psi(theta_1, x_1)/q(theta_1, x_1)
  logWeights <- log.obs.density(y[1], particlesI[,1], lambda, p)
  max_logWeight <- max(logWeights, na.rm = TRUE)
  Weights[,1] <-  exp(logWeights - max_logWeight)
  normalisedWeights[, 1] <- Weights[,1]/sum(Weights[,1])
  
  # Compute marginal likelihood
  logLikelihood[1,1] <- log(sum(Weights[,1])) + max_logWeight - log(N)

  
  # For t=2,...,T
  while (t <= T){
    
    t <- t+1
    
    # Draw M ancestor indices using weights as the probability, and replicate them K times to obtain N ancestor indices
    ind <- rep(sample(1:N, M, replace = TRUE, prob = normalisedWeights[, t-1]), K)
    
    # Propagation step, x_t is deterministically controlled
    particlesX[ , t] <- rep(regimes, each = M)
    
    # Generate theta_t according to q(θ(i)t|θ(a(i)t)0:t−1,x(a(i)t)0:t−1,x(i)t)
    for (i in 1:N){
      
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
    
    # Replace zero values generated from Dirichlet distribution
    for (i in 1:N){
      particlesS[i, ] <- replace.zero(particlesS[i,])
      particlesE[i, ] <- replace.zero(particlesE[i,])
      particlesI[i, ] <- replace.zero(particlesI[i,])
      particlesR[i, ] <- replace.zero(particlesR[i,])
    }
    
    # Set the m-th particle to be reference path
    m <- X_ref[t]*M
    
    particlesS[m, t] <- as.vector(S_ref)[t]
    particlesE[m, t] <- as.vector(E_ref)[t]
    particlesI[m, t] <- as.vector(I_ref)[t]
    particlesR[m, t] <- as.vector(R_ref)[t]
    particlesX[m, t] <- as.vector(X_ref)[t]
    
    # Ancestor sampling step: Draw b_{t-1} 

    ## f(y_t|theta_t^{(b_t)}, x_t^{(b_t)})
    # log.fy <- dbeta(y[t], shape1 = lambda*p*particlesI[,t], shape2 = lambda*(1-p*particlesI[,t]), log = TRUE)      
    log.fy <- dbeta(y[t], shape1 = lambda*p*I_ref[t], shape2 = lambda*(1-p*I_ref[t]), log = TRUE)      
    
    ## gψ(θt|θ(bt−1)t−1,x(bt)t)
    log.gtheta <- matrix(NA, nrow=N, ncol=1)
    for (i in 1:N){
      
      # runge.kutta <- update.SEIR.rk4(particlesS[i,t-1],
      #                                particlesE[i,t-1],
      #                                particlesI[i,t-1],
      #                                particlesR[i,t-1],
      #                                alpha, beta, gamma,
      #                                f[particlesX[i,t]],
      #                                pop.size)
      # 
      # log.gtheta[i,1] <- DirichletReg::ddirichlet(matrix(c(S_ref[t], E_ref[t], I_ref[t], R_ref[t]), nrow=1),
      #                                      kappa * c(runge.kutta$S.new, runge.kutta$E.new,
      #                                                runge.kutta$I.new, runge.kutta$R.new), log = TRUE)

      runge.kutta <- update.SEIR.rk4(S_ref[t-1],
                                     E_ref[t-1],
                                     I_ref[t-1],
                                     R_ref[t-1],
                                     alpha, beta, gamma,
                                     f[X_ref[t]],
                                     pop.size)

      log.gtheta[i,1] <- DirichletReg::ddirichlet(matrix(c(particlesS[i,t], particlesI[i,t],
                                                           particlesE[i,t], particlesR[i,t]), nrow=1),
                                                  kappa * c(runge.kutta$S.new, runge.kutta$E.new,
                                                            runge.kutta$I.new, runge.kutta$R.new), log = TRUE)
    }
    
    ## gψ(x(bt)t|x(bt−1)t−1)
    # log.gx <- log(Px[particlesX[,t-1], X_ref[t]])
    log.gx <- log(Px[X_ref[t-1], X_ref[t]])
    
    ## f(y_t|theta_t^{(b_t)}, x_t^{(b_t)}) x gψ(θt|θ(bt−1)t−1,x(bt)t) x gψ(x(bt)t|x(bt−1)t−1)x w^(bt−1)_(t-1)
    log_w_as <- log(normalisedWeights[, t-1]) + log.fy + log.gtheta + log.gx
    # log_w_as <- log.fy + log.gtheta + log.gx 
    # log_w_as <- log.gtheta + log.gx + log(normalisedWeights[, t-1])
    normalize_w_as <- normalize_weight(log_w_as)  # Save the normalized weights

    # Replace a_t^(m) with b_{t-1} where m=x_t^(b_t)M
    m <- X_ref[t]*M
    ind[m] <- which(runif(1) < cumsum(normalize_w_as))[1]
    # ind[m] <- sample(1:N, 1, prob = normalize_w_as)
    AL[ ,t-1] <- ind  # Store indices of ancestral particles at time t-1
  
    
    
    # Weighting step (avoid numerical instability)
    ## Get target transition probability
    target_transition_probs <- c()
    for (i in 1:N) {
      state_t_1 <- particlesX[ind[i],t-1]
      state_t <- particlesX[i,t]
      target_transition_probs <- c(target_transition_probs, 
                                   Px[state_t_1, state_t])
    }
    
    logWeights <- log.obs.density(y[t], particlesI[,t], lambda, p) + log(target_transition_probs)
    max_logWeight <- max(logWeights, na.rm = TRUE)
    Weights[,t] <- exp(logWeights - max_logWeight)
    normalisedWeights[, t] <- Weights[, t] /sum(Weights[, t])

    # Compute marginal log-likelihood log[pψ(y_t|y_{1:t−1})] = log[1/N ∑_{i=1}^N w^(i)]
    logLikelihood[1,t]  <- log(sum(Weights[,t])) + max_logWeight - log(N)
    
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


