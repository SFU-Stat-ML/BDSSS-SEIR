
# Set the working directory to current path 
library("rstudioapi")
cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

# Load R scripts for well-defined functions
source("SMC - BDSSSM-SEIR.R")
source("CSMC-AS-Rep-M_BDSSSM-SEIR.R") # conditional SMC with replicator M
source("Update_alpha.R")
source("Update_beta.R")
source("Update_gamma.R")
source("Update_lambda.R")
source("Update_kappa.R")
source("Update_p.R")
source("Update_Px.R")
source("Update_f.R")
source("InitialFun.R")  # Define initial values of latent variables here
source("ReplaceZeroInTrajectory.R")


############################ PG sampler function #############################


PG.CSMC.AS <- function(y, regimes, M, niter, hyperparams, pop.size=1){
  
  T <- length(y)-1
  lenXset <- length(regimes)
  N <- M*lenXset 


  # Set up matrices to store particles
  SsampleMat.CSMC.AS.repM <- matrix(0, niter, T + 1, byrow = TRUE)
  EsampleMat.CSMC.AS.repM <- matrix(0, niter, T + 1, byrow = TRUE)
  IsampleMat.CSMC.AS.repM <- matrix(0, niter, T + 1, byrow = TRUE)
  RsampleMat.CSMC.AS.repM <- matrix(0, niter, T + 1, byrow = TRUE)
  XsampleMat.CSMC.AS.repM <- matrix(0, niter, T + 1, byrow = TRUE)

  parameters.CSMC.AS.repM <-
    list(
      "alpha" = matrix(0, 1, niter, byrow = TRUE),
      "beta" = matrix(0, 1, niter, byrow = TRUE),
      "gamma" = matrix(0, 1, niter, byrow = TRUE),
      "lambda" = matrix(0, 1, niter, byrow = TRUE),
      "kappa" = matrix(0, 1, niter, byrow = TRUE),
      "Px" = vector(mode = "list", length = niter),
      "f" = matrix(0, lenXset, niter, byrow = TRUE),
      "p" = matrix(0, 1, niter, byrow = TRUE)
    )
  
  marginalLogLik.CSMC.AS.repM <- matrix(0, 1, niter, byrow = TRUE)
  
  
  # Prior distributions for ψ = {α, β, γ, λ, κ, PX , fxt }
  
  ### Hyperparameters setting
  m.alpha <- hyperparams$m.alpha
  sigma.alpha <- hyperparams$sigma.alpha
  m.beta <- hyperparams$m.beta
  sigma.beta <- hyperparams$sigma.beta
  m.gamma <- hyperparams$m.gamma
  sigma.gamma <- hyperparams$sigma.gamma
  a.lambda <- hyperparams$a.lambda
  b.lambda <- hyperparams$b.lambda
  a.kappa <- hyperparams$a.kappa
  b.kappa <- hyperparams$b.kappa
  a.p <- hyperparams$a.p
  b.p <- hyperparams$b.p
  m.p <- hyperparams$m.p
  sigma.p <- hyperparams$sigma.p
  delta.mat <- hyperparams$delta.mat
  a.f <- hyperparams$a.f
  b.f <- hyperparams$b.f


  # (1) Initialisation, r=0
  # Choose psi arbitrarily, and draw {theta^(b_0:T)_0:T,x^(b_0:T)_0:T} from {Theta_0:T,X_0:T,A_1:T}
  r = 1
  print(paste('PG-CSMC-AS-repM: iteration ', r))
  
  ## Arbitrarily generate psi from prior distributions
  ### SEIR model parameters: alpha, beta, gamma
  (parameters.CSMC.AS.repM$alpha[1, r] <- rtruncnorm(1, a = 0, b = Inf, mean = m.alpha, sd = sigma.alpha))
  (parameters.CSMC.AS.repM$beta[1, r] <- rtruncnorm(1, a = 0, b = Inf, mean = m.beta, sd = sigma.beta))
  (parameters.CSMC.AS.repM$gamma[1, r] <- rtruncnorm(1, a = 0, b = Inf, mean = m.gamma, sd = sigma.gamma))
  
  
  ### Precision parameters: lambda, kappa
  (parameters.CSMC.AS.repM$lambda[1, r] <-
      rgamma(1, shape = a.lambda, rate = b.lambda))  # mean = a.lambda/b.lambda
  (parameters.CSMC.AS.repM$kappa[1, r] <-
      rgamma(1, shape = a.kappa, rate = b.kappa))

  
  ### Transition probability matrix: Px
  library("gtools")
  (parameters.CSMC.AS.repM$Px[[r]] <-
      matrix(
        t(DirichletReg::rdirichlet(lenXset, delta.mat)),
        nrow = lenXset,
        ncol = lenXset,
        byrow = TRUE
      ))
  
  
  ### Transmission Rate modifier: f_{x_t}
  (parameters.CSMC.AS.repM$f[, r] <- c(1,
                                       runif(1, min = a.f[1], max = b.f[1]),
                                       runif(1, min = a.f[2], max = b.f[2]),
                                       runif(1, min = a.f[3], max = b.f[3])))
  
  
  ### Detection Rate: p
  (parameters.CSMC.AS.repM$p[1, r] <- rtruncnorm(1, a=a.p, b=b.p, mean=m.p, sd=sigma.p))
  

  ### Draw {theta^(b_0:T)_0:T,x^(b_0:T)_0:T} from {THETA_0:T,X_0:T,A_1:T} by running one iteration of SMC
  ptm <- proc.time()
  SMC.results <- SMC(
    y,
    regimes,
    N,
    parameters.CSMC.AS.repM$Px[[r]],
    parameters.CSMC.AS.repM$alpha[1, r],
    parameters.CSMC.AS.repM$beta[1, r],
    parameters.CSMC.AS.repM$gamma[1, r],
    parameters.CSMC.AS.repM$lambda[1, r],
    parameters.CSMC.AS.repM$kappa[1, r],
    parameters.CSMC.AS.repM$p[1, r],
    parameters.CSMC.AS.repM$f[, r],
    pop.size
  )
  proc.time() - ptm
  
  
  marginalLogLik.CSMC.AS.repM[r] <- sum(SMC.results$logLikY)
  
  ### Sample the reference trajectory
  # Draw L in {1:N}
  L <- sample(1:N, 1, prob = SMC.results$normalisedWeights[, T + 1])
  
  # Create reference trajectory
  RefParticleS <- rep(0, T + 1)
  RefParticleE <- rep(0, T + 1)
  RefParticleI <- rep(0, T + 1)
  RefParticleR <- rep(0, T + 1)
  RefParticleX <- rep(0, T + 1)
  
  # Recover b_0:T
  b <- rep(0, T + 1)
  b[T + 1] <- L
  
  RefParticleS[T + 1] <- SMC.results$particlesS[b[T + 1], T + 1]
  RefParticleE[T + 1] <- SMC.results$particlesE[b[T + 1], T + 1]
  RefParticleI[T + 1] <- SMC.results$particlesI[b[T + 1], T + 1]
  RefParticleR[T + 1] <- SMC.results$particlesR[b[T + 1], T + 1]
  RefParticleX[T + 1] <- SMC.results$particlesX[b[T + 1], T + 1]
  
  for (t in T:1) {
    b[t] <- SMC.results$AncestorLineage[b[t + 1], t]
    
    RefParticleS[t] <- SMC.results$particlesS[b[t], t]
    RefParticleE[t] <- SMC.results$particlesE[b[t], t]
    RefParticleI[t] <- SMC.results$particlesI[b[t], t]
    RefParticleR[t] <- SMC.results$particlesR[b[t], t]
    RefParticleX[t] <- SMC.results$particlesX[b[t], t]
    
  }
  
  SsampleMat.CSMC.AS.repM[r, ] <- replace.zero(RefParticleS)
  EsampleMat.CSMC.AS.repM[r, ] <- replace.zero(RefParticleE)
  IsampleMat.CSMC.AS.repM[r, ] <- replace.zero(RefParticleI)
  RsampleMat.CSMC.AS.repM[r, ] <- replace.zero(RefParticleR)
  XsampleMat.CSMC.AS.repM[r, ] <- replace.zero(RefParticleX)
  
  
  
  # Compute acceptance rate
  accept.kappa <- c()
  accept.lambda <- c()
  accept.alpha <- c()
  accept.beta <- c()
  accept.gamma <- c()
  accept.p <- c()
  accept.pi.1 <- c()
  accept.pi.2 <- c()
  accept.pi.3 <- c()
  accept.pi.4 <- c()
  accept.f <- matrix(NA, nrow = lenXset, ncol = niter - 1)
  
  
  # (2) For iteration at r = 1, ..., R
  ptm <- proc.time()
  
  pb <- txtProgressBar(min = 0, max = niter + 1, style = 3) # Show progress bar
  while (r < niter) {
    
  r = r+1
  
  # i) Run CSMC-AS with replicator M conditional on {theta(b_0:T)_0:T,x(b_0:T)_0:T} to obtain {Θ^1:N_0:T , X^1:N_0:T , A^1:N_1:T }
  CSMC.AS.repMresults <- CSMC.AS.repM(y, 
                                      regimes, 
                                      RefParticleS, RefParticleE, RefParticleI, RefParticleR, 
                                      RefParticleX,
                                      M, 
                                      parameters.CSMC.AS.repM$Px[[r-1]], 
                                      parameters.CSMC.AS.repM$alpha[1,r-1], 
                                      parameters.CSMC.AS.repM$beta[1,r-1], 
                                      parameters.CSMC.AS.repM$gamma[1,r-1], 
                                      parameters.CSMC.AS.repM$lambda[1,r-1], 
                                      parameters.CSMC.AS.repM$kappa[1,r-1],
                                      parameters.CSMC.AS.repM$p[1,r-1], 
                                      parameters.CSMC.AS.repM$f[,r-1],
                                      pop.size)                               # f_{x_t}
  
  
  # Get marginal likelihood
  marginalLogLik.CSMC.AS.repM[r] <- sum(CSMC.AS.repMresults$logLikY)
  
  ### ii) Sample the reference trajectory
  ### Sample the reference trajectory
  # Draw L in {1:N}
  L <- sample(1:N, 1, prob = CSMC.AS.repMresults$normalisedWeights[, T+1])
  
  # Create reference trajectory
  RefParticleS <- rep(0, T+1)
  RefParticleE <- rep(0, T+1)
  RefParticleI <- rep(0, T+1)
  RefParticleR <- rep(0, T+1)
  RefParticleX <- rep(0, T+1)
  
  # Recover b_0:T
  b <- rep(0, T+1) 
  b[T+1] <- L
  
  RefParticleS[T+1] <- CSMC.AS.repMresults$particlesS[b[T+1],T+1]
  RefParticleE[T+1] <- CSMC.AS.repMresults$particlesE[b[T+1],T+1]
  RefParticleI[T+1] <- CSMC.AS.repMresults$particlesI[b[T+1],T+1]
  RefParticleR[T+1] <- CSMC.AS.repMresults$particlesR[b[T+1],T+1]
  RefParticleX[T+1] <- CSMC.AS.repMresults$particlesX[b[T+1],T+1]
  
  for (t in T:1){
    b[t] <- CSMC.AS.repMresults$AncestorLineage[b[t+1], t]
    
    RefParticleS[t] <- CSMC.AS.repMresults$particlesS[b[t],t]
    RefParticleE[t] <- CSMC.AS.repMresults$particlesE[b[t],t]
    RefParticleI[t] <- CSMC.AS.repMresults$particlesI[b[t],t]
    RefParticleR[t] <- CSMC.AS.repMresults$particlesR[b[t],t]
    RefParticleX[t] <- CSMC.AS.repMresults$particlesX[b[t],t]
    
  }
  
  SsampleMat.CSMC.AS.repM[r,] <- replace.zero(RefParticleS)
  EsampleMat.CSMC.AS.repM[r,] <- replace.zero(RefParticleE)
  IsampleMat.CSMC.AS.repM[r,] <- replace.zero(RefParticleI)
  RsampleMat.CSMC.AS.repM[r,] <- replace.zero(RefParticleR) 
  XsampleMat.CSMC.AS.repM[r,] <- replace.zero(RefParticleX) 
  
  
  
  # iii) Draw \psi from p(\psi|\theta_(b_{0:T}), x_(b_{0:T}), b_0:T)
  
  # Latency rate: alpha
  mh.alpha.update <- update.alpha(y,  
                                  RefParticleX,
                                  RefParticleS, RefParticleE, RefParticleI, RefParticleR,
                                  parameters.CSMC.AS.repM$alpha[1,r-1], m.alpha, sigma.alpha,
                                  parameters.CSMC.AS.repM$beta[1,r-1], m.beta, sigma.beta,
                                  parameters.CSMC.AS.repM$gamma[1, r-1], m.gamma, sigma.gamma,
                                  parameters.CSMC.AS.repM$kappa[1, r-1], a.kappa, b.kappa,
                                  parameters.CSMC.AS.repM$lambda[1, r-1], a.lambda, b.lambda,
                                  parameters.CSMC.AS.repM$p[1, r-1], a.p, b.p,
                                  parameters.CSMC.AS.repM$Px[[r-1]], delta.mat,
                                  parameters.CSMC.AS.repM$f[, r-1], a.f, b.f,
                                  pop.size,
                                  0.02)
  parameters.CSMC.AS.repM$alpha[1,r] <- mh.alpha.update$new.alpha
  accept.alpha <- c(accept.alpha, mh.alpha.update$indicator)
  
  # Transmission rate: beta
  mh.beta.update <- update.beta(y,
                                RefParticleX,
                                RefParticleS, RefParticleE, RefParticleI, RefParticleR,
                                parameters.CSMC.AS.repM$alpha[1, r], m.alpha, sigma.alpha,
                                parameters.CSMC.AS.repM$beta[1, r -1], m.beta, sigma.beta,
                                parameters.CSMC.AS.repM$gamma[1, r -1], m.gamma, sigma.gamma,
                                parameters.CSMC.AS.repM$kappa[1, r -1], a.kappa, b.kappa,
                                parameters.CSMC.AS.repM$lambda[1, r -1], a.lambda, b.lambda,
                                parameters.CSMC.AS.repM$p[1, r-1], a.p, b.p,
                                parameters.CSMC.AS.repM$Px[[r - 1]], delta.mat,
                                parameters.CSMC.AS.repM$f[, r-1],a.f, b.f,
                                pop.size,
                                0.02)
  parameters.CSMC.AS.repM$beta[1,r] <- mh.beta.update$new.beta
  accept.beta <- c(accept.beta, mh.beta.update$indicator)
  

  
  
  # # Recovery rate: gamma
  mh.gamma.update<- update.gamma(y,
                                 RefParticleX,
                                 RefParticleS, RefParticleE, RefParticleI, RefParticleR,
                                 parameters.CSMC.AS.repM$alpha[1, r], m.alpha, sigma.alpha,
                                 parameters.CSMC.AS.repM$beta[1, r], m.beta, sigma.beta,
                                 parameters.CSMC.AS.repM$gamma[1, r -1], m.gamma, sigma.gamma,
                                 parameters.CSMC.AS.repM$kappa[1, r -1], a.kappa, b.kappa,
                                 parameters.CSMC.AS.repM$lambda[1, r -1], a.lambda, b.lambda,
                                 parameters.CSMC.AS.repM$p[1, r -1], a.p, b.p,
                                 parameters.CSMC.AS.repM$Px[[r - 1]], delta.mat,
                                 parameters.CSMC.AS.repM$f[, r -1], a.f, b.f,
                                 pop.size,
                                 0.02)
  parameters.CSMC.AS.repM$gamma[1,r] <- mh.gamma.update$new.gamma
  accept.gamma <- c(accept.gamma, mh.gamma.update$indicator)
  

  # # # Precision parameter: kappa
  mh.kappa.update <- update.kappa(y,
                                  RefParticleX,
                                  RefParticleS, RefParticleE, RefParticleI, RefParticleR,
                                  parameters.CSMC.AS.repM$alpha[1, r], m.alpha, sigma.alpha,
                                  parameters.CSMC.AS.repM$beta[1, r], m.beta, sigma.beta,
                                  parameters.CSMC.AS.repM$gamma[1, r], m.gamma, sigma.gamma,
                                  parameters.CSMC.AS.repM$kappa[1, r -1], a.kappa, b.kappa,
                                  parameters.CSMC.AS.repM$lambda[1, r-1], a.lambda, b.lambda,
                                  parameters.CSMC.AS.repM$p[1, r - 1], a.p, b.p,
                                  parameters.CSMC.AS.repM$Px[[r - 1]], delta.mat,
                                  parameters.CSMC.AS.repM$f[, r - 1],a.f, b.f,
                                  pop.size,
                                  500)
  parameters.CSMC.AS.repM$kappa[1,r] <- mh.kappa.update$new.kappa
  accept.kappa <- c(accept.kappa, mh.kappa.update$indicator)

  # # Precision parameter: lambda
  mh.lambda.update <- update.lambda(y,
                                    RefParticleX,
                                    RefParticleS, RefParticleE, RefParticleI, RefParticleR,
                                    parameters.CSMC.AS.repM$alpha[1, r], m.alpha, sigma.alpha,
                                    parameters.CSMC.AS.repM$beta[1, r], m.beta, sigma.beta,
                                    parameters.CSMC.AS.repM$gamma[1, r], m.gamma, sigma.gamma,
                                    parameters.CSMC.AS.repM$kappa[1, r], a.kappa, b.kappa,
                                    parameters.CSMC.AS.repM$lambda[1, r -1], a.lambda, b.lambda,
                                    parameters.CSMC.AS.repM$p[1, r -1], a.p, b.p,
                                    parameters.CSMC.AS.repM$Px[[r - 1]], delta.mat,
                                    parameters.CSMC.AS.repM$f[, r -1], a.f, b.f,
                                    pop.size,
                                    100)
  parameters.CSMC.AS.repM$lambda[1,r] <- mh.lambda.update$new.lambda
  accept.lambda <- c(accept.lambda, mh.lambda.update$indicator)

  
  # p
  mh.p.update <- update.p(y,
                          RefParticleX,
                          RefParticleS, RefParticleE, RefParticleI, RefParticleR,
                          parameters.CSMC.AS.repM$alpha[1, r], m.alpha, sigma.alpha,
                          parameters.CSMC.AS.repM$beta[1, r], m.beta, sigma.beta,
                          parameters.CSMC.AS.repM$gamma[1, r], m.gamma, sigma.gamma,
                          parameters.CSMC.AS.repM$kappa[1, r], a.kappa, b.kappa,
                          parameters.CSMC.AS.repM$lambda[1, r], a.lambda, b.lambda,
                          parameters.CSMC.AS.repM$p[1, r - 1], a.p, b.p,
                          parameters.CSMC.AS.repM$Px[[r - 1]], delta.mat,
                          parameters.CSMC.AS.repM$f[, r - 1], a.f, b.f,
                          pop.size,
                          0.02)
  parameters.CSMC.AS.repM$p[1, r] <- mh.p.update$new.p
  accept.p <- c(accept.p, mh.p.update$indicator)

  # pi.k in Px
  mh.pi.k.update <- update.pi.k(y,
                                RefParticleX,
                                RefParticleS, RefParticleE, RefParticleI, RefParticleR,
                                parameters.CSMC.AS.repM$alpha[1, r], m.alpha, sigma.alpha,
                                parameters.CSMC.AS.repM$beta[1, r], m.beta, sigma.beta,
                                parameters.CSMC.AS.repM$gamma[1, r], m.gamma, sigma.gamma,
                                parameters.CSMC.AS.repM$kappa[1, r], a.kappa, b.kappa,
                                parameters.CSMC.AS.repM$lambda[1, r], a.lambda, b.lambda,
                                parameters.CSMC.AS.repM$p[1, r], a.p, b.p,
                                parameters.CSMC.AS.repM$Px[[r - 1]], delta.mat,
                                parameters.CSMC.AS.repM$f[, r - 1], a.f, b.f,
                                pop.size,
                                0.02)
  parameters.CSMC.AS.repM$Px[[r]]<- mh.pi.k.update$newPx
  if (mh.pi.k.update$k == 1){
    accept.pi.1 <- c(accept.pi.1, mh.pi.k.update$indicator)
  }else if (mh.pi.k.update$k == 2){
    accept.pi.2 <- c(accept.pi.2, mh.pi.k.update$indicator)
  }else if (mh.pi.k.update$k == 3){
    accept.pi.3 <- c(accept.pi.3, mh.pi.k.update$indicator)
  }else if (mh.pi.k.update$k == 4){
    accept.pi.4 <- c(accept.pi.4, mh.pi.k.update$indicator)
  }
  
  
  # f_{x_t}: f1=1, update f2, f3 and f4
  parameters.CSMC.AS.repM$f[1,r] <- 1
  
  mh.f2.update  <- update.f2(y, RefParticleX,
                           RefParticleS, RefParticleE, RefParticleI, RefParticleR,
                           parameters.CSMC.AS.repM$alpha[1, r], m.alpha, sigma.alpha,
                           parameters.CSMC.AS.repM$beta[1, r], m.beta, sigma.beta,
                           parameters.CSMC.AS.repM$gamma[1, r], m.gamma, sigma.gamma,
                           parameters.CSMC.AS.repM$kappa[1, r], a.kappa, b.kappa,
                           parameters.CSMC.AS.repM$lambda[1, r], a.lambda, b.lambda,
                           parameters.CSMC.AS.repM$p[1, r], a.p, b.p,
                           parameters.CSMC.AS.repM$Px[[r]], delta.mat,
                           parameters.CSMC.AS.repM$f[, r - 1], a.f, b.f,
                           pop.size,
                           0.05)
  parameters.CSMC.AS.repM$f[2,r] <- mh.f2.update$new.f2
  accept.f[2,r-1] <-  mh.f2.update$indicator
  
  mh.f3.update  <- update.f3(y, RefParticleX,
                             RefParticleS, RefParticleE, RefParticleI, RefParticleR,
                             parameters.CSMC.AS.repM$alpha[1, r], m.alpha, sigma.alpha,
                             parameters.CSMC.AS.repM$beta[1, r], m.beta, sigma.beta,
                             parameters.CSMC.AS.repM$gamma[1, r], m.gamma, sigma.gamma,
                             parameters.CSMC.AS.repM$kappa[1, r], a.kappa, b.kappa,
                             parameters.CSMC.AS.repM$lambda[1, r], a.lambda, b.lambda,
                             parameters.CSMC.AS.repM$p[1, r], a.p, b.p,
                             parameters.CSMC.AS.repM$Px[[r]], delta.mat,
                             c(parameters.CSMC.AS.repM$f[1:2, r], parameters.CSMC.AS.repM$f[3:4, r-1]), a.f, b.f,
                             pop.size,
                             0.05)
  parameters.CSMC.AS.repM$f[3,r] <- mh.f3.update$new.f3
  accept.f[3,r-1] <-  mh.f3.update$indicator
  
  mh.f4.update  <- update.f4(y, RefParticleX,
                             RefParticleS, RefParticleE, RefParticleI, RefParticleR,
                             parameters.CSMC.AS.repM$alpha[1, r], m.alpha, sigma.alpha,
                             parameters.CSMC.AS.repM$beta[1, r], m.beta, sigma.beta,
                             parameters.CSMC.AS.repM$gamma[1, r], m.gamma, sigma.gamma,
                             parameters.CSMC.AS.repM$kappa[1, r], a.kappa, b.kappa,
                             parameters.CSMC.AS.repM$lambda[1, r], a.lambda, b.lambda,
                             parameters.CSMC.AS.repM$p[1, r], a.p, b.p,
                             parameters.CSMC.AS.repM$Px[[r]], delta.mat,
                             c(parameters.CSMC.AS.repM$f[1:3, r], parameters.CSMC.AS.repM$f[4, r-1]), a.f, b.f,
                             pop.size,
                             0.05)
  parameters.CSMC.AS.repM$f[4,r] <- mh.f4.update$new.f4
  accept.f[4,r-1] <-  mh.f4.update$indicator
  
  
  setTxtProgressBar(pb, r)
  }
  
  close(pb)
  proc.time() - ptm 
  
  ## Acceptance rate
  accept.rate <- data.frame("Parameters"=c("alpha", "beta", "gamma", "lambda", "kappa", "p", 
                                           "pi.1", "pi.2", "pi.3", "pi.4",  "f2", "f3", "f4"),
                            "AcceptanceRate"=c(sum(accept.alpha)/length(accept.alpha), 
                                               sum(accept.beta)/length(accept.beta), 
                                               sum(accept.gamma)/length(accept.gamma), 
                                               sum(accept.lambda)/length(accept.lambda),
                                               sum(accept.kappa)/length(accept.kappa),
                                               sum(accept.p)/length(accept.p),
                                               sum(accept.pi.1)/length(accept.pi.1),
                                               sum(accept.pi.2)/length(accept.pi.2),
                                               sum(accept.pi.3)/length(accept.pi.3),
                                               sum(accept.pi.4)/length(accept.pi.4),
                                               sum(accept.f[2,])/length(accept.f[2,]),
                                               sum(accept.f[3,])/length(accept.f[3,]),
                                               sum(accept.f[4,])/length(accept.f[4,])))
  
  
  return(list("SsampleMat.CSMC.AS.repM"=SsampleMat.CSMC.AS.repM,
              "EsampleMat.CSMC.AS.repM"=EsampleMat.CSMC.AS.repM,
              "IsampleMat.CSMC.AS.repM"=IsampleMat.CSMC.AS.repM,
              "RsampleMat.CSMC.AS.repM"=RsampleMat.CSMC.AS.repM,
              "XsampleMat.CSMC.AS.repM"=XsampleMat.CSMC.AS.repM,
              "parameters.CSMC.AS.repM"=parameters.CSMC.AS.repM,
              "marginalLogLik.CSMC.AS.repM"=marginalLogLik.CSMC.AS.repM,
              "acceptance.rate" = accept.rate))

}
