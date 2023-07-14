
# Log Likelihood for BDSSSM-SEIR
source("Update_SEIR_RK4.R")
source("ObservationDensity.R")
library("HiddenMarkov")

# For two regimes
log.full.conditional <- function(y, x, 
                       S, E, I, R, 
                       alpha, m.alpha, sigma.alpha,
                       beta, m.beta, sigma.beta,
                       gamma, m.gamma, sigma.gamma,
                       kappa, a.kappa, b.kappa, 
                       lambda, a.lambda, b.lambda,
                       p, a.p, b.p, m.p, sigma.p,
                       Px, delta.mat,
                       f, a.f, b.f, 
                       pop.size){

  # mu(theta_1|x_1)mu(x1)π(α)π(β)π(γ)π(p)π(λ)π(κ)π(PX )π(f)
  log.targ  <- log(dtruncnorm(alpha, a=0, b=Inf, mean = m.alpha, sd = sigma.alpha)) +
    log(dtruncnorm(beta, a=0, b=Inf, mean = m.beta, sd = sigma.beta)) + 
    log(dtruncnorm(gamma, a=0, b=Inf, mean = m.gamma, sd = sigma.gamma)) + 
    dgamma(kappa, shape = a.kappa, rate = b.kappa, log=TRUE) +
    dgamma(lambda, shape = a.lambda , rate = b.lambda, log=TRUE) +
    log(dtruncnorm(p, a=a.p, b=b.p, mean = m.p, sd = sigma.p)) + 
    sum(DirichletReg::ddirichlet(Px, alpha=delta.mat,log=TRUE)) +
    sum(dunif(f, min=a.f, max=b.f, log=TRUE)) +
    DirichletReg::ddirichlet(matrix(c(S[1], E[1], I[1], R[1]), nrow=1, ncol=4), alpha=c(100,1,1,1),log=TRUE)+
    log(1/nrow(f))
  
  
  # gψ (xt|xt−1) for t=2,...,T
  T <- length(y)
  df <- data.frame("t0" = x[1:T-1],
                   "t1" = x[2:T],
                   "trans.prob" = rep(0, T-1))
  for (i in 1:nrow(df)){
    df[i, "trans.prob"] <- Px[df[i, "t0"], df[i, "t1"]]
  }
  log.targ <- log.targ + sum(log(df$trans.prob))
  
  # hψ(yt|θt, xt) from t=1,...T; gψ(θt|θt−1, xt) from t=2,...,T
  log.targ <- log.targ + log.obs.density(y[1], I[1], lambda, p)
  runge.kutta <- matrix(NA, nrow=4, ncol=T-1)
  log.dirichlet.density <- matrix(NA, nrow=1, ncol=T-1)
  for (t in 2:T){
    
    runge.kutta[,t-1] <- unlist(update.SEIR.rk4(S[t-1], E[t-1], I[t-1], R[t-1],
                                   alpha, beta, gamma,
                                   f[x[t]], pop.size))
    log.dirichlet.density[1,t-1] <- DirichletReg::ddirichlet(matrix(c(S[t], E[t], I[t], R[t]), nrow = 1),
                                                  alpha = kappa*runge.kutta[,t-1],
                                                  log = TRUE) 
    log.targ <- log.targ + log.obs.density(y[t], I[t], lambda, p) + log.dirichlet.density[1,t-1]
  }
  
  return(log.targ)
  
}


