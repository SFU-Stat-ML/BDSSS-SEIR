
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
                       p1, a.p1, b.p1, m.p1, sigma.p1,
                       p2, a.p2, b.p2, m.p2, sigma.p2,
                       Px, delta.mat,
                       f, a.f, b.f, 
                       pop.size){
  # truncated_beta_density <- function(x, a, b, lower, upper) {
  #   dbeta(x, a, b) / (pbeta(upper, a, b) - pbeta(lower, a, b))
  # }
  # 
  # π(α)π(β)π(γ)π(λ)π(κ)π(PX )π(theta_0)π(f)
  log.targ  <- log(dtruncnorm(alpha, a=0, b=Inf, mean = m.alpha, sd = sigma.alpha)) +
    log(dtruncnorm(beta, a=0, b=Inf, mean = m.beta, sd = sigma.beta)) + 
    log(dtruncnorm(gamma, a=0, b=Inf, mean = m.gamma, sd = sigma.gamma)) + 
    dgamma(kappa, shape = a.kappa, rate = b.kappa, log=TRUE) +
    dgamma(lambda, shape = a.lambda , rate = b.lambda, log=TRUE) +
    log(dtruncnorm(p1, a=a.p1, b=b.p1, mean = m.p1, sd = sigma.p1)) + 
    log(dtruncnorm(p2, a=a.p2, b=b.p2, mean = m.p2, sd = sigma.p2)) + 
    sum(DirichletReg::ddirichlet(Px, alpha=delta.mat,log=TRUE)) +
    dunif(f[2], min=a.f[1], max=b.f[1], log=TRUE) +
    dunif(f[3], min=a.f[2], max=b.f[2], log=TRUE) +
    dunif(f[4], min=a.f[3], max=b.f[3], log=TRUE) +
    DirichletReg::ddirichlet(matrix(c(S[1], E[1], I[1], R[1]), nrow=1, ncol=4), alpha=c(100,1,1,1),log=TRUE)+
    log(1/length(f))
  
  
    # gψ (xt|xt−1) for t=2,...,T
    T <- length(y)
    df <- data.frame("t0" = x[1:T-1],
                   "t1" = x[2:T],
                   "trans.prob" = rep(0, T-1))
    for (i in 1:nrow(df)){
      df[i, "trans.prob"] <- Px[df[i, "t0"], df[i, "t1"]]
    }
    log.targ <- log.targ + sum(log(df$trans.prob))
  
    # fψ(yt|θt, xt)π(\tilde{p}_t) from t=1,...T; gψ(θt|θt−1, xt) from t=2,...,T
    log.targ <- log.targ + log.obs.density(y[1], I[1], lambda, p1)
    runge.kutta <- matrix(NA, nrow=4, ncol=T-1)
    log.dirichlet.density <- matrix(NA, nrow=1, ncol=T-1)
    for (t in 2:T){
      if (t <= T_star){
        runge.kutta[,t-1] <- unlist(update.SEIR.rk4(S[t-1], E[t-1], I[t-1], R[t-1],
                                                    alpha, beta, gamma,
                                                    f[x[t]], pop.size))
        log.dirichlet.density[1,t-1] <- DirichletReg::ddirichlet(matrix(c(S[t], E[t], I[t], R[t]), nrow = 1),
                                                                 alpha = kappa*runge.kutta[,t-1],
                                                                 log = TRUE) 
        log.targ <- log.targ + log.obs.density(y[t], I[t], lambda, p1) + log.dirichlet.density[1,t-1]
      }else{
        runge.kutta[,t-1] <- unlist(update.SEIR.rk4(S[t-1], E[t-1], I[t-1], R[t-1],
                                                    alpha, beta, gamma,
                                                    f[x[t]], pop.size))
        log.dirichlet.density[1,t-1] <- DirichletReg::ddirichlet(matrix(c(S[t], E[t], I[t], R[t]), nrow = 1),
                                                                 alpha = kappa*runge.kutta[,t-1],
                                                                 log = TRUE) 
        log.targ <- log.targ + log.obs.density(y[t], I[t], lambda, p2) + log.dirichlet.density[1,t-1]
      }
 
  }
  
  return(log.targ)
  
}


