
# Log Likelihood for BDSSSM-SEIR
source("Update_SEIR_RK4.R")
library("HiddenMarkov")

log.full.conditional <- function(y, x,             # y_1:T x_0:T
                                 S, E, I, R, 
                                 alpha, m.alpha, sigma.alpha,
                                 beta, m.beta, sigma.beta,
                                 gamma, m.gamma, sigma.gamma,
                                 kappa, a.kappa, b.kappa, 
                                 lambda, a.lambda, b.lambda,
                                 p, a.p, b.p,
                                 Px, delta.mat,
                                 f, a.f, b.f, 
                                 pop.size){
  
  # π(α)π(β)π(γ)π(p)π(λ)π(κ)π(PX )π(f)
  log.targ  <- log(dtruncnorm(alpha, a=0, b=Inf, mean = m.alpha, sd = sigma.alpha)) +
    log(dtruncnorm(beta, a=0, b=Inf, mean = m.beta, sd = sigma.beta)) + 
    log(dtruncnorm(gamma, a=0, b=Inf, mean = m.gamma, sd = sigma.gamma)) + 
    dgamma(kappa, shape = a.kappa, rate = b.kappa, log=TRUE) +
    dgamma(lambda, shape = a.lambda , rate = b.lambda, log=TRUE) +
    log(dtruncnorm(p, a=a.p, b=b.p, mean=m.p, sd=sigma.p)) +
    sum(DirichletReg::ddirichlet(Px, alpha=delta.mat,log=TRUE)) +
    dunif(f[2], min=a.f[1], max=b.f[1], log=TRUE) +
    dunif(f[3], min=a.f[2], max=b.f[2], log=TRUE) +
    dunif(f[4], min=a.f[3], max=b.f[3], log=TRUE) +
    DirichletReg::ddirichlet(matrix(c(S[1], E[1], I[1], R[1]), nrow=1, ncol=4), alpha=c(100,1,1,1),log=TRUE) +
    log(1/length(f))
  
  
  # gψ (xt|xt−1)
  T <- length(y)
  df <- data.frame("t0" = x[1:T-1],
                   "t1" = x[2:T],
                   "trans.prob" = rep(0, T-1))
  for (i in 1:nrow(df)){
    df[i, "trans.prob"] <- Px[df[i, "t0"], df[i, "t1"]]
  }
  log.targ <- log.targ + sum(log(df$trans.prob))
  
  # fψ (yt|θ0:t, x1:t)gψ (θt|θ0:t−1, x1:t)gψ (xt|xt−1)from t=1,...T1
  # t=1
  log.targ <- log.targ + dbeta(y[1], lambda*p*I[1], lambda*(1-p*I[1]), log=TRUE)
  
  for (t in 2:T){
    
    runge.kutta <- update.SEIR.rk4(S[t-1], E[t-1], I[t-1], R[t-1],
                                   alpha, beta, gamma,
                                   f[x[t]], pop.size)
    
    log.targ <- log.targ + dbeta(y[t], lambda*p*I[t], lambda*(1-p*I[t]), log=TRUE) +
      DirichletReg::ddirichlet(matrix(c(S[t], E[t], I[t], R[t]), nrow = 1),
                               alpha = kappa*c(runge.kutta$S.new, runge.kutta$E.new,
                                               runge.kutta$I.new, runge.kutta$R.new),
                               log = TRUE) 
    
  }
  
  return(log.targ)
  
}
