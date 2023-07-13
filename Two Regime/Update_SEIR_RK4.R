

##### Goal: Update the SEIR states at time t

update.SEIR.rk4 <- function(S.old, E.old, I.old, R.old, alpha, beta, gamma, f, pop.size){
  
  k.S1.old <- -f*beta/pop.size*S.old*I.old
  k.E1.old <- f*beta/pop.size*S.old*I.old - alpha*E.old
  k.I1.old <- alpha*E.old - gamma*I.old
  k.R1.old <- gamma*I.old
  
  k.S2.old <- -f*beta/pop.size*(S.old + 0.5*k.S1.old)*(I.old + 0.5*k.I1.old)
  k.E2.old <- f*beta/pop.size*(S.old + 0.5*k.S1.old)*(I.old + 0.5*k.I1.old)-alpha*(E.old + 0.5*k.E1.old)
  k.I2.old <- alpha*(E.old + 0.5*k.E1.old) - gamma*(I.old + 0.5*k.I1.old)
  k.R2.old <- gamma*(I.old + 0.5*k.I1.old)
  
  k.S3.old <- -f*beta/pop.size*(S.old + 0.5*k.S2.old)*(I.old + 0.5*k.I2.old)
  k.E3.old <- f*beta/pop.size*(S.old + 0.5*k.S2.old)*(I.old + 0.5*k.I2.old)-alpha*(E.old + 0.5*k.E2.old)
  k.I3.old <- alpha*(E.old + 0.5*k.E2.old) - gamma*(I.old + 0.5*k.I2.old)
  k.R3.old <- gamma*(I.old + 0.5*k.I2.old)
  
  k.S4.old <- -f*beta/pop.size*(S.old + 0.5*k.S3.old)*(I.old + 0.5*k.I3.old)
  k.E4.old <- f*beta/pop.size*(S.old + 0.5*k.S3.old)*(I.old + 0.5*k.I3.old)-alpha*(E.old + 0.5*k.E3.old)
  k.I4.old <- alpha*(E.old + 0.5*k.E3.old) - gamma*(I.old + 0.5*k.I3.old)
  k.R4.old <- gamma*(I.old + 0.5*k.I3.old)
  
  # S.new <- S.old + (k.S1.old + 2*k.S2.old + 2*k.S3.old + k.S4.old)/6
  # E.new <- E.old + (k.E1.old + 2*k.E2.old + 2*k.E3.old + k.E4.old)/6
  # I.new <- I.old + (k.I1.old + 2*k.I2.old + 2*k.I3.old + k.I4.old)/6
  # R.new <- R.old + (k.R1.old + 2*k.R2.old + 2*k.R3.old + k.R4.old)/6
  
  # Implement a non-negativity constraint
  S.new <- max(S.old + (k.S1.old + 2*k.S2.old + 2*k.S3.old + k.S4.old)/6, 1e-6)
  E.new <- max(E.old + (k.E1.old + 2*k.E2.old + 2*k.E3.old + k.E4.old)/6, 1e-6)
  I.new <- max(I.old + (k.I1.old + 2*k.I2.old + 2*k.I3.old + k.I4.old)/6, 1e-6)
  R.new <- max(R.old + (k.R1.old + 2*k.R2.old + 2*k.R3.old + k.R4.old)/6, 1e-6)
  
  return((list("S.new" = S.new,
               "E.new" = E.new,
               "I.new" = I.new,
               "R.new" = R.new
  )))
}



update.SEIR.rk5 <- function(S.old, E.old, I.old, R.old, alpha, beta, gamma, f, pop.size){
  
  k.S1.old <- -f*beta/pop.size*S.old*I.old
  k.E1.old <- f*beta/pop.size*S.old*I.old - alpha*E.old
  k.I1.old <- alpha*E.old - gamma*I.old
  k.R1.old <- gamma*I.old
  
  k.S2.old <- -f*beta/pop.size*(S.old + 0.25*k.S1.old)*(I.old + 0.25*k.I1.old)
  k.E2.old <- f*beta/pop.size*(S.old + 0.25*k.S1.old)*(I.old + 0.25*k.I1.old)-alpha*(E.old + 0.25*k.E1.old)
  k.I2.old <- alpha*(E.old + 0.25*k.E1.old) - gamma*(I.old + 0.25*k.I1.old)
  k.R2.old <- gamma*(I.old + 0.25*k.I1.old)
  
  k.S3.old <- -f*beta/pop.size*(S.old + (k.S1.old + k.S2.old)/8)*(I.old + (k.I1.old + k.I2.old)/8)
  k.E3.old <- f*beta/pop.size*(S.old + (k.S1.old + k.S2.old)/8)*(I.old + (k.I1.old + k.I2.old)/8)-alpha*(E.old + (k.E1.old + k.E2.old)/8)
  k.I3.old <- alpha*(E.old + (k.E1.old + k.E2.old)/8) - gamma*(I.old + (k.I1.old + k.I2.old)/8)
  k.R3.old <- gamma*(I.old + (k.I1.old + k.I2.old)/8)
  
  k.S4.old <- -f*beta/pop.size*(S.old + 0.5*k.S2.old + k.S3.old)*(I.old + 0.5*k.I2.old + k.I3.old)
  k.E4.old <- f*beta/pop.size*(S.old + 0.5*k.S2.old+ k.S3.old)*(I.old + 0.5*k.I2.old + k.I3.old)-alpha*(E.old + 0.5*k.E2.old + k.E3.old)
  k.I4.old <- alpha*(E.old + 0.5*k.E2.old + k.E3.old) - gamma*(I.old + 0.5*k.I2.old + k.I3.old)
  k.R4.old <- gamma*(I.old + 0.5*k.I2.old + k.I3.old)
  
  k.S5.old <- -f*beta/pop.size*(S.old + 3*k.S1.old/16 + 9*k.S4.old/16)*(I.old + 3*k.I1.old/16 + 9*k.I4.old/16)
  k.E5.old <- f*beta/pop.size*(S.old + 3*k.S1.old/16 + 9*k.S4.old/16)*(I.old + 3*k.I1.old/16 + 9*k.I4.old/16)-alpha*(E.old + 3*k.E1.old/16 + 9*k.E4.old/16)
  k.I5.old <- alpha*(E.old + 3*k.E1.old/16 + 9*k.E4.old/16) - gamma*(I.old + 3*k.I1.old/16 + 9*k.I4.old/16)
  k.R5.old <- gamma*(I.old + 3*k.I1.old/16 + 9*k.I4.old/16)
  
  k.S6.old <- -f*beta/pop.size*(S.old + 3*k.S1.old/7 + 2*k.S2.old/7 + 12*k.S3.old/7 - 12*k.S4.old/7 + 8*k.S5.old/7)*(I.old + 3*k.I1.old/7 + 2*k.I2.old/7 + 12*k.I3.old/7 - 12*k.I4.old/7 + 8*k.I5.old/7)
  k.E6.old <- f*beta/pop.size*(S.old + 3*k.S1.old/7 + 2*k.S2.old/7 + 12*k.S3.old/7 - 12*k.S4.old/7 + 8*k.S5.old/7)*(I.old + 3*k.I1.old/7 + 2*k.I2.old/7 + 12*k.I3.old/7 - 12*k.I4.old/7 + 8*k.I5.old/7)-alpha*(E.old + 3*k.E1.old/7 + 2*k.E2.old/7 + 12*k.E3.old/7 - 12*k.E4.old/7 + 8*k.E5.old/7)
  k.I6.old <- alpha*(E.old + 3*k.E1.old/7 + 2*k.E2.old/7 + 12*k.E3.old/7 - 12*k.E4.old/7 + 8*k.E5.old/7) - gamma*(I.old + 3*k.I1.old/7 + 2*k.I2.old/7 + 12*k.I3.old/7 - 12*k.I4.old/7 + 8*k.I5.old/7)
  k.R6.old <- gamma*(I.old + 3*k.I1.old/7 + 2*k.I2.old/7 + 12*k.I3.old/7 - 12*k.I4.old/7 + 8*k.I5.old/7)
  
  S.new <- S.old + (7*k.S1.old + 32*k.S3.old + 12*k.S4.old + 32*k.S5.old + 7*k.S6.old)/90
  E.new <- E.old + (7*k.E1.old + 32*k.E3.old + 12*k.E4.old + 32*k.E5.old + 7*k.E6.old)/90
  I.new <- I.old + (7*k.I1.old + 32*k.I3.old + 12*k.I4.old + 32*k.I5.old + 7*k.I6.old)/90
  R.new <- R.old + (7*k.R1.old + 32*k.R3.old + 12*k.R4.old + 32*k.R5.old + 7*k.R6.old)/90
  
  return((list("S.new" = S.new,
               "E.new" = E.new,
               "I.new" = I.new,
               "R.new" = R.new
  )))
}
