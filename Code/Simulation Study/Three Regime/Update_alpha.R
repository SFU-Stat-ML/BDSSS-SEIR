

## Update the latency parameter alpha 
source("log_full_conditional.R")

update.alpha <- function(y, x,           
                         S, E, I, R, 
                         alpha, m.alpha, sigma.alpha,
                         beta, m.beta, sigma.beta,
                         gamma, m.gamma, sigma.gamma,
                         kappa, a.kappa, b.kappa,
                         lambda, a.lambda, b.lambda,
                         p, a.p, b.p, m.p, sigma.p,
                         Px, delta.mat,
                         f, a.f, b.f,
                         pop.size,
                         step.size){

    new.alpha <- rtruncnorm(1,a=0, b=Inf, mean=alpha, sd=step.size)


    log.r <- min(c(log(dtruncnorm(alpha,  a=0, b=Inf, mean = new.alpha, sd = step.size)) +
                     log.full.conditional(y, x,             
                                          S, E, I, R, 
                                          new.alpha, m.alpha, sigma.alpha,
                                          beta, m.beta, sigma.beta,
                                          gamma, m.gamma, sigma.gamma,
                                          kappa, a.kappa, b.kappa,
                                          lambda, a.lambda, b.lambda,
                                          p, a.p, b.p, m.p, sigma.p,
                                          Px, delta.mat,
                                          f, a.f, b.f,
                                          pop.size)-
                     log(dtruncnorm(new.alpha,  a=0, b=Inf, mean = alpha, sd = step.size)) -
                     log.full.conditional(y, x,             
                                          S, E, I, R, 
                                          alpha, m.alpha, sigma.alpha,
                                          beta, m.beta, sigma.beta,
                                          gamma, m.gamma, sigma.gamma,
                                          kappa, a.kappa, b.kappa,
                                          lambda, a.lambda, b.lambda,
                                          p, a.p, b.p, m.p, sigma.p,
                                          Px, delta.mat,
                                          f, a.f, b.f,
                                          pop.size), log(1)))

    if(log(runif(1)) < log.r){
      new.alpha = new.alpha       # accept move with probability min(1,r)
      indicator = 1               # indicator of acceptance
    } else{
      new.alpha = alpha           # otherwise "reject" move, and stay where we are
      indicator = 0
    }


  return(list("new.alpha" = new.alpha,
              "indicator" = indicator))

}

