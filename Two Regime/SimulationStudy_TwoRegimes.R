
#------------------------ Clear memory and graphics ---------------------------#
rm(list=ls())
graphics.off()

options(digits = 16)

# ----------------Set the working directory to current path -------------------#
library("rstudioapi")
cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

# --------- I. Simulate data under switching state-space SEIR model -----------#

# yt|θt,ψ ∼ Beta(λpIt, λ(1 − pIt))
# θt|θt−1, ψ ∼ Dir(κg(θt−1; α, β, γ, fxt ))

#-------------------------------- Two Regimes ----------------------------------#
source("Update_SEIR_RK4.R") # Rk4 ODE solver

seed=1580
set.seed(seed)

# # Set up true model parameters  
pop.size = 1
alpha = 1/3               # the incubation rate, the average incubation period is 1/alpha (Use week as unit)
beta = 0.39                # The transmission rate
gamma = 0.18               # the recovery rate, time to recovery or death is 1/gamma (Use week as unit)
lambda = 2500             # parameter in Beta(λpIt,λ(1−pIt))
kappa = 5500          # parameter in Dir(κf_xt(θt−1,α,β,γ))
p = 0.25                   # parameter in Beta(λpIt,λ(1−pIt))


# Set up regimes
regimes <- c(1, 2)
lenXset <- length(regimes)

# Set up transition probability matrix, row sum=1
Px <- matrix(c(0.9, 0.1,
               0.1, 0.9),
             nrow = 2,
             ncol = 2,
             byrow = TRUE)


# Set up transmission rate modifier
f <- matrix(0, nrow=length(regimes), ncol=1)
f[1,1] <- 1                    # Transmission rate modifier under regime 1
f[2,1] <- 0.1

# Data record of length
T <- 149

# Set up matrices to store simulated data
theta <- matrix(0, nrow=4, ncol=T+1)
rownames(theta) <- c("S", "E", "I", "R")
y     <- matrix(0, nrow=1, ncol=T+1)
x     <- matrix(0, nrow=1, ncol=T+1)

# At t=0

# Starting values
theta[1,1] <- 0.99
theta[2,1] <- 0.001
theta[3,1] <- 0.003
theta[4,1] <- 1-theta[1,1]-theta[2,1]-theta[3,1]

# Simulate discrete Markov chain for switch variable Xt
x[1] <- 1

# Simulate y_t
y[1] <- rbeta(1, shape1 = lambda*p*theta[3,1], shape2 = lambda*(1-p*theta[3,1]))

# At time 1, 2, 3, ..., T
library("DirichletReg")
for (t in 2:(T+1)) {
  
  x[t] <-  which(rmultinom(1, 1, Px[x[t-1],]) == 1)
  
  # Simulate theta_t
  runge.kutta <- update.SEIR.rk4(theta[1,t-1],
                                 theta[2,t-1],
                                 theta[3,t-1],
                                 theta[4,t-1],
                                 alpha, 
                                 beta, 
                                 gamma,
                                 f[x[t],1], 
                                 pop.size)
  
  theta[,t] <- rdirichlet(1, alpha = 
                            kappa*c(runge.kutta$S.new, 
                                    runge.kutta$E.new, 
                                    runge.kutta$I.new, 
                                    runge.kutta$R.new))
  
  
  # Simulate y_t
  y[t] <- rbeta(1, lambda*p*theta[3,t], lambda*(1-p*theta[3,t]))
  
}

library(ggplot2)
library(RColorBrewer)
display.brewer.all()  # check color pallette

##### Draw simulated Y_t and I_t
pdf(paste0("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Simulation Study/Two Regime/TwoRegimeSimulatedData_seed", seed, "_T",T+1, ".pdf"), width = 12, height = 5)
par(mfrow=c(1,1))
plot(as.vector(y), 
     type = "p", 
     lty=1, 
     col = "gray50", 
     pch = 20, 
     ylab = "Proportion Infected", 
     xlab = "Time",
     cex = 1.5,
     cex.lab = 1.2)
lines(theta[3,]*p, type = "l", col = "gray20", lwd=2)

# Identify different regimes using colored background

## Observed true regimes
library("R.utils")
fromto.x1 <- seqToIntervals(which(x==1))
for (i in 1:nrow(fromto.x1)){
  rect(fromto.x1[i,1], -1, fromto.x1[i,2], 1, 
       col = rgb(1,1,1,1/4), border = rgb(1,1,1,1/4))
}

fromto.x2 <- seqToIntervals(which(x==2))
for (i in 1:nrow(fromto.x2)){
  rect(fromto.x2[i,1], -1, fromto.x2[i,2], 1, 
       col = rgb(0.9, 0.6, 0.6, 1/3), border = rgb(0.9, 0.6, 0.6, 1/3))
}  


# Add a legend to the plot
legend("topleft", 
       legend=c(expression(y[t]), expression(pI[t])), 
       col=c("gray50", "gray20"),
       lty=c(NA, 1), 
       lwd=c(NA, 2),
       pch=c(20, NA),
       cex=1.2)

dev.off()


##### Draw S_t, E_t, I_t, R_t
pdf(paste0("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Simulation Study/Two Regime/TwoRegimeSimulatedSEIR_seed", seed, "_T",T+1, ".pdf"), width = 10, height = 7)

# display.brewer.all() 
SEIR_colors <- c("#80B1D3", "#B3DE69", "#FB8072", "#BEBADA")  # set the number of colors and the palette name

# Convert theta to a data frame
theta_df <- data.frame(
  time = 1:length(theta[1,]),
  S = theta[1,],
  E = theta[2,],
  I = theta[3,],
  R = theta[4,]
)

# Convert the data frame to a long format
theta_long <- tidyr::gather(theta_df, key = "variable", value = "value", -time)
var_order <- c("S", "E", "I", "R")
theta_long$variable <- factor(theta_long$variable, levels = var_order)

# Create a plot with four facets
ggplot(theta_long, aes(x = time, y = value, color = variable)) +
  geom_line(size=0.8) +
  xlab("Time") + ylab("Proportion") +
  facet_wrap(~variable, nrow = 2, scales = "free_y") +
  scale_color_manual(values = SEIR_colors) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  ggtitle("Simulated Dynamics of SEIR") +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
  coord_cartesian(ylim = c(0, 1))  # Set the vertical limits here

dev.off()

#################### Save simulated data and model parameters ####################
# Save simulated data and model parameters
sim_data <- list(seed = seed, theta = theta, y = y, x = x, regimes = regimes,
                 pop.size = pop.size, alpha = alpha, beta = beta, gamma = gamma,
                 lambda = lambda, kappa = kappa, p = p, f = f, Px = Px, T = T)

saveRDS(sim_data, 
        file = paste0("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Data/Simulation Data/Two Regime/simulated_data_seed", seed, "_T", T+1,"_K",lenXset, ".RDS"))


# -------------------------- II. Run PG-CSMC-AS  ---------------------------------#

library(foreach)
library(doParallel)
detectCores()
cl <- makeCluster(6)
registerDoParallel(cl)

# Read in data
sim_data <- readRDS(file = "~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Data/Simulation Data/Two Regime/simulated_data_seed1580_T150_K2.RDS")
y <- sim_data$y
regimes <- sim_data$regimes
lenXset <- length(regimes)

# Number of MCMC chains
nchain  <- 2    

# Number of MCMC iterations
burnin <- 1000
niter   <- 10000 + burnin      # number of MCMC iterations

# Number of theta_t for each x_t 
M <- 50                  

# Set up hyperparameters
hyperparams <- list(
  m.alpha = 0.3,
  sigma.alpha = 0.1,
  m.beta = 0.4,
  sigma.beta = 0.1,
  m.gamma = 0.2,
  sigma.gamma = 0.1,
  a.lambda = 2,
  b.lambda = 0.001,
  a.kappa = 200,
  b.kappa = 0.01,
  a.p = 0.1,
  b.p = 0.4,
  m.p = 0.25,
  sigma.p = 0.05,
  delta.mat = matrix(c(10,1,1,10), nrow = lenXset, ncol = lenXset),
  a.f = 0,
  b.f = 1) 

# Load the script that contains the PG.CSMC.AS() function
source("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Code/Two Regime/ParticleGibbs_TwoRegime.R")

ptm <- proc.time() 
PG.results <- foreach(i = 1:nchain, .combine = "list", .packages = c("truncnorm", "DirichletReg")) %dopar% {
  # Run PG sampler
  PG.CSMC.AS(y, regimes, M, niter, hyperparams, pop.size=1)
  
} 
proc.time()-ptm 

stopCluster(cl)

# Show contents of PG.results
str(PG.results)

# Save results to local file
setwd("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Output/Simulation Study/Two Regime")
saveRDS(hyperparams, paste0("PG_results_hyperparams_pUnknown_seed", sim_data$seed, "_niter", niter,"_M", M, "_K", length(regimes), ".rds"))
saveRDS(PG.results, paste0("PG_results_pUnknown_seed", sim_data$seed, "_niter", niter, "_M", M, "_K", length(regimes), ".rds"))

# # Read results
# hyperparams <- readRDS("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Output/Simulation Study/Two Regime/PG_results_hyperparams_seed1580_niter11000_M50_K2.rds")
PG.results <- readRDS("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Output/Simulation Study/Two Regime/PG_results_pUnknown_seed1580_niter11000_M50_K2.rds")

# Extract results from each chain 
MCMC.chain.1 <- PG.results[[1]]
MCMC.chain.2 <- PG.results[[2]]

# Chain 1
SsampleMat.CSMC.AS.repM <- MCMC.chain.1$SsampleMat.CSMC.AS.repM
EsampleMat.CSMC.AS.repM <- MCMC.chain.1$EsampleMat.CSMC.AS.repM
IsampleMat.CSMC.AS.repM <- MCMC.chain.1$IsampleMat.CSMC.AS.repM
RsampleMat.CSMC.AS.repM <- MCMC.chain.1$RsampleMat.CSMC.AS.repM
XsampleMat.CSMC.AS.repM <- MCMC.chain.1$XsampleMat.CSMC.AS.repM
parameters.CSMC.AS.repM <- MCMC.chain.1$parameters.CSMC.AS.repM
marginalLogLik.CSMC.AS.repM <- MCMC.chain.1$marginalLogLik.CSMC.AS.repM

# Chain 2
SsampleMat.CSMC.AS.repM <- MCMC.chain.2$SsampleMat.CSMC.AS.repM
EsampleMat.CSMC.AS.repM <- MCMC.chain.2$EsampleMat.CSMC.AS.repM
IsampleMat.CSMC.AS.repM <- MCMC.chain.2$IsampleMat.CSMC.AS.repM
RsampleMat.CSMC.AS.repM <- MCMC.chain.2$RsampleMat.CSMC.AS.repM
XsampleMat.CSMC.AS.repM <- MCMC.chain.2$XsampleMat.CSMC.AS.repM
parameters.CSMC.AS.repM <- MCMC.chain.2$parameters.CSMC.AS.repM
marginalLogLik.CSMC.AS.repM <- MCMC.chain.2$marginalLogLik.CSMC.AS.repM


############################ Data Visualization #######################################
sim_data <- readRDS("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Data/Simulation Data/Two Regime/simulated_data_seed1580_T150_K2.RDS")
seed <- sim_data$seed
alpha <- sim_data$alpha
beta <- sim_data$beta
gamma <- sim_data$gamma
lambda <- sim_data$lambda
kappa <- sim_data$kappa
p <- sim_data$p
f <- sim_data$f
Px <- sim_data$Px
y <- sim_data$y
theta <- sim_data$theta
x <- sim_data$x
T <- sim_data$T
regimes <- sim_data$regimes
lenXset <- length(sim_data$regimes)

# Trace plot of parameters'
pdf(paste0("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Simulation Study/Two Regime/", Sys.Date(), " SimulationTwoRegimeTracePlot_seed", seed, "_T", T+1, ".pdf"), width = 12, height = 10)
par(mfrow=c(4,2))
plot(parameters.CSMC.AS.repM$alpha[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(alpha), panel.first=abline(h = alpha, col = "red"))
plot(parameters.CSMC.AS.repM$beta[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(beta), panel.first=abline(h = beta, col = "red"))
plot(parameters.CSMC.AS.repM$gamma[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(gamma), panel.first=abline(h = gamma, col = "red"))
plot(parameters.CSMC.AS.repM$kappa[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(kappa), panel.first=abline(h = kappa, col = "red"), ylim=c(1000, 12000))
plot(parameters.CSMC.AS.repM$lambda[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(lambda), panel.first=abline(h = lambda, col = "red"))
plot(parameters.CSMC.AS.repM$p[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(p), panel.first=abline(h = p, col = "red"))
# plot(parameters.CSMC.AS.repM$f[1, ], type="l", xlab="Iterations", ylab = expression(f[1]), panel.first=abline(h = f[1,1], col = "red"))
plot(parameters.CSMC.AS.repM$f[2, burnin:niter], type="l", xlab="Iterations", ylab = expression(f[2]), panel.first=abline(h = f[2,1], col = "red"))

post.pi.k1 <- matrix(0, niter, lenXset, byrow=TRUE)
post.pi.k2 <- matrix(0, niter, lenXset, byrow=TRUE)
for (i in 2:(niter)){
  post.pi.k1[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][1,]
  post.pi.k2[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][2,]
}
post.pi.k1
post.pi.k2

# pi_{k1}
plot(post.pi.k1[burnin:(niter-1),1],
     type="l",
     xlab="Iterations",
     ylab = expression(pi[11]),
     # ylim = c(0.6, 1),
     panel.first=abline(h = Px[1,1], col = "red")) # pi_11
# plot(post.pi.k1[1:(niter-1),2],
#      type="l",
#      xlab="Iterations",
#      ylab = expression(pi[12]),
#      #  ylim = c(0, 0.5),
#      panel.first=abline(h = Px[1,2], col = "red")) # pi_12


# pi_{k2}
# plot(post.pi.k2[1:(niter-1),1],
#      type="l",
#      xlab="Iterations",
#      ylab = expression(pi[21]),
#      panel.first=abline(h = Px[2,1], col = "red")) # pi_21
plot(post.pi.k2[burnin:(niter-1),2],
     type="l",
     xlab="Iterations",
     ylab = expression(pi[22]),
     panel.first=abline(h = Px[2,2], col = "red"), ylim=c(0,1))

dev.off()


# Histogram of parameters
pdf(paste0("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Simulation Study/Two Regime/", Sys.Date(), " SimulationTwoRegimeHistogram_seed", seed, "_T", T+1, ".pdf"), width = 10, height = 7)
par(mfrow=c(3,3))

# Define colors
blue <- "#0072B2"
red <- "#D55E00"

# Set up the plot
par(mar = c(4, 4, 1, 1))
hist(
  parameters.CSMC.AS.repM$alpha[burnin:niter],
  xlab = expression(alpha),
  main = "",
  col = "lightgrey",
  border = "white",
  freq = FALSE
)
axis(1, col = "gray70")
axis(2, col = "gray70")
abline(v = alpha, col = red, lwd = 2)
# abline(v = median(parameters.CSMC.AS.repM$alpha[burnin:niter]), col = blue, lwd = 2)
abline(v = mean(parameters.CSMC.AS.repM$alpha[burnin:niter]), col = blue, lwd = 2)
alpha.CI <- quantile(parameters.CSMC.AS.repM$alpha[burnin:niter], c(0.025, 0.975))
abline(v = alpha.CI, col = blue, lty = 2, lwd = 2)


# beta
hist(
  parameters.CSMC.AS.repM$beta[burnin:niter],
  xlab = expression(beta),
  main = "",
  col = "lightgrey",
  border = "white",
  freq = FALSE
)
axis(1, col = "gray70")
axis(2, col = "gray70")
abline(v = beta, col = red, lwd = 2)
abline(v = median(parameters.CSMC.AS.repM$beta[burnin:niter]), col = blue, lwd = 2)
beta.CI <- quantile(parameters.CSMC.AS.repM$beta[burnin:niter], c(0.025, 0.975))
abline(v = beta.CI, col = blue, lty = 2, lwd = 2)

# gamma
hist(
  parameters.CSMC.AS.repM$gamma[burnin:niter],
  xlab = expression(gamma),
  main = "",
  col = "lightgrey",
  border = "white",
  freq = FALSE
)
axis(1, col = "gray70")
axis(2, col = "gray70")
abline(v = gamma, col = red, lwd = 2)
# abline(v = median(parameters.CSMC.AS.repM$gamma[burnin:niter]), col = blue, lwd = 2)
abline(v = mean(parameters.CSMC.AS.repM$gamma[burnin:niter]), col = blue, lwd = 2)
gamma.CI <- quantile(parameters.CSMC.AS.repM$gamma[burnin:niter], c(0.025, 0.975))
abline(v = gamma.CI, col = blue, lty = 2, lwd = 2)


# kappa
hist(
  parameters.CSMC.AS.repM$kappa[burnin:niter],
  xlab = expression(kappa),
  xlim = c(600, 10200),
  main = "",
  col = "lightgrey",
  border = "white",
  freq = FALSE
)
axis(1, col = "gray70")
axis(2, col = "gray70")
abline(v = kappa, col = red, lwd = 2)
# abline(v = median(parameters.CSMC.AS.repM$kappa[burnin:niter]), col = blue, lwd = 2)
abline(v = mean(parameters.CSMC.AS.repM$kappa[burnin:niter]), col = blue, lwd = 2)
kappa.CI <- quantile(parameters.CSMC.AS.repM$kappa[burnin:niter], c(0.025, 0.975))
abline(v = kappa.CI, col = blue, lty = 2, lwd = 2)

# lambda
hist(
  parameters.CSMC.AS.repM$lambda[burnin:niter],
  xlab = expression(lambda),
  main = "",
  col = "lightgrey",
  border = "white",
  freq = FALSE
)
axis(1, col = "gray70")
axis(2, col = "gray70")
abline(v = lambda, col = red, lwd = 2)
# abline(v = median(parameters.CSMC.AS.repM$lambda[burnin:niter]), col = blue, lwd = 2)
abline(v = mean(parameters.CSMC.AS.repM$lambda[burnin:niter]), col = blue, lwd = 2)
lambda.CI <- quantile(parameters.CSMC.AS.repM$lambda[burnin:niter], c(0.025, 0.975))
abline(v = lambda.CI, col = blue, lty = 2, lwd = 2)

# p
hist(
  parameters.CSMC.AS.repM$p[burnin:niter],
  xlab = expression(p),
  main = "",
  col = "lightgrey",
  border = "white",
  freq = FALSE
)
axis(1, col = "gray70")
axis(2, col = "gray70")
abline(v = p, col = red, lwd = 2)
# abline(v = median(parameters.CSMC.AS.repM$p[burnin:niter]), col = blue, lwd = 2)
abline(v = mean(parameters.CSMC.AS.repM$p[burnin:niter]), col = blue, lwd = 2)
p.CI <- quantile(parameters.CSMC.AS.repM$p[burnin:niter], c(0.025, 0.975))
abline(v = p.CI, col = blue, lty = 2, lwd = 2)

# f_2
hist(
  parameters.CSMC.AS.repM$f[2,burnin:niter],
  xlab = expression(f[2]),
  main = "",
  col = "lightgrey",
  border = "white",
  freq = FALSE
)
axis(1, col = "gray70")
axis(2, col = "gray70")
abline(v = f[2,1], col = red, lwd = 2)
# abline(v = median(parameters.CSMC.AS.repM$f[2,burnin:niter]), col = blue, lwd = 2)
abline(v = mean(parameters.CSMC.AS.repM$f[2,burnin:niter]), col = blue, lwd = 2)
f2.CI <- quantile(parameters.CSMC.AS.repM$f[2,burnin:niter], c(0.025, 0.975))
abline(v = f2.CI, col = blue, lty = 2, lwd = 2)

# pi_11
hist(
  post.pi.k1[burnin:(niter-1),1],
  xlab = expression(pi[11]),
  main = "",
  col = "lightgrey",
  border = "white",
  freq = FALSE
)
axis(1, col = "gray70")
axis(2, col = "gray70")
abline(v = Px[1,1], col = red, lwd = 2)
# abline(v = median(post.pi.k1[burnin:(niter-1),1]), col = blue, lwd = 2)
abline(v = mean(post.pi.k1[burnin:(niter-1),1]), col = blue, lwd = 2)
pi11.CI <- quantile(post.pi.k1[burnin:(niter-1),1], c(0.025, 0.975))
abline(v = pi11.CI, col = blue, lty = 2, lwd = 2)


# pi_22
hist(
  post.pi.k2[burnin:(niter-1),2],
  xlab = expression(pi[22]),
  main = "",
  col = "lightgrey",
  border = "white",
  freq = FALSE
)
axis(1, col = "gray70")
axis(2, col = "gray70")
abline(v = Px[2,2], col = red, lwd = 2)
# abline(v = median(post.pi.k2[burnin:(niter-1),2]), col = blue, lwd = 2)
abline(v = mean(post.pi.k2[burnin:(niter-1),2]), col = blue, lwd = 2)
pi22.CI <- quantile(post.pi.k2[burnin:(niter-1),2], c(0.025, 0.975))
abline(v = pi22.CI, col = blue, lty = 2, lwd = 2)

dev.off()



## 95% Credible Interval of posterior samples
S.CredInterval <- t(apply(SsampleMat.CSMC.AS.repM[burnin:(niter-1),], 2, function(x) quantile(x, c(0.025, 0.975))))
S.CredInterval.LL <- S.CredInterval[,1]
S.CredInterval.UL <- S.CredInterval[,2]
E.CredInterval <- t(apply(EsampleMat.CSMC.AS.repM[burnin:(niter-1),], 2, function(x) quantile(x, c(0.025, 0.975))))
E.CredInterval.LL <- E.CredInterval[,1]
E.CredInterval.UL <- E.CredInterval[,2]
I.CredInterval <- t(apply(IsampleMat.CSMC.AS.repM[burnin:(niter-1),], 2, function(x) quantile(x, c(0.025, 0.975))))
I.CredInterval.LL <- I.CredInterval[,1]
I.CredInterval.UL <- I.CredInterval[,2]
R.CredInterval <- t(apply(RsampleMat.CSMC.AS.repM[burnin:(niter-1),], 2, function(x) quantile(x, c(0.025, 0.975))))
R.CredInterval.LL <- R.CredInterval[,1]
R.CredInterval.UL <- R.CredInterval[,2]


pdf(paste0("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Simulation Study/Two Regime/", Sys.Date(), " TwoRegimeEstimatedSEIR_seed", seed, "_T", T+1, ".pdf"), width = 10, height = 10)

par(mfrow=c(2,2), mar=c(4, 4, 2, 2), mgp=c(2, 0.7, 0))


# # S plot
# plot(colMeans(SsampleMat.CSMC.AS.repM),
#      xlab="Time",
#      type = "l",
#      ylab = "Proportion",
#      col= "darkgrey",
#      lwd = 1.5,
#      ylim=c(0,1),
#      main="Susceptible (S)")
# lines(as.vector(theta[1,]), col= "black", lwd = 1.5, lty = 1)
# lines(S.CredInterval.LL, col= "darkgrey", lwd = 1.5, lty = 2)
# lines(S.CredInterval.UL, col= "darkgrey", lwd = 1.5, lty = 2)
# 
# # E plot
# plot(colMeans(EsampleMat.CSMC.AS.repM),
#      xlab="Time",
#      type = "l",
#      ylab = "Proportion",
#      col="darkgrey",
#      lwd = 1.5,
#      ylim=c(0, 0.1),
#      main="Exposed (E)")
# lines(as.vector(theta[2,]), col= "black", lwd = 1.5, lty = 1)
# lines(E.CredInterval.LL, col= "darkgrey", lwd = 1.5, lty = 2)
# lines(E.CredInterval.UL, col= "darkgrey", lwd = 1.5, lty = 2)
# 
# # I plot
# plot(colMeans(IsampleMat.CSMC.AS.repM),
#      xlab="Time",
#      type = "l",
#      ylab = "Proportion",
#      col="darkgrey",
#      ylim=c(0,0.1),
#      main="Infected (I)")
# lines(as.vector(theta[3,]), col= "black", lwd = 1.5, lty = 1)
# lines(I.CredInterval.LL, col= "darkgrey", lwd = 1.5, lty = 2)
# lines(I.CredInterval.UL, col= "darkgrey", lwd = 1.5, lty = 2)
# 
# # R plot
# plot(colMeans(RsampleMat.CSMC.AS.repM),
#      xlab="Time",
#      type = "l",
#      ylab = "Proportion",
#      col="darkgrey",
#      ylim=c(0,1),
#      main="Recovered (R)")
# lines(as.vector(theta[4,]), col= "black", lwd = 1.5, lty = 1)
# lines(R.CredInterval.LL, col= "darkgrey", lwd = 1.5, lty = 2)
# lines(R.CredInterval.UL, col= "darkgrey", lwd = 1.5, lty = 2)
# 
# # Increase font size of plot labels and titles
# par(cex.lab=1.2, cex.main=1.5)
# 
# # Add legend
# legend("topright",
#        legend = c("True Value", "Posterior Mean", "95% CI"),
#        lty = c(1, 2, 2),
#        lwd = c(1.5, 1.5, 1.5),
#        col = c("black", "darkgrey", "darkgrey"),
#        bty = "n",
#        cex = 1.2)

## S plot
# Compute the x and y coordinates for the shaded area
xx <- 1:length(S.CredInterval.LL)
y_lower <- S.CredInterval.LL
y_upper <- S.CredInterval.UL

# Plot the mean trajectory
plot(colMeans(SsampleMat.CSMC.AS.repM[burnin:(niter-1),]),
     xlab = "Time",
     type = "l",
     ylab = "Proportion",
     col = "#969696",
     lty = 2,
     lwd = 2,
     ylim = c(0, 1),
     main = "Susceptible (S)",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 1.5)

# Plot the mean trajectory line
lines(as.vector(theta[1,]), col = "black", lwd = 1.5, lty = 1)

# Plot the shaded area for the credible interval
polygon(c(xx, rev(xx)), c(y_lower, rev(y_upper)), col = rgb(0.5, 0.5, 0.5, alpha = 1/4), border = NA)


## E plot
# Compute the x and y coordinates for the shaded area
xx <- 1:length(E.CredInterval.LL)
y_lower <- E.CredInterval.LL
y_upper <- E.CredInterval.UL

# Plot the mean trajectory
plot(colMeans(EsampleMat.CSMC.AS.repM[burnin:(niter-1),]),
     xlab = "Time",
     type = "l",
     ylab = "Proportion",
     col = "#969696",
     lty = 2,
     lwd = 2,
     ylim = c(0, 0.1),
     main = "Exposed (E)",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 1.5)

# Plot the mean trajectory line
lines(as.vector(theta[2,]), col = "black", lwd = 1.5, lty = 1)

# Plot the shaded area for the credible interval
polygon(c(xx, rev(xx)), 
        c(y_lower, rev(y_upper)), 
        col = rgb(0.5, 0.5, 0.5, alpha = 1/4), 
        border = NA)


## I plot
# Compute the x and y coordinates for the shaded area
xx <- 1:length(I.CredInterval.LL)
y_lower <- I.CredInterval.LL
y_upper <- I.CredInterval.UL

# Plot the mean trajectory
plot(colMeans(IsampleMat.CSMC.AS.repM[burnin:(niter-1),]),
     xlab = "Time",
     type = "l",
     ylab = "Proportion",
     col = "#969696",
     lty = 2,
     lwd = 2,
     ylim = c(0, 0.1),
     main = "Infected (I)",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 1.5)

# Plot the mean trajectory line
lines(as.vector(theta[3,]), col = "black", lwd = 1.5, lty = 1)

# Plot the shaded area for the credible interval
polygon(c(xx, rev(xx)), 
        c(y_lower, rev(y_upper)), 
        col = rgb(0.5, 0.5, 0.5, alpha = 1/4), 
        border = NA)


## R plot
# Compute the x and y coordinates for the shaded area
xx <- 1:length(R.CredInterval.LL)
y_lower <- R.CredInterval.LL
y_upper <- R.CredInterval.UL

# Plot the mean trajectory
plot(colMeans(RsampleMat.CSMC.AS.repM[burnin:(niter-1),]),
     xlab = "Time",
     type = "l",
     ylab = "Proportion",
     col = "#969696",
     lty = 2,
     lwd = 2,
     ylim = c(0, 1),
     main = "Recovered (R)",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 1.5)

# Plot the mean trajectory line
lines(as.vector(theta[4,]), col = "black", lwd = 1.5, lty = 1)

# Plot the shaded area for the credible interval
polygon(c(xx, rev(xx)), 
        c(y_lower, rev(y_upper)), 
        col = rgb(0.5, 0.5, 0.5, alpha = 1/4), 
        border = NA)


# Add legend
legend("topright",
       legend = c("True Value", "Posterior Mean"),
       lty = c(1, 2),
       lwd = c(1.5, 1.5),
       col = c("black", "#969696"),
       bty = "n",
       cex = 1.2)

dev.off()



# Posterior I_t
pdf(paste0("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Two Regime/", Sys.Date(), " TwoRegimeEstimatedIt_seed", seed, "_T", T+1, ".pdf"), width = 10, height = 7)
par(mfrow=c(1,1))
plot(1:length(y), y, type = "p", col = "grey", pch=20, ylab = "Proportion Infected", ylim=c(0,0.03))
# lines(1:(T+1), theta[3,]*p, col="red")
lines(1:length(y), apply(IsampleMat.CSMC.AS.repM*mean(parameters.CSMC.AS.repM$p), 2, mean), col="grey")
lines(I.CredInterval.LL*mean(parameters.CSMC.AS.repM$p), col="grey", lty = 2)
lines(I.CredInterval.UL*mean(parameters.CSMC.AS.repM$p), col="grey", lty = 2)
dev.off()



## Observed true regimes
pdf(paste0("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Simulation Study/Two Regime/", Sys.Date(), " TwoRegimeEstimatedXt.pdf"), width = 12, height = 8)
# par(mfrow=c(2,1))
# 
# plot(1:length(y), y, type = "p", 
#      lty = 1, col = "#BDBDBD", pch = 20, 
#      ylim = c(0,0.02), xlab = "Time",
#      main = "True Regimes", cex = 1.5)
# lines(theta[3,]*p, type = "l", col = "indianred", lwd=2)
# 
# library("R.utils")
# fromto.x1 <- seqToIntervals(which(x==1))
# for (i in 1:nrow(fromto.x1)){
#   rect(fromto.x1[i,1], -1, fromto.x1[i,2], 1, 
#        col = rgb(1,1,1,1/4), border = rgb(1,1,1,1/4))
# }
# 
# fromto.x2 <- seqToIntervals(which(x==2))
# for (i in 1:nrow(fromto.x2)){
#   rect(fromto.x2[i,1], -1, fromto.x2[i,2], 1, 
#        col = rgb(0.5,0.5,0.5,1/3), border = rgb(0.5,0.5,0.5,1/3))
# }  
# 
# # Estimated X_t
# plot(1:length(y), y, type = "p", lty = 1, 
#      col = "#BDBDBD", pch = 20,  ylim = c(0,0.02), 
#      xlab = "Time", main = "Estimated Regimes", cex = 1.5)
# lines(1:(T+1), apply(IsampleMat.CSMC.AS.repM*mean(parameters.CSMC.AS.repM$p), 2, median), col="black")
# lines(I.CredInterval.LL*mean(parameters.CSMC.AS.repM$p), col="black", lty = 2)
# lines(I.CredInterval.UL*mean(parameters.CSMC.AS.repM$p), col="black", lty = 2)
# 
# est.prob.x <- data.frame("prob.x1" = colSums(XsampleMat.CSMC.AS.repM == 1)/niter,
#                          "prob.x2" = colSums(XsampleMat.CSMC.AS.repM == 2)/niter)
# est.x <- apply(est.prob.x, 1, which.max)
# 
# library("R.utils")
# fromto.x1 <- seqToIntervals(which(est.x==1))
# for (i in 1:nrow(fromto.x1)){
#   rect(fromto.x1[i,1], -1, fromto.x1[i,2], 1, 
#        col = rgb(1,1,1,1/4), border = rgb(1,1,1,1/4))
# }
# 
# fromto.x2 <- seqToIntervals(which(est.x==2))
# for (i in 1:nrow(fromto.x2)){
#   rect(fromto.x2[i,1], -1, fromto.x2[i,2], 1, 
#        col = rgb(0.5,0.5,0.5,1/3), border = rgb(0.5,0.5,0.5,1/3))
# }  


# Create a new plot with two subplots side by side
par(mfrow = c(2, 1))

# Plot True Regimes
plot(1:length(y), y, type = "p", 
     lty = 1, col = "gray50", pch = 20, 
     ylim = c(0, 0.02), xlab="",
     main = "True Regimes vs. Estimated Regimes", cex = 1.5)
lines(theta[3,]*p, type = "l", col = "gray20", lwd = 2)

library("R.utils")
fromto.x1 <- seqToIntervals(which(x == 1))
for (i in 1:nrow(fromto.x1)){
  rect(fromto.x1[i, 1], -1, fromto.x1[i, 2], 1, 
       col = rgb(1, 1, 1, 1/4), border = rgb(1, 1, 1, 1/4))
}

fromto.x2 <- seqToIntervals(which(x == 2))
for (i in 1:nrow(fromto.x2)){
  rect(fromto.x2[i, 1], -1, fromto.x2[i, 2], 1, 
       col = rgb(0.9, 0.6, 0.6, 1/3), border = rgb(0.9, 0.6, 0.6, 1/3))
}  

# Plot Estimated Regimes
plot(1:length(y), y, type = "p", lty = 1, 
     col = "gray50", pch = 20,  ylim = c(0, 0.02), 
     xlab = "Time", main = "", cex = 1.5)
lines(1:(T+1), apply(IsampleMat.CSMC.AS.repM*mean(parameters.CSMC.AS.repM$p), 2, median), col = "#969696", lwd=2.5, lty=1)
# lines(I.CredInterval.LL*mean(parameters.CSMC.AS.repM$p), col = "#969696", lty = 2, lwd=1.5)
# lines(I.CredInterval.UL*mean(parameters.CSMC.AS.repM$p), col = "#969696", lty = 2, lwd=1.5)

# Compute the x and y coordinates for the shaded area
xx <- 1:length(apply(IsampleMat.CSMC.AS.repM*mean(parameters.CSMC.AS.repM$p), 2, median))
y_lower <- I.CredInterval.LL*mean(parameters.CSMC.AS.repM$p)
y_upper <- I.CredInterval.UL*mean(parameters.CSMC.AS.repM$p)

# Plot the shaded area for the credible interval
polygon(c(xx, rev(xx)), 
        c(y_lower, rev(y_upper)), 
        col = rgb(0.5, 0.5, 0.5, alpha = 1/4), 
        border = NA)

# Plot regimes
est.prob.x <- data.frame("prob.x1" = colSums(XsampleMat.CSMC.AS.repM == 1)/niter,
                         "prob.x2" = colSums(XsampleMat.CSMC.AS.repM == 2)/niter)
est.x <- apply(est.prob.x, 1, which.max)

library("R.utils")
fromto.x1 <- seqToIntervals(which(est.x == 1))
for (i in 1:nrow(fromto.x1)){
  rect(fromto.x1[i, 1], -1, fromto.x1[i, 2], 1, 
       col = rgb(1, 1, 1, 1/4), border = rgb(1, 1, 1, 1/4))
}

fromto.x2 <- seqToIntervals(which(est.x == 2))
for (i in 1:nrow(fromto.x2)){
  rect(fromto.x2[i, 1], -1, fromto.x2[i, 2], 1, 
       col = rgb(0.9, 0.6, 0.6, 1/3), border = rgb(0.9, 0.6, 0.6, 1/3))
}  
                 
                 
# # Estimated X_t
# plot(colSums(XsampleMat.CSMC.AS.repM == 1)/niter,
#      type = "l",
#      # xlim = c(20, T),
#      ylim = c(0, 1),            # Adjust accordingly!
#      xlab = paste("Data of length T =", T+1),
#      ylab = "Probability",
#      main = expression('Estimated P(X'[t]*'|y'[1:T]*')'),
#      col = "grey")
# lines(colSums(XsampleMat.CSMC.AS.repM == 2)/niter, col= blue)
# 
# # Add legend
# legend("topright",
#        legend=c(expression('P(X'[t]*'=1|y'[1:T]*')'), expression('P(X'[t]*'=2|y'[1:T]*')')),
#        col=c("grey", blue),
#        lty=1,
#        bty="n",
#        cex=0.8)
dev.off()


## Marginal likelihood
mean(marginalLogLik.CSMC.AS.repM)
plot(1:niter, marginalLogLik.CSMC.AS.repM, type="l")

