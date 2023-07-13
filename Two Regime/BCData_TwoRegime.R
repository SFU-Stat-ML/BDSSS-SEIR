
#------------------------ Clear memory and graphics ---------------------------#
rm(list=ls())
graphics.off()

options(digits = 16)

# ----------------Set the working directory to current path -------------------#
library("rstudioapi")
cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

################################ Daily active counts ###########################################

bc_data <- read.csv(file = "~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Data/Real Data/BC COVID CASES - Daily Cases.csv",
                    header = TRUE,
                    sep = ",",
                    stringsAsFactors = FALSE)
bc_data$Date <- as.Date(bc_data$Date, format="%Y-%m-%d")
bc_data$Active.cases <- as.integer(bc_data$Active.cases)
str(bc_data)
head(bc_data)

bc_data_active <- data.frame("Date" = bc_data$Date, "ActiveCases" = bc_data$Active.cases)
bc_data_active <- na.omit(bc_data_active) # Remove the rows with missing values
str(bc_data_active)
head(bc_data_active)
summary(bc_data_active$ActiveCases)

library("ggplot2")
g <- ggplot(bc_data_active, aes(x = Date, y = ActiveCases)) +
  geom_line(color = "steelblue") + 
  geom_vline(xintercept = as.numeric(as.Date("2021-11-26")), 
             linetype = "dashed", color = "red", size = 1) +
  geom_point(size=0.5) +
  xlab("") +
  ylab("Active Case counts in B.C") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.8, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  scale_x_date(limit = c(as.Date(min(bc_data_active$Date)), as.Date(max(bc_data_active$Date)))) +
  annotate("text", x = as.Date("2021-11-30"), y = max(bc_data_active$ActiveCases), 
           label = "1st case of the COVID-19 omicron variant", vjust = -0.3, hjust = 0.8,
           size = 5, color = "red")
g

ggsave(file="~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Real Data/BCDailyActiveCases.png", g)


################################ Weekly active counts ###########################################

# Get weekly active counts
bc_data_active$Week <- as.Date(cut(bc_data_active$Date,
                                   breaks = "week",
                                   start.on.monday = TRUE))

bc_data_active[1:200,]

# Extract data
library(dplyr)
# start <- which(bc_data_active$Date == "2020-07-20")
# end <- which(bc_data_active$Date == "2021-08-03")
new_data_weekly <- as.data.frame(bc_data_active%>%
                                   group_by(Week) %>%
                                   dplyr::summarise(weekly.n = sum(ActiveCases)))
new_data_weekly

g.weekly <- ggplot(new_data_weekly, aes(x = Week, y = weekly.n)) +
  geom_line(color="grey") + 
  geom_vline(xintercept = as.numeric(as.Date("2021-11-26")), 
             linetype = "dashed", color = "red", size = 0.8) +
  geom_point(color="gray50") +
  xlab("") +
  ylab("Weekly Active Cases in B.C") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_x_date(limit = c(as.Date(min(bc_data_active$Date)), as.Date(max(bc_data_active$Date)))) +
  annotate("text", x = as.Date("2021-11-30"), y = max(bc_data_active$ActiveCases), 
           label = "1st case of the COVID-19 Omicron variant", vjust = -20, hjust = 1.1,
           size = 5, color = "red")
g.weekly

ggsave(file="~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Real Data/BCWeeklyActiveCases.png", g.weekly, width = 7, height = 5)

# scale_x_date(limit=c(as.Date("2020-07-20"),as.Date("2021-07-02"))) # extract certain period

############################ Run PG Sampler ##############################

library(foreach)
library(doParallel)
detectCores()
cl <- makeCluster(6)
registerDoParallel(cl)

## Daily active cases (before omicron starts)
# y <- bc_data_active$ActiveCases[1:438]/5070000
# time <- bc_data_active$Date[1:438]

## Weekly active cases
y <- new_data_weekly$weekly.n[1:95]/5070000
time <- new_data_weekly$Week[1:95]

## Regimes
regimes <- c(1, 2)
lenXset <- length(regimes)

## Number of MCMC chains
nchain  <- 2    

## Number of MCMC iterations
burnin <- 1000
niter   <- 20000 + burnin      # number of MCMC iterations

## Number of theta_t for each x_t 
M <- 50                  

## Set up hyperparameters
hyperparams <- list(
  m.alpha = 0.2*7,
  sigma.alpha = 0.5,
  m.beta = 0.4*7,
  sigma.beta = 0.5,
  m.gamma = 0.2*7,
  sigma.gamma = 0.5,
  a.lambda = 20,
  b.lambda = 0.01,
  a.kappa = 200,
  b.kappa = 0.01,
  a.p = 0.1,
  b.p = 0.4,
  m.p = 0.25,
  sigma.p = 0.05,
  delta.mat = matrix(c(10,1,1,10), nrow = lenXset, ncol = lenXset),
  a.f = 0,
  b.f = 1) 

# Call PG sampler
source("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Code/Two Regime/ParticleGibbs_TwoRegime.R")

ptm <- proc.time()
PG.results <- PG.CSMC.AS(y, regimes, M, niter, hyperparams, pop.size=1)
proc.time()-ptm

ptm <- proc.time() 
PG.results <- foreach(i = 1:nchain, .combine = "list", .packages = c("truncnorm", "DirichletReg")) %dopar% {
  # Run PG sampler
  PG.CSMC.AS(y, regimes, M, niter, hyperparams, pop.size=1)
  
} 
proc.time()-ptm 

stopCluster(cl)


# Save results to local file
setwd("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Output/Real Data")
saveRDS(hyperparams, paste0("PG_results_final_hyperparams_niter", niter,"_M", M, "_K", length(regimes), ".rds"))
saveRDS(PG.results, paste0("PG_results_final_niter", niter, "_M", M, "_K", length(regimes), ".rds"))

# Read results
# hyperparams <- readRDS("~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Output/Real Data/PG_results_hyperparams_BCWeekly_niter11000_M50_K2.rds")
PG.results <- readRDS("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Output/Real Data/PG_results_final_niter21000_M50_K2.rds")

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
setwd("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Real Data")

# Trace plot of parameters
pdf(paste("(Weekly)TwoRegimeTracePlot.pdf"), width = 10, height = 7)
par(mfrow=c(3,3))
plot(parameters.CSMC.AS.repM$alpha[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(alpha), cex=2)
plot(parameters.CSMC.AS.repM$beta[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(beta), cex=2)
plot(parameters.CSMC.AS.repM$gamma[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(gamma), cex=2)
plot(parameters.CSMC.AS.repM$kappa[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(kappa), cex=2)
plot(parameters.CSMC.AS.repM$lambda[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(lambda), cex=2)
plot(parameters.CSMC.AS.repM$p[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(p))
# plot(parameters.CSMC.AS.repM$f[1, 1:niter], type="l", xlab="Iterations", ylab = expression(f[1]))
plot(parameters.CSMC.AS.repM$f[2, burnin:niter], type="l", xlab="Iterations", ylab = expression(f[2]), cex=2)
# plot(parameters.CSMC.AS.repM$f[3, 1:niter], type="l", xlab="Iterations", ylab = expression(f[3]))


# Posterior Px
post.pi.k1 <- matrix(0, niter, lenXset, byrow=TRUE)
post.pi.k2 <- matrix(0, niter, lenXset, byrow=TRUE)
for (i in 1:(niter)){
  post.pi.k1[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][1,]
  post.pi.k2[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][2,]
}
post.pi.k1
post.pi.k2

# pi_{k1}
plot(post.pi.k1[burnin:(niter-1),1], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[11]), cex=2) # pi_11


plot(post.pi.k2[burnin:(niter-1),2], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[22]), cex=2)

dev.off()


# Histogram of parameters
pdf(paste("(Weekly)TwoRegimeHistogram.pdf"), width = 12, height = 7)
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
# abline(v = median(parameters.CSMC.AS.repM$beta[burnin:niter]), col = blue, lwd = 2)
abline(v = mean(parameters.CSMC.AS.repM$beta[burnin:niter]), col = blue, lwd = 2)
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
# abline(v = median(parameters.CSMC.AS.repM$gamma[burnin:niter]), col = blue, lwd = 2)
abline(v = mean(parameters.CSMC.AS.repM$gamma[burnin:niter]), col = blue, lwd = 2)
gamma.CI <- quantile(parameters.CSMC.AS.repM$gamma[burnin:niter], c(0.025, 0.975))
abline(v = gamma.CI, col = blue, lty = 2, lwd = 2)
  
  
# kappa
hist(
    parameters.CSMC.AS.repM$kappa[burnin:niter],
    xlab = expression(kappa),
    main = "",
    xlim = c(3000, 15000),
    col = "lightgrey",
    border = "white",
    freq = FALSE
)
axis(1, col = "gray70")
axis(2, col = "gray70")
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
# abline(v = median(parameters.CSMC.AS.repM$lambda[burnin:niter]), col = blue, lwd = 2)
abline(v = mean(parameters.CSMC.AS.repM$lambda[burnin:niter]), col = blue, lwd = 2)
lambda.CI <- quantile(parameters.CSMC.AS.repM$lambda[burnin:niter], c(0.025, 0.975))
abline(v = lambda.CI, col = blue, lty = 2, lwd = 2)
  
# p
hist(
    parameters.CSMC.AS.repM$p[1,burnin:niter],
    xlab = expression(p),
    main = "",
    col = "lightgrey",
    border = "white",
    freq = FALSE
)
axis(1, col = "gray70")
axis(2, col = "gray70")
# abline(v = median(parameters.CSMC.AS.repM$p[1,burnin:niter]), col = blue, lwd = 2)
abline(v = mean(parameters.CSMC.AS.repM$p[1,burnin:niter]), col = blue, lwd = 2)
p.CI <- quantile(parameters.CSMC.AS.repM$p[1,1:niter], c(0.025, 0.975))
abline(v = p.CI, col = blue, lty = 2, lwd = 2)
  
# f2
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
# abline(v = median(parameters.CSMC.AS.repM$f[2,burnin:niter]), col = blue, lwd = 2)
abline(v = mean(parameters.CSMC.AS.repM$f[2,burnin:niter]), col = blue, lwd = 2)
f2.CI <- quantile(parameters.CSMC.AS.repM$f[2,1:niter], c(0.025, 0.975))
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
# abline(v = median(post.pi.k2[burnin:(niter-1),2]), col = blue, lwd = 2)
abline(v = mean(post.pi.k2[burnin:(niter-1),2]), col = blue, lwd = 2)
pi22.CI <- quantile(post.pi.k2[burnin:(niter-1),2], c(0.025, 0.975))
abline(v = pi22.CI, col = blue, lty = 2, lwd = 2)
  
dev.off()
  
  
  
  
  
  ## 95% Credible Interval of posterior samples
  S.CredInterval <- t(apply(SsampleMat.CSMC.AS.repM[burnin:niter,], 2, function(x) quantile(x, c(0.025, 0.975))))
  S.CredInterval.LL <- S.CredInterval[,1]
  S.CredInterval.UL <- S.CredInterval[,2]
  E.CredInterval <- t(apply(EsampleMat.CSMC.AS.repM[burnin:niter,], 2, function(x) quantile(x, c(0.025, 0.975))))
  E.CredInterval.LL <- E.CredInterval[,1]
  E.CredInterval.UL <- E.CredInterval[,2]
  I.CredInterval <- t(apply(IsampleMat.CSMC.AS.repM[burnin:niter,], 2, function(x) quantile(x, c(0.025, 0.975))))
  I.CredInterval.LL <- I.CredInterval[,1]
  I.CredInterval.UL <- I.CredInterval[,2]
  R.CredInterval <- t(apply(RsampleMat.CSMC.AS.repM[burnin:niter,], 2, function(x) quantile(x, c(0.025, 0.975))))
  R.CredInterval.LL <- R.CredInterval[,1]
  R.CredInterval.UL <- R.CredInterval[,2]
  
  pdf(paste("(Weekly)TwoRegimeEstimatedSEIR.pdf"), width = 12, height = 10)
  par(mfrow=c(2,2))
  
  # plot(as.Date(time), apply(SsampleMat.CSMC.AS.repM[burnin:niter,], 2, mean),
  #      xlab = "",
  #      type = "l",
  #      ylab = "S",
  #      col = "grey",
  #      ylim = c(0,1),
  #      xaxt = "n",
  #      lwd = 2)
  # lines(time, S.CredInterval.LL, col= "darkgrey", lty = 2, lwd = 2)
  # lines(time, S.CredInterval.UL, col= "darkgrey", lty = 2, lwd = 2)
  # axis(1, at = time, lwd = 0, labels = format(time, "%b\n%Y"))
  # 
  # plot(time, apply(EsampleMat.CSMC.AS.repM[burnin:niter,], 2, mean),
  #      xlab = "",
  #      type = "l",
  #      ylab = "E",
  #      col = "grey",
  #      ylim = c(0,0.1),
  #      xaxt = "n",
  #      lwd = 2)
  # lines(time, E.CredInterval.LL, col= "darkgrey", lty = 2, lwd = 2)
  # lines(time, E.CredInterval.UL, col= "darkgrey", lty = 2, lwd = 2)
  # axis(1, at = time, lwd = 0, labels = format(time, "%b\n%Y"))
  # 
  # plot(time, apply(IsampleMat.CSMC.AS.repM[burnin:niter,], 2, mean),
  #      xlab = "",
  #      type = "l",
  #      ylab = "I",
  #      col = "grey",
  #      ylim=c(0,0.1),
  #      xaxt = "n",
  #      lwd = 2)
  # lines(time, I.CredInterval.LL, col= "darkgrey", lty = 2, lwd = 2)
  # lines(time, I.CredInterval.UL, col= "darkgrey", lty = 2, lwd = 2)
  # axis(1, at = time, lwd = 0, labels = format(time, "%b\n%Y"))
  # 
  # plot(time, apply(RsampleMat.CSMC.AS.repM[burnin:niter,], 2, mean),
  #      xlab = "",
  #      type = "l",
  #      ylab = "R",
  #      col="grey",
  #      ylim = c(0,1),
  #      xaxt = "n",
  #      lwd = 2)
  # lines(time, R.CredInterval.LL, col= "darkgrey", lty = 2, lwd = 2)
  # lines(time, R.CredInterval.UL, col= "darkgrey", lty = 2, lwd = 2)
  # axis(1, at = time, lwd = 0, labels = format(time, "%b\n%Y"))
  
  
  ## S plot
  # Compute the x and y coordinates for the shaded area
  xx <- 1:length(S.CredInterval.LL)
  y_lower <- S.CredInterval.LL
  y_upper <- S.CredInterval.UL
  
  # Plot the mean trajectory
  plot(colMeans(SsampleMat.CSMC.AS.repM[burnin:niter,]),
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
  
  # Plot the shaded area for the credible interval
  polygon(c(xx, rev(xx)), c(y_lower, rev(y_upper)), col = rgb(0.5, 0.5, 0.5, alpha = 1/4), border = NA)
  
  
  ## E plot
  # Compute the x and y coordinates for the shaded area
  xx <- 1:length(E.CredInterval.LL)
  y_lower <- E.CredInterval.LL
  y_upper <- E.CredInterval.UL
  
  # Plot the mean trajectory
  plot(colMeans(EsampleMat.CSMC.AS.repM[burnin:niter,]),
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
  plot(colMeans(IsampleMat.CSMC.AS.repM[burnin:niter,]),
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
  plot(colMeans(RsampleMat.CSMC.AS.repM[burnin:niter,]),
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
  
  # Plot the shaded area for the credible interval
  polygon(c(xx, rev(xx)), 
          c(y_lower, rev(y_upper)), 
          col = rgb(0.5, 0.5, 0.5, alpha = 1/4), 
          border = NA)
  
  dev.off()
  
  
  # Posterior I_t
  pdf(paste("(Weekly)TwoRegimeEstimatedIt.pdf"), width = 10, height = 7)
  par(mfrow=c(1,1))
  plot(1:length(y), y, type = "p", col = "grey", pch=20, ylab = "Proportion Infected", xlab = "Time", ylim=c(0,0.05))
  # lines(1:T, theta[3,]*p, col="red")
  lines(1:length(y), apply(IsampleMat.CSMC.AS.repM*mean(parameters.CSMC.AS.repM$p), 2, median), col="grey")
  lines(I.CredInterval.LL*mean(parameters.CSMC.AS.repM$p), col="grey", lty = 2)
  lines(I.CredInterval.UL*mean(parameters.CSMC.AS.repM$p), col="grey", lty = 2)
  dev.off()
  
  
  # Estimated X_t
  pdf(paste("~/Dropbox/(local) Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Real Data/(Weekly)TwoRegimeEstimatedXt.pdf"), width = 10, height = 5)
  par(mfrow=c(1,1),  mar=c(5, 4, 4, 6) + 0.1) 
  
  plot(as.Date(time), y, type = "b", col = "gray50", pch=20, lty=1,
       xlab="",ylab = "Proportion Infected", main="Estimated Regimes",
       ylim=c(0,0.015), xaxt = "n")
  lines(as.Date(time), apply(IsampleMat.CSMC.AS.repM[burnin:niter,]*mean(parameters.CSMC.AS.repM$p), 2, mean), col="gray20")
  # lines(as.Date(time), I.CredInterval.LL*mean(parameters.CSMC.AS.repM$p), col="grey", lty = 2)
  # lines(as.Date(time), I.CredInterval.UL*mean(parameters.CSMC.AS.repM$p), col="grey", lty = 2)
  
  
  # Compute the x and y coordinates for the shaded area
  xx <- as.Date(time)
  y_lower <- I.CredInterval.LL*mean(parameters.CSMC.AS.repM$p)
  y_upper <- I.CredInterval.UL*mean(parameters.CSMC.AS.repM$p)
  
  # Plot the shaded area for the credible interval
  polygon(c(xx, rev(xx)), 
          c(y_lower, rev(y_upper)), 
          col = rgb(0.5, 0.5, 0.5, alpha = 1/4), 
          border = NA)
  
  # Add month labels to x-axis
  axis(1, at = time, lwd = 0, labels = format(time, "%b\n%Y"))
  
  # Obtain the regime at each time point
  est.prob.x <- data.frame("prob.x1" = colSums(XsampleMat.CSMC.AS.repM[burnin:niter,] == 1)/(niter-burnin),
                     "prob.x2" = colSums(XsampleMat.CSMC.AS.repM[burnin:niter,] == 2)/(niter-burnin))
  est.x <- apply(est.prob.x, 1, which.max)
  
  library("R.utils")
  fromto.x1 <- seqToIntervals(which(est.x==1))
  for (i in 1:nrow(fromto.x1)){
    rect(time[fromto.x1[i,1]], -1, time[fromto.x1[i,2]], 1, 
         col = rgb(1,1,1,1/4), border = rgb(1,1,1,1/4))
  }
  fromto.x2 <- seqToIntervals(which(est.x==2))
  for (i in 1:nrow(fromto.x2)){
    rect(time[fromto.x2[i,1]], -1, time[fromto.x2[i,2]], 1, 
         col = rgb(0.9, 0.6, 0.6, 1/3), border = rgb(0.9, 0.6, 0.6, 1/3))
  }  
  
  # Add a vertical line to indicate the new intervention time
  abline(v = as.Date("2020-3-27"), lty = 3, lwd=3, col = red, cex=1.5) # Distancing: Recommendation to avoid gatherings of any size issued
  abline(v = as.Date("2020-4-8"), lty = 3, lwd=2, col = red, cex=1.5) # Closures/Openings: Provincial parks closed
  # abline(v = as.Date("2020-9-8"), lty = 3, lwd=2, col = red, cex=1.5) # Closures/Openings: All nightclubs and stand-alone banquet halls closed: Liquor sales restricted after 10 p.m. at all bars, pubs and restaurants
  abline(v = as.Date("2020-11-19"), lty = 3, lwd=2, col = red) # Distancing: All social gatherings outside household bubbles prohibited
  abline(v = as.Date("2020-12-2"), lty = 3, lwd=2, col = red, cex=1.5) # Closures/Openings: Indoor fitness and team sports prohibited
  abline(v = as.Date("2021-4-5"), lty = 3, lwd=2, col = red, cex=1.5) # Vaccine: B.C. entered Phase 3 of the Immunization Plan. 3rd phase of vaccination announced to target 1) people age 79 to 60; 2) people age 69 to 16 who are clinically extremely vulnerable
  # abline(v = as.Date("2021-9-24"), lty = 3, lwd=2, col = red) # 80% of B.C. residents now fully vaccinated.
  
  # # Add a vertical line to indicate the eased intervention time
  # abline(v = as.Date("2020-6-24"), lty = 3, lwd=1.5, col = blue) # Phase 3 of BC's reopening starts
  # abline(v = as.Date("2021-2-16"), lty = 3, lwd=1.5, col = blue) # As hospitalizations slowly decline, Dr. Henry loosens COVID restrictions, allowing indoor organized events as long as participants are masked and show their vaccine passports. That includes theatres, sporting events and movie theatres. Fitness centres and swimming pools are allowed to go back to pre-pandemic levels. Restaurants and night clubs go back to full capacity with mingling and dancing allowed if masked and vaccinated.
  # abline(v = as.Date("2021-7-1"), lty = 3, lwd=1.5, col = blue) # Most COVID restrictions removed as outdoor gatherings of up to 5,000 people allowed, limits are taken off the number of diners in restaurants but they still cannot socialize between tables, masks no longer mandatory indoors and recreational travel outside the province can resume.
  # 
  
  # legend("topright",
  #        legend=c("Intervention", "Regime Two"),
  #        col = c(red, rgb(0.5, 0.55, 0.8, 1/4)),
  #        lty = c(2,1),
  #        bty = "n",
  #        cex = 0.8)
  
  # plot(time, colSums(XsampleMat.CSMC.AS.repM == 1)/niter, 
  #      type = "l", 
  #      # xlim = c(20, T),
  #      ylim = c(0, 1),            # Adjust accordingly!
  #      xlab = paste("Data record of length T =", T+1),
  #      ylab = "Probability",
  #      main = expression('PG-CSMC-AS-repM'),
  #      col="grey")
  # lines(time, colSums(XsampleMat.CSMC.AS.repM == 2)/niter, col="blue")
  
  # # Add legend
  # legend("topright",
  #        legend=c(expression('Estimated P(X'[t]*'= 1|y'[1:T]*')'), expression('Estimated P(X'[t]*'= 2|y'[1:T]*')')),
  #        col = c("grey", "blue"),
  #        lty = 1,
  #        bty = "n", 
  #        cex = 0.8)
  
  dev.off()
  

  # Acceptance rate
  MCMC.chain.1$acceptance.rate
  MCMC.chain.2$acceptance.rate
  
  # Check time points with regime 1
  time[which(colSums(XsampleMat.CSMC.AS.repM[burnin:niter,] == 1)/(niter-burnin)-colSums(XsampleMat.CSMC.AS.repM[burnin:niter,] == 2)/(niter-burnin)>0)]
  
  # Check time points with regime 2
  time[which(colSums(XsampleMat.CSMC.AS.repM[burnin:niter,] == 2)/(niter-burnin)-colSums(XsampleMat.CSMC.AS.repM[burnin:niter,] == 1)/(niter-burnin)>0)]
  
  # Marginal likelihood
  mean(marginalLogLik.CSMC.AS.repM[burnin:niter])
  sd(marginalLogLik.CSMC.AS.repM[burnin:niter])
  
  ## Timeline of COVID-19 in BC (https://valandos.lt/en/timeline-how-covid-changed-our-lives-over-the-past-three-years)
  # 2020 April 7
  # No new COVID-19 cases reported in the province in 72 hours.
  
  # 2020 April 12
  # B.C. reports hunting licence sales have nearly doubled.
  # 
  # 2020 April 15
  # B.C. extends state of emergency for another two weeks. There have been 75 COVID-19 related deaths.
  # 
  # 2020 April 15
  # Interior Health region records its first COVID-19 death, a man in his 60s.
  # 
  # 2020 April 24
  # Two hundred pairs of a special $399 Dr. Bonnie Henry shoe designed by John Fluevog sell out as soon as they go on sale. Proceeds go to Food Banks B.C.
  
  # 2020 June 1
  # School children return to school using a hybrid model for learning that combines in-class and remote learning.
  # 
  # 2020 June 24
  # Premier John Horgan announces that Phase 3 of B.C.’s reopening can start, meaning people are encouraged to travel within B.C. Hotels, movie theatres, parks and the film industry gradually re-opening. At that time, B.C. had recorded fewer than 3,000 cases.
  # 
  # 2020 July 1
  # Kelowna's Bernard Avenue is closed to traffic so restaurants can expand their patios onto the street.
  # 
  # 2020 July 13
  # After going five weeks with no more than one new COVID-19 case per day in the Interior Health region, new cases start to climb as the result of parties held in Kelowna around the July 1 holiday. That becomes known as the Kelowna Cluster and leads to at least 130 positive tests with hundreds more people going into self-isolation.
  # 
  # 2020 July 27
  # The Kelowna Cluster leads to Dr. Henry ordering a clampdown on the number of people who can be in rental accommodation and each unit is allowed only five visitors.
  # 
  # 2020 July 30
  # A UBCO students puts a magnetic notice board on his truck explaining why he has Washington licence plates following hostile reactions to out of province visitors.
  # 
  # 2020 Aug. 7
  # A woman has her Alberta licence plates stolen while stopping overnight in Kamloops.
  
  # 2020 Aug. 24
  # B.C. case total passes the 5,000 mark.
  
  # 2020 Sept. 8
  # Nightclubs and banquet halls ordered closed.
  
  # 2020 Sept. 25
  # Outbreak at Calvary Chapel in Kelowna leads to seven cases.
  
  # 2020 Oct. 8
  # B.C. passes 10,000 COVID-19 cases.
  
  # 2020 Dec. 7
  # Lock down on socializing in B.C. extended for another month.
  # 
  # 2020 Dec. 8
  # B.C. announces a one-time, tax-free $1,000 per family B.C. Recovery Benefit.
  # 
  # 2020 Dec. 13
  # The organizer of a Kelowna protest against COVID-19 restrictions is fined $2,300 in Kelowna. More fines to come.
  # 
  # 2020 Dec. 14
  # Two employees at Big White Ski Resort test positive for COVID-19. That cluster grows to 203 by Jan. 19, 2021. Most are people who work and live on the mountain. Big White cancels reservations for out of region visitors.
  # 
  # 2020 Dec. 20
  # B.C.’s socializing restrictions over Christmas are some of the toughest in Canada.
  
  # 2021 Jan. 22
  # It’s announced that, in mid-March, people can start registering to get their first vaccine, based on age.
  
  # 2021 March 1
  # B.C.’s vaccination program is in Phase 2, meaning those over 80 living in the community can get first doses. Phase 3 starts in April for those 75-79 and continues by age increments with those aged 60 to 64 expected to get their first doses in June.
  # 
  # 2021 March 11
  # Outdoor gatherings of up to 10 people now allowed.
  # 
  # 2021 March 14
  # As many as 250 rally outside the Harvest Church in Kelowna protesting restrictions on faith services. The church had already received several $2,300 fines for holding services contrary to public health orders.
  # 
  # 2021 March 23
  # Ban on outdoor faith services lifted.
  # 
  # 2021 March 25
  # Dr. Henry announces that indoor church services can be held from March 28 to May 13, with safety plans in place.
  # 
  # 2021 March 29
  # Indoor church services banned, as are gyms, fitness clubs and indoor dining in restaurants. Students in Grades 4 to 12 have to wear masks as COVID case numbers spike.
  # 
  # 2021 April 1
  # Visitors allowed back into long term care homes.
  # 
  # 2021 April 5
  # Houseboat party in Kelowna shut down as parties abound in the city of the Easter weekend.
  # 
  # 2021 April 8
  # Record high single daily count of new COVID cases in B.C. at 1,293.
  # 
  # 2021 April 10
  # Pharmacies start providing AstraZeneca vaccines.
  # 
  # 2021 April 12
  # Another downtown Kelowna rally, this time in favour of reopening churches, draws 200 protesters and no fines.
  # 
  # Those 40 and over now eligible for first vaccine shots.
  # 
  # 2021 April 19
  # Those over the age of 18 can start registering for vaccines.
  
  # 2021 June 15
  # Further easing of COVID restrictions.
  
  # 2021 July 1
  # Most COVID restrictions removed as outdoor gatherings of up to 5,000 people allowed, limits are taken off 
  # the number of diners in restaurants but they still cannot socialize between tables, masks no longer mandatory 
  # indoors and recreational travel outside the province can resume.
  
  # 2021 July 1
  # Most COVID restrictions removed as outdoor gatherings of up to 5,000 people allowed, limits are taken off the number of diners in restaurants but they still cannot socialize between tables, masks no longer mandatory indoors and recreational travel outside the province can resume.
  
  # 2021 Sept. 24
  # 80% of B.C. residents now fully vaccinated.
  # 
  # 2021 Sept. 28
  # “Pandemic of the unvaccinated” is the new catchphrase but older, vaccinated people are disproportionately dying in the pandemic.
  # 
  # 2021 Oct. 4
  # All students in schools have to wear masks.
  # 
  # 2021 Oct. 5
  # Provincial government announces all its employees must be fully vaccinated by Nov. 22 as will all visitors to long term care and assisted living homes by Oct. 12. Visitors to acute care hospitals have until Oct. 26 to be vaccinated.
  
  # 2021 Oct. 19
  # Sports events, weddings and concerts return to 100% capacity.
  
  # 2021 Nov. 26
  # Omicron variant sparks travel ban to Canada from a select few countries. 
  # Omicron makes it to Canada anyways.
  
  
  # 2021 Nov. 30
  # International flights resume at Kelowna International Airport, long after some smaller Canadian airports got that right back after the global lockdown of air travel in 2020.
  # 
  # Stricter lockdown rules relaxed for the Interior Health region, such as indoor venues going to 100% occupancy.
  # 
  # 2021 Dec. 7
  # First five cases of Omicron appear in B.C.
  # 
  # 2021 Dec. 20
  # Rising Omicron case counts trigger new restrictions, including the cancellation of New Year’s Eve celebrations, limiting personal indoor gatherings to 10 people, venues with 1,000 or more seats limited to 50% capacity and everyone attending indoor events must be vaccinated. All sports tournaments and travel associated with them are cancelled.
  
  
  