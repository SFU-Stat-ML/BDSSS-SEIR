################################################################################
# Author: Jingxue (Grace) Feng
#         Simon Fraser University, Burnaby, BC, Canada
#         Email: jingxuef@sfu.ca 
################################################################################

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
  geom_line( color="steelblue") + 
  geom_vline(xintercept = as.numeric(as.Date("2021-11-26")), 
             linetype = "dashed", color = "red", size = 1) +
  geom_point() +
  xlab("") +
  ylab("Weekly Active Cases in B.C") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) + 
  scale_x_date(limit = c(as.Date(min(bc_data_active$Date)), as.Date(max(bc_data_active$Date)))) +
  annotate("text", x = as.Date("2021-11-30"), y = max(bc_data_active$ActiveCases), 
           label = "1st case of the COVID-19 Omicron variant", vjust = -30, hjust = 1.1,
           size = 5, color = "red")
g.weekly

############################ Run PG Sampler ##############################

library(foreach)
library(doParallel)
detectCores()
cl <- makeCluster(2)
registerDoParallel(cl)

## Daily active cases (before omicron starts)
# y <- bc_data_active$ActiveCases[1:438]/5070000
# time <- bc_data_active$Date[1:438]

## Weekly active cases
y <- new_data_weekly$weekly.n[1:95]/5070000
time <- new_data_weekly$Week[1:95]

## Regimes
regimes <- c(1, 2, 3, 4)
lenXset <- length(regimes)

## Number of MCMC chains
nchain  <- 2    

## Number of MCMC iterations
burnin <- 1000
niter   <- 10000 + burnin      # number of MCMC iterations

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
  delta.mat = matrix(c(10, 1, 1, 1,
                       1, 10, 1, 1, 
                       1, 1, 10, 1,
                       1, 1, 1, 10), nrow = lenXset, ncol = lenXset),
  a.f = c(2/3, 1/3, 0),
  b.f = c(1, 2/3, 1/3)) 


# Call PG sampler
source("ParticleGibbs_FourRegimes.R")

ptm <- proc.time() 
PG.results <- foreach(i = 1:nchain, .combine = "list", .packages = c("truncnorm", "DirichletReg")) %dopar% {
  # Run PG sampler
  PG.CSMC.AS(y, regimes, M, niter, hyperparams, pop.size=1)
  
} 
proc.time()-ptm 

stopCluster(cl)

# 
# # Save results to local file
# setwd("~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Output/Real Data")
# saveRDS(hyperparams, paste0("PG_results_hyperparams_BCWeekly_niter", niter,"_M", M, "_K", length(regimes), ".rds"))
# saveRDS(PG.results, paste0("PG_results_BCWeekly_niter", niter, "_M", M, "_K", length(regimes), ".rds"))
# 
# # Read results
# # hyperparams <- readRDS("~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Output/Real Data/PG_results_hyperparams_BCWeekly_niter11000_M50_K3.rds")
# PG.results <- readRDS("~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Output/Real Data/PG_results_BCWeekly_niter11000_M50_K3.rds")

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



############################## Data Visualization ##################################

# Trace plot of parameters
par(mfrow=c(5,2))
plot(parameters.CSMC.AS.repM$alpha[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(alpha), cex=2)
plot(parameters.CSMC.AS.repM$beta[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(beta), cex=2)
plot(parameters.CSMC.AS.repM$gamma[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(gamma), cex=2)
plot(parameters.CSMC.AS.repM$kappa[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(kappa), cex=2)
plot(parameters.CSMC.AS.repM$lambda[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(lambda), cex=2)
plot(parameters.CSMC.AS.repM$p[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(p))
plot(parameters.CSMC.AS.repM$f[2, burnin:niter], type="l", xlab="Iterations", ylab = expression(f[2]), cex=2)
plot(parameters.CSMC.AS.repM$f[3, burnin:niter], type="l", xlab="Iterations", ylab = expression(f[3]))
plot(parameters.CSMC.AS.repM$f[4, burnin:niter], type="l", xlab="Iterations", ylab = expression(f[3]))



# Posterior Px
post.pi.k1 <- matrix(0, niter, lenXset, byrow=TRUE)
post.pi.k2 <- matrix(0, niter, lenXset, byrow=TRUE)
post.pi.k3 <- matrix(0, niter, lenXset, byrow=TRUE)
post.pi.k4 <- matrix(0, niter, lenXset, byrow=TRUE)
for (i in burnin:niter){
  post.pi.k1[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][1,]
  post.pi.k2[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][2,]
  post.pi.k3[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][3,]
  post.pi.k4[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][4,]
}
post.pi.k1
post.pi.k2
post.pi.k3
post.pi.k4

# pi_{k1}
plot(post.pi.k1[,1], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[11]), cex=2) # pi_11
plot(post.pi.k1[,2], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[12]), cex=2) # pi_12
plot(post.pi.k1[,3],
     type="l",
     xlab="Iterations",
     ylab = expression(pi[13]), cex=2) # pi_13


# pi_{k2}
plot(post.pi.k2[,1], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[21]), cex=2) # pi_21
plot(post.pi.k2[,2], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[22]), cex=2) # pi_22
plot(post.pi.k2[,3],
     type="l",
     xlab="Iterations",
     ylab = expression(pi[23]), cex=2) # pi_23

# pi_{k3}
plot(post.pi.k3[,1], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[31]), cex=2) # pi_31
plot(post.pi.k3[,2], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[32]), cex=2) # pi_32
plot(post.pi.k3[,3],
     type="l",
     xlab="Iterations",
     ylab = expression(pi[33]), cex=2) # pi_33

# pi_{k4}
plot(post.pi.k4[,1], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[41]), cex=2) # pi_31
plot(post.pi.k4[,2], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[42]), cex=2) # pi_32
plot(post.pi.k4[,3],
     type="l",
     xlab="Iterations",
     ylab = expression(pi[43]), cex=2) # pi_33




# Histogram of parameters
par(mfrow=c(5,3))
hist(parameters.CSMC.AS.repM$alpha[burnin:niter], xlab="", main = expression(alpha))
abline(v=mean(parameters.CSMC.AS.repM$alpha[burnin:niter]), col="red", lwd=2)
hist(parameters.CSMC.AS.repM$beta[burnin:niter], xlab="",main = expression(beta))
abline(v=mean(parameters.CSMC.AS.repM$beta[burnin:niter]), col="red", lwd=2)
hist(parameters.CSMC.AS.repM$gamma[burnin:niter], xlab="",main = expression(gamma))
abline(v=mean(parameters.CSMC.AS.repM$gamma[burnin:niter]), col="red", lwd=2)
hist(parameters.CSMC.AS.repM$kappa[burnin:niter], xlab="",main = expression(kappa))
abline(v=mean(parameters.CSMC.AS.repM$kappa[burnin:niter]), col="red", lwd=2)
hist(parameters.CSMC.AS.repM$lambda[burnin:niter], xlab="",main = expression(lambda))
abline(v=mean(parameters.CSMC.AS.repM$lambda[burnin:niter]), col="red", lwd=2)
hist(parameters.CSMC.AS.repM$f[2,burnin:niter], xlab="",main = expression(f[2]))
abline(v=mean(parameters.CSMC.AS.repM$f[2, burnin:niter]), col="red", lwd=2)
hist(parameters.CSMC.AS.repM$f[3,burnin:niter], xlab="",main = expression(f[3]))
abline(v=mean(parameters.CSMC.AS.repM$f[3, burnin:niter]), col="red", lwd=2)
hist(post.pi.k1[,1], xlab="",main = expression(pi[11]))
abline(v=mean(post.pi.k1[,1]), col="red", lwd=2)
hist(post.pi.k1[,2], xlab="",main = expression(pi[12]))
abline(v=mean(post.pi.k1[,2]), col="red", lwd=2)
hist(post.pi.k1[,3], xlab="",main = expression(pi[13]))
abline(v=mean(post.pi.k1[,3]), col="red", lwd=2)
hist(post.pi.k2[,1], xlab="",main = expression(pi[21]))
abline(v=mean(post.pi.k2[,1]), col="red", lwd=2)
hist(post.pi.k2[,2], xlab="",main = expression(pi[22]))
abline(v=mean(post.pi.k2[,2]), col="red", lwd=2)
hist(post.pi.k2[,3], xlab="",main = expression(pi[23]))
abline(v=mean(post.pi.k2[,3]), col="red", lwd=2)
hist(post.pi.k3[,1], xlab="",main = expression(pi[31]))
abline(v=mean(post.pi.k3[,1]), col="red", lwd=2)
hist(post.pi.k3[,2], xlab="",main = expression(pi[32]))
abline(v=mean(post.pi.k3[,2]), col="red", lwd=2)
hist(post.pi.k3[,3], xlab="",main = expression(pi[33]))
abline(v=mean(post.pi.k3[,3]), col="red", lwd=2)



## 95% Credible Interval of posterior samples
S.CredInterval <- t(apply(SsampleMat.CSMC.AS.repM[burnin:iter,], 2, function(x) quantile(x, c(0.025, 0.975))))
S.CredInterval.LL <- S.CredInterval[,1]
S.CredInterval.UL <- S.CredInterval[,2]
E.CredInterval <- t(apply(EsampleMat.CSMC.AS.repM[burnin:iter,], 2, function(x) quantile(x, c(0.025, 0.975))))
E.CredInterval.LL <- E.CredInterval[,1]
E.CredInterval.UL <- E.CredInterval[,2]
I.CredInterval <- t(apply(IsampleMat.CSMC.AS.repM[burnin:iter,], 2, function(x) quantile(x, c(0.025, 0.975))))
I.CredInterval.LL <- I.CredInterval[,1]
I.CredInterval.UL <- I.CredInterval[,2]
R.CredInterval <- t(apply(RsampleMat.CSMC.AS.repM[burnin:iter,], 2, function(x) quantile(x, c(0.025, 0.975))))
R.CredInterval.LL <- R.CredInterval[,1]
R.CredInterval.UL <- R.CredInterval[,2]

par(mfrow=c(2,2))

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
     ylim = c(0, 0.2),
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


# Add legend
legend("topright",
       legend = c("Posterior Mean"),
       lty = c(2),
       lwd = c(1.5),
       col = c("#969696"),
       bty = "n",
       cex = 1.2)





# Estimated X_t
par(mfrow=c(2,1))
plot(time, y, type = "b", col = "grey", pch=20, lty=1,
     ylab = "Proportion Infected", ylim=c(0,0.05))
plot(time, colSums(XsampleMat.CSMC.AS.repM == 1)/niter, 
     type = "l", 
     # xlim = c(20, T),
     ylim = c(0, 1),            # Adjust accordingly!
     xlab = paste("Data record of length T =", T+1),
     ylab = expression('Estimated P(X'[t]*'= k|y'[1:T]*')'),
     main = expression('PG-CSMC-AS-repM'),
     col="grey")
lines(time, colSums(XsampleMat.CSMC.AS.repM == 2)/niter, col="blue")
lines(time, colSums(XsampleMat.CSMC.AS.repM == 3)/niter, col="darkblue")
lines(time, colSums(XsampleMat.CSMC.AS.repM == 4)/niter, col="black")

# Add legend
legend("topright",
       legend=c("k=1", "k=2", "k=3", "k=4"),
       col = c("grey", "blue", "darkblue", "black"),
       lty = 1,
       bty = "n", 
       cex = 0.8)



## Marginal likelihood of the model
mean(marginalLogLik.CSMC.AS.repM[burnin:iter])
sd(marginalLogLik.CSMC.AS.repM[burnin:iter])