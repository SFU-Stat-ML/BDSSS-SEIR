
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

bc_data <- read.csv(file = "../../Data/Real Data/BC COVID CASES - Daily Cases.csv",
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

############################ Run PG Sampler ##############################

library(foreach)
library(doParallel)
detectCores()
cl <- makeCluster(6)
registerDoParallel(cl)


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
source("ParticleGibbs_TwoRegime.R")

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

str(PG.results)

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

### Trace plot of parameters
par(mfrow=c(3,3))
plot(parameters.CSMC.AS.repM$alpha[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(alpha), cex=2)
plot(parameters.CSMC.AS.repM$beta[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(beta), cex=2)
plot(parameters.CSMC.AS.repM$gamma[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(gamma), cex=2)
plot(parameters.CSMC.AS.repM$kappa[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(kappa), cex=2)
plot(parameters.CSMC.AS.repM$lambda[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(lambda), cex=2)
plot(parameters.CSMC.AS.repM$p[1, burnin:niter], type="l", xlab="Iterations", ylab = expression(p))
plot(parameters.CSMC.AS.repM$f[2, burnin:niter], type="l", xlab="Iterations", ylab = expression(f[2]), cex=2)


# Posterior Px
post.pi.k1 <- matrix(0, niter, lenXset, byrow=TRUE)
post.pi.k2 <- matrix(0, niter, lenXset, byrow=TRUE)
for (i in 1:(niter)){
  post.pi.k1[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][1,]
  post.pi.k2[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][2,]
}
post.pi.k1
post.pi.k2

# pi_{11}
plot(post.pi.k1[burnin:(niter-1),1], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[11]), cex=2) 

# pi_{22}
plot(post.pi.k2[burnin:(niter-1),2], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[22]), cex=2)




### Histogram of parameters
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



# Estimated X_t
par(mfrow=c(1,1),  mar=c(5, 4, 4, 6) + 0.1) 

plot(as.Date(time), y, type = "b", col = "gray50", pch=20, lty=1,
     xlab="",ylab = "Proportion Infected", main="Estimated Regimes",
     ylim=c(0,0.015), xaxt = "n")
lines(as.Date(time), apply(IsampleMat.CSMC.AS.repM[burnin:niter,]*mean(parameters.CSMC.AS.repM$p[burnin:niter]), 2, mean), col="gray20")

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
abline(v = as.Date("2020-11-19"), lty = 3, lwd=2, col = red) # Distancing: All social gatherings outside household bubbles prohibited
abline(v = as.Date("2020-12-2"), lty = 3, lwd=2, col = red, cex=1.5) # Closures/Openings: Indoor fitness and team sports prohibited
abline(v = as.Date("2021-4-5"), lty = 3, lwd=2, col = red, cex=1.5) # Vaccine: B.C. entered Phase 3 of the Immunization Plan. 3rd phase of vaccination announced to target 1) people age 79 to 60; 2) people age 69 to 16 who are clinically extremely vulnerable


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

