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

bc_data <- read.csv(file = "../../../Data/Real Data/BC COVID CASES - Daily Cases.csv",
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

# ggsave(file="~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Real Data/BCDailyActiveCases.png", g)


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

# ggsave(file="~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Real Data/BCWeeklyActiveCases.png", g.weekly)

# scale_x_date(limit=c(as.Date("2020-07-20"),as.Date("2021-07-02"))) # extract certain period

############################ Run PG Sampler, Part I#######################
####### Time period: the week of 2020-01-27 ~ 2020-10-05 #################

library(foreach)
library(doParallel)
detectCores()
cl <- makeCluster(3)
registerDoParallel(cl)

## Daily active cases (before omicron starts)
# y <- bc_data_active$ActiveCases[1:438]/5070000
# time <- bc_data_active$Date[1:438]

## Weekly active cases
y <- new_data_weekly$weekly.n[1:95]/5070000
time <- new_data_weekly$Week[1:95]

## Regimes
regimes <- c(1)

## Number of MCMC iterations
burnin <- 1000
niter   <- 10000 + burnin      # number of MCMC iterations

## Number of theta_t for each x_t 
M <- 50

## Set up hyperparameters
hyperparams <- list(
  m.alpha = 0.2*7,
  sigma.alpha = 0.5,
  m.beta = 0.2*7,
  sigma.beta = 0.5,
  m.gamma = 0.2*7,
  sigma.gamma = 0.5,
  a.lambda = 20,
  b.lambda = 0.01,
  a.kappa = 200,
  b.kappa = 0.01,
  a.p1 = 0.1,
  b.p1 = 0.4,
  m.p1 = 0.2,
  sigma.p1 = 0.05,
  a.p2 = 0.2,
  b.p2 = 1,
  m.p2 = 0.3,
  sigma.p2 = 0.05,
  delta.mat = matrix(c(1), nrow = 1, ncol = 1),
  a.f = 1,
  b.f = 1) 


# Time of introducing antigen test
T_star <- 37

# Call PG sampler
source("ParticleGibbs_OneRegime.R")

# Run one MCMC chain
set.seed(123)
ptm <- proc.time()
PG.results <- PG.CSMC.AS(y, regimes, M, niter, hyperparams, pop.size=1, num.mh.updates=5, T_star)
proc.time()-ptm

# PG.results$acceptance.rate

# # Extract results
SsampleMat.CSMC.AS.repM <- PG.results$SsampleMat.CSMC.AS.repM
EsampleMat.CSMC.AS.repM <- PG.results$EsampleMat.CSMC.AS.repM
IsampleMat.CSMC.AS.repM <- PG.results$IsampleMat.CSMC.AS.repM
RsampleMat.CSMC.AS.repM <- PG.results$RsampleMat.CSMC.AS.repM
XsampleMat.CSMC.AS.repM <- PG.results$XsampleMat.CSMC.AS.repM
parameters.CSMC.AS.repM <- PG.results$parameters.CSMC.AS.repM
marginalLogLik.CSMC.AS.repM <- PG.results$marginalLogLik.CSMC.AS.repM

############################ Data Visualization #######################################

# Trace plot of parameters
par(mfrow=c(3,3))
plot(parameters.CSMC.AS.repM$alpha[1, ], type="l", xlab="Iterations", ylab = expression(alpha))
plot(parameters.CSMC.AS.repM$beta[1, ], type="l", xlab="Iterations", ylab = expression(beta))
plot(parameters.CSMC.AS.repM$gamma[1, ], type="l", xlab="Iterations", ylab = expression(gamma))
plot(parameters.CSMC.AS.repM$kappa[1, ], type="l", xlab="Iterations", ylab = expression(kappa))
plot(parameters.CSMC.AS.repM$lambda[1, ], type="l", xlab="Iterations", ylab = expression(lambda))
plot(parameters.CSMC.AS.repM$p1[1, ], type="l", xlab="Iterations", ylab = expression(p[1]))
plot(parameters.CSMC.AS.repM$p2[1, ], type="l", xlab="Iterations", ylab = expression(p[2]))

# Histogram of parameters
par(mfrow=c(3,3))
hist(parameters.CSMC.AS.repM$alpha[1:niter], xlab="", main = expression(alpha))
hist(parameters.CSMC.AS.repM$beta[1:niter], xlab="",main = expression(beta))
hist(parameters.CSMC.AS.repM$gamma[1:niter], xlab="",main = expression(gamma))
hist(parameters.CSMC.AS.repM$kappa[1:niter], xlab="",main = expression(kappa))
hist(parameters.CSMC.AS.repM$lambda[1:niter], xlab="",main = expression(lambda))

# Estimated S_t, E_t, I_t, R_t
## 95% Credible Interval of posterior samples
S.CredInterval <- t(apply(SsampleMat.CSMC.AS.repM, 2, function(x) quantile(x, c(0.025, 0.975))))
S.CredInterval.LL <- S.CredInterval[,1]
S.CredInterval.UL <- S.CredInterval[,2]
E.CredInterval <- t(apply(EsampleMat.CSMC.AS.repM, 2, function(x) quantile(x, c(0.025, 0.975))))
E.CredInterval.LL <- E.CredInterval[,1]
E.CredInterval.UL <- E.CredInterval[,2]
I.CredInterval <- t(apply(IsampleMat.CSMC.AS.repM, 2, function(x) quantile(x, c(0.025, 0.975))))
I.CredInterval.LL <- I.CredInterval[,1]
I.CredInterval.UL <- I.CredInterval[,2]
R.CredInterval <- t(apply(RsampleMat.CSMC.AS.repM, 2, function(x) quantile(x, c(0.025, 0.975))))
R.CredInterval.LL <- R.CredInterval[,1]
R.CredInterval.UL <- R.CredInterval[,2]


par(mfrow=c(2,2))
plot(colMeans(SsampleMat.CSMC.AS.repM),
     xlab="Time",
     type = "l",
     ylab = "S",
     col="grey",
     ylim=c(0,1))
lines(S.CredInterval.LL, col= "grey", lty = 2)
lines(S.CredInterval.UL, col= "grey", lty = 2)

plot(colMeans(EsampleMat.CSMC.AS.repM),
     xlab="Time",
     type = "l",
     ylab = "E",
     col="grey",
     ylim=c(0,0.1))
lines(E.CredInterval.LL, col= "grey", lty = 2)
lines(E.CredInterval.UL, col= "grey", lty = 2)

plot(colMeans(IsampleMat.CSMC.AS.repM),
     xlab="Time",
     type = "l",
     ylab = "I",
     col="grey",
     ylim=c(0,0.1))
lines(I.CredInterval.LL, col= "grey", lty = 2)
lines(I.CredInterval.UL, col= "grey", lty = 2)

plot(colMeans(RsampleMat.CSMC.AS.repM),
     xlab="Time",
     type = "l",
     ylab = "R",
     col="grey",
     ylim=c(0,1))
lines(R.CredInterval.LL, col= "grey", lty = 2)
lines(R.CredInterval.UL, col= "grey", lty = 2)


# Posterior I_t
par(mfrow=c(1,1))
plot(1:length(y), y, type = "p", col = "grey", pch=20, ylab = "Proportion Infected", ylim=c(0,0.03))
#lines(1:(T+1), apply(IsampleMat.CSMC.AS.repM*mean(parameters.CSMC.AS.repM$p), 2, median), col="grey")
#lines(I.CredInterval.LL*mean(parameters.CSMC.AS.repM$p), col="grey", lty = 2)
#lines(I.CredInterval.UL*mean(parameters.CSMC.AS.repM$p), col="grey", lty = 2)

# Extimated X_t
plot(colSums(XsampleMat.CSMC.AS.repM == 1)/niter, 
     type = "l", 
     xlim = c(0, T),
     ylim = c(0, 1),            # Adjust accordingly!
     xlab = paste("Data record of length T =", T),
     ylab = expression('Estimated P(X'[t]*'= 1|y'[1:T]*')'),
     main = expression('PG-CSMC-AS with replicator M'))

# Marginal log likelihood
mean(marginalLogLik.CSMC.AS.repM[burnin:niter])
plot(burnin:niter, marginalLogLik.CSMC.AS.repM[burnin:niter])

