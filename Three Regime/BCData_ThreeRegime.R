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

ggsave(file="~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Real Data/BCWeeklyActiveCases.png", g.weekly)

# scale_x_date(limit=c(as.Date("2020-07-20"),as.Date("2021-07-02"))) # extract certain period

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
regimes <- c(1, 2, 3)
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
  a.p = 0.25,
  b.p = 0.25,
  delta.mat = matrix(c(10, 1, 1, 1, 10, 1, 1, 1, 10), nrow = lenXset, ncol = lenXset),
  a.f = c(0.5, 0),
  b.f = c(1, 0.5)) 

# Call PG sampler
source("~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Code/Three Regime/ParticleGibbs_ThreeRegimes.R")

ptm <- proc.time() 
PG.results <- foreach(i = 1:nchain, .combine = "list", .packages = c("truncnorm", "DirichletReg")) %dopar% {
  # Run PG sampler
  PG.CSMC.AS(y, regimes, M, niter, hyperparams, pop.size=1)
  
} 
proc.time()-ptm 

stopCluster(cl)


# Save results to local file
setwd("~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Output/Real Data")
saveRDS(hyperparams, paste0("PG_results_hyperparams", "_niter", niter,"_M", M, "_K", length(regimes), ".rds"))
saveRDS(PG.results, paste0("PG_results", "_niter", niter, "_M", M, "_K", length(regimes), ".rds"))

# Read results
# hyperparams <- readRDS("~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Output/Real Data/PG_results_hyperparams_BCWeekly_niter11000_M50_K3.rds")
PG.results <- readRDS("~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Output/Real Data/PG_results_BCWeekly_niter11000_M50_K3.rds")

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
setwd("~/Dropbox/Beta-Dirichlet-Time-Series-Model/pMCMC - BDSSSM/Figures/Real Data")

# Trace plot of parameters
pdf(paste("ThreeRegimeTracePlot.pdf"), width = 10, height = 20)
par(mfrow=c(5,2))
plot(parameters.CSMC.AS.repM$alpha[1, 1:niter], type="l", xlab="Iterations", ylab = expression(alpha), cex=2)
plot(parameters.CSMC.AS.repM$beta[1, 1:niter], type="l", xlab="Iterations", ylab = expression(beta), cex=2)
plot(parameters.CSMC.AS.repM$gamma[1, 1:niter], type="l", xlab="Iterations", ylab = expression(gamma), cex=2)
plot(parameters.CSMC.AS.repM$kappa[1, 1:niter], type="l", xlab="Iterations", ylab = expression(kappa), cex=2)
plot(parameters.CSMC.AS.repM$lambda[1, 1:niter], type="l", xlab="Iterations", ylab = expression(lambda), cex=2)
# plot(parameters.CSMC.AS.repM$p[1, 1:niter], type="l", xlab="Iterations", ylab = expression(p))
# plot(parameters.CSMC.AS.repM$f[1, 1:niter], type="l", xlab="Iterations", ylab = expression(f[1]))
plot(parameters.CSMC.AS.repM$f[2, 1:niter], type="l", xlab="Iterations", ylab = expression(f[2]), cex=2)
plot(parameters.CSMC.AS.repM$f[3, 1:niter], type="l", xlab="Iterations", ylab = expression(f[3]))


# Posterior Px
post.pi.k1 <- matrix(0, niter, lenXset, byrow=TRUE)
post.pi.k2 <- matrix(0, niter, lenXset, byrow=TRUE)
post.pi.k3 <- matrix(0, niter, lenXset, byrow=TRUE)
for (i in 1:(niter)){
  post.pi.k1[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][1,]
  post.pi.k2[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][2,]
  post.pi.k3[i-1,] <- parameters.CSMC.AS.repM$Px[[i]][3,]
}
post.pi.k1
post.pi.k2
post.pi.k3


# pi_{k1}
plot(post.pi.k1[,1], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[11]), cex=2) # pi_11
plot(post.pi.k1[,2], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[12]), cex=2) # pi_12
# plot(post.pi.k1[,3],
#      type="l",
#      xlab="Iterations",
#      ylab = expression(pi[13]), cex=2) # pi_13


# pi_{k2}
plot(post.pi.k2[,1], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[21]), cex=2) # pi_21
plot(post.pi.k2[,2], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[22]), cex=2) # pi_22
# plot(post.pi.k2[,3],
#      type="l",
#      xlab="Iterations",
#      ylab = expression(pi[23]), cex=2) # pi_23

# pi_{k3}
plot(post.pi.k3[,1], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[31]), cex=2) # pi_31
plot(post.pi.k3[,2], 
     type="l", 
     xlab="Iterations", 
     ylab = expression(pi[32]), cex=2) # pi_32
# plot(post.pi.k3[,3],
#      type="l",
#      xlab="Iterations",
#      ylab = expression(pi[33]), cex=2) # pi_33

dev.off()


# Histogram of parameters
pdf(paste("ThreeRegimeHistogram.pdf"), width = 10, height = 20)
par(mfrow=c(5,3))
hist(parameters.CSMC.AS.repM$alpha[1:niter], xlab="", main = expression(alpha))
abline(v=mean(parameters.CSMC.AS.repM$alpha[1:niter]), col="red", lwd=2)
hist(parameters.CSMC.AS.repM$beta[1:niter], xlab="",main = expression(beta))
abline(v=mean(parameters.CSMC.AS.repM$beta[1:niter]), col="red", lwd=2)
hist(parameters.CSMC.AS.repM$gamma[1:niter], xlab="",main = expression(gamma))
abline(v=mean(parameters.CSMC.AS.repM$gamma[1:niter]), col="red", lwd=2)
hist(parameters.CSMC.AS.repM$kappa[1:niter], xlab="",main = expression(kappa))
abline(v=mean(parameters.CSMC.AS.repM$kappa[1:niter]), col="red", lwd=2)
hist(parameters.CSMC.AS.repM$lambda[1:niter], xlab="",main = expression(lambda))
abline(v=mean(parameters.CSMC.AS.repM$lambda[1:niter]), col="red", lwd=2)
# hist(parameters.CSMC.AS.repM$p[1:niter], xlab="",main = "p")
# hist(parameters.CSMC.AS.repM$f[1,1:niter], xlab="",main = expression(f[1]))
hist(parameters.CSMC.AS.repM$f[2,1:niter], xlab="",main = expression(f[2]))
abline(v=mean(parameters.CSMC.AS.repM$f[2, 1:niter]), col="red", lwd=2)
hist(parameters.CSMC.AS.repM$f[3,1:niter], xlab="",main = expression(f[3]))
abline(v=mean(parameters.CSMC.AS.repM$f[3, 1:niter]), col="red", lwd=2)
hist(post.pi.k1[,1], xlab="",main = expression(pi[11]))
abline(v=mean(post.pi.k1[,1]), col="red", lwd=2)
hist(post.pi.k1[,2], xlab="",main = expression(pi[12]))
abline(v=mean(post.pi.k1[,2]), col="red", lwd=2)
# hist(post.pi.k1[,3], xlab="",main = expression(pi[13]))
# abline(v=mean(post.pi.k1[,3]), col="red", lwd=2)
hist(post.pi.k2[,1], xlab="",main = expression(pi[21]))
abline(v=mean(post.pi.k2[,1]), col="red", lwd=2)
hist(post.pi.k2[,2], xlab="",main = expression(pi[22]))
abline(v=mean(post.pi.k2[,2]), col="red", lwd=2)
# hist(post.pi.k2[,3], xlab="",main = expression(pi[23]))
# abline(v=mean(post.pi.k2[,3]), col="red", lwd=2)
hist(post.pi.k3[,1], xlab="",main = expression(pi[31]))
abline(v=mean(post.pi.k3[,1]), col="red", lwd=2)
hist(post.pi.k3[,2], xlab="",main = expression(pi[32]))
abline(v=mean(post.pi.k3[,2]), col="red", lwd=2)
# hist(post.pi.k3[,3], xlab="",main = expression(pi[33]))
# abline(v=mean(post.pi.k3[,3]), col="red", lwd=2)
dev.off()



## 95% Credible Interval of posterior samples
S.CredInterval.LL <- c()
S.CredInterval.UL <- c()
E.CredInterval.LL <- c()
E.CredInterval.UL <- c()
I.CredInterval.LL <- c()
I.CredInterval.UL <- c()
R.CredInterval.LL <- c()
R.CredInterval.UL <- c()

for (i in 1:T){
  S.CredInterval.LL <- c(S.CredInterval.LL, quantile(SsampleMat.CSMC.AS.repM[,i], 0.01))
  S.CredInterval.UL <- c(S.CredInterval.UL, quantile(SsampleMat.CSMC.AS.repM[,i], 0.99))
  E.CredInterval.LL <- c(E.CredInterval.LL, quantile(EsampleMat.CSMC.AS.repM[,i], 0.01))
  E.CredInterval.UL <- c(E.CredInterval.UL, quantile(EsampleMat.CSMC.AS.repM[,i], 0.99))
  I.CredInterval.LL <- c(I.CredInterval.LL, quantile(IsampleMat.CSMC.AS.repM[,i], 0.01))
  I.CredInterval.UL <- c(I.CredInterval.UL, quantile(IsampleMat.CSMC.AS.repM[,i], 0.99))
  R.CredInterval.LL <- c(R.CredInterval.LL, quantile(RsampleMat.CSMC.AS.repM[,i], 0.01))
  R.CredInterval.UL <- c(R.CredInterval.UL, quantile(RsampleMat.CSMC.AS.repM[,i], 0.99))
}


setwd("/Users/Grace/Dropbox/Bayesian-Inference-on-SSSM/BC Cases")
pdf(paste(Sys.Date(), "ThreeRegimeEstimatedSEIR.pdf"), width = 10, height = 10)

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

dev.off()


# Posterior I_t
pdf(paste(Sys.Date(), "ThreeRegimeEstimatedIt.pdf"), width = 10, height = 7)
par(mfrow=c(1,1))
plot(1:length(y), y, type = "p", col = "grey", pch=20, ylab = "Proportion Infected", ylim=c(0,0.05))
# lines(1:T, theta[3,]*p, col="red")
lines(1:(T+1), apply(IsampleMat.CSMC.AS.repM*mean(parameters.CSMC.AS.repM$p), 2, mean), col="grey")
lines(I.CredInterval.LL*mean(parameters.CSMC.AS.repM$p), col="grey", lty = 2)
lines(I.CredInterval.UL*mean(parameters.CSMC.AS.repM$p), col="grey", lty = 2)
dev.off()

# Estimated X_t
pdf(paste(Sys.Date(), "ThreeRegimeEstimatedXt.pdf"), width = 10, height = 7)
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

# Add legend
legend("topright",
       legend=c("k=1", "k=2", "k=3"),
       col = c("grey", "blue", "darkblue"),
       lty = 1,
       bty = "n", 
       cex = 0.8)
dev.off()



## Timeline of COVID-19 in BC (https://valandos.lt/en/timeline-how-covid-changed-our-lives-over-the-past-three-years)
# 2020 June 1
# School children return to school using a hybrid model for learning that combines in-class and remote learning.

# 2020 June 24
# Premier John Horgan announces that Phase 3 of B.C.’s reopening can start, meaning people are encouraged to 
# travel within B.C. Hotels, movie theatres, parks and the film industry gradually re-opening. 
# At that time, B.C. had recorded fewer than 3,000 cases.

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
# 2021Oct. 5
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


