rm(list = ls())
setwd("E:/RR/Save2022_10_14/Fan1")
t1 <- Sys.time()
p <- 2000; n <- 600 ; K <- 3; lb <- 0.1; ub <- 0.4; v <- 4
# p <- 2000; n <- 50 ; K <- 3; lb <- 0.3; ub <- 0.6; v <- 4
PP <- 0.2
gama<-0.05
source("use.R");source("core.R");library(MASS)
iteraton <- 100
meann <- rep(0,p)
hh <- ttrue_pai1 <- ppai1_hat <- FFDP <- ttrue_FDP <- NNDR <- FFNR <- Ft <- St <- RR <- Recive <- ppower2 <- rep(0,iteraton)
for (jj in 1:iteraton){
  berlii <- rbinom(p, 1, PP)
  index0 <- which(berlii == 0); index1 <- which(berlii == 1)
  true_pai1 <- length(index1)/p
  mu <- matrix(rep(0, 1*p), nrow=p)
  mu[index1] <- runif(length(index1), min=lb, max=ub)
  
  X <- matrix(rep(0, n*p), nrow=p)
  B <- matrix(runif(K*p, min=-1, max=1), nrow=p)
  
  sigma <- diag(p)
  sigma[index0,index0] <- diag(runif(length(index0), 1, 1))
  sigma[index1,index1] <- diag(c(2 + mu[index1]))
  mydata <- mvrnorm(n, meann, sigma)#n by p
  for (i in 1:n){
    f <- matrix(rnorm(K*1), nrow=K)
    X[,i] <- mu + B %*% f + mydata[i,]########
  }
  
  X <- t(X)
  RE2 <- pfa.test(X, tval= "pval",mat_est="sample", plot = "linear")#poet
  summary(RE2)
  cat(jj)
  hh[jj] <- RE2$K
  ttrue_pai1[jj] <- true_pai1
  ppai1_hat[jj] <- 1 - RE2$pi0
  
  index_FDP <- which(RE2$FDP$FDP<=gama)
  index_FDP <- index_FDP[length(index_FDP)]
  if (RR[jj] == 0){
    FFDP[jj] <- 0
  } else {
    FFDP[jj] <- RE2$FDP$FDP[index_FDP]
  }
  
  
  t <- RE2$FDP$t[index_FDP]
  if (length(RE2$FDP$rejects[index_FDP]) == 0){
    RR[jj] <- 0 
    St[jj] <- 0
    ttrue_FDP[jj] <- 0
  } else{
    RR[jj] <- RE2$FDP$rejects[index_FDP]
    St[jj] <- sum(RE2$FDP$pvaluess[index1]<=t)
    ttrue_FDP[jj] <- (RR[jj]-St[jj])/RR[jj]
  }
  Recive[jj] <- p - RR[jj]
  
  Ft[jj] <- sum(RE2$FDP$pvaluess[index1]>t)
  
  NNDR[jj] <- St[jj]/length(index1)
  
}
(true_FDP <- mean(ttrue_FDP[1:jj], na.rm = T))
(power2 <- mean(NNDR[1:jj], na.rm = T))

t=qt((1-0.05)/2 + .5, length(ttrue_FDP[1:jj])-1)
se_FDP = sd(ttrue_FDP[1:jj], na.rm = T) / sqrt(length(ttrue_FDP[1:jj]))
(CI_FDP <- t*se_FDP)

se_power = sd(NNDR[1:jj], na.rm = T) / sqrt(length(NNDR[1:jj]))
(CI_power <- t*se_power)
t2 <- Sys.time()
(t2-t1)