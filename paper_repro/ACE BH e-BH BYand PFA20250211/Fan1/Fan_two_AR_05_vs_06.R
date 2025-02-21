rm(list = ls())
setwd("E:/RR/Save2022_10_14/Fan1")
t1 <- Sys.time()
# p <- 2000; n1 <- 600 ; n2 <- 800; K1 <- 1; K2 <- 2; lb <- 0.1; ub <- 0.4; v<-4
p <- 2000; n1 <- 50 ; n2 <- 100; K1 <- 1; K2 <- 2; lb <- 0.6; ub <- 0.9; v<-4
PP <- 0.2
gama<-0.05
sigmaa <- read.csv("sigma_e.csv"); sigmaa <- sigmaa[1:p,1:p]
sigmaa2 <- read.csv("sigma_e_06.csv"); sigmaa2 <- sigmaa2[1:p,1:p]
source("use.R");source("core.R");library(MASS)
iteraton <- 100
ttrue_pai1 <- ppai1_hat <- FFDP <- ttrue_FDP <- NNDR <- FFNR <- Ft <- St <- 
  RR <- Recive <- ppower2 <- rep(0,iteraton)
for (jj in 1:iteraton){
  berlii <- rbinom(p, 1, PP)
  index0 <- which(berlii == 0); index1 <- which(berlii == 1)
  true_pai1 <- length(index1)/p
  
  mu <- matrix(rep(0, 1*p), nrow=p)
  mu[index0] <- runif(length(index0), 1, 4)##########
  mu[index1] <- runif(length(index1), 1, 4)##########
  mu2 <- mu;
  mu[index1] <- runif(length(index1), min=lb, max=ub) + mu[index1];
  
  X <- matrix(rep(0, n1*p), nrow=p)
  B <- matrix(runif(K1*p, min=-1, max=1), nrow=p)
  Y <- matrix(rep(0, n2*p), nrow=p)
  B2 <- matrix(runif(K2*p, min=-1, max=1), nrow=p)
  
  mydata1 <- mvrnorm(n1, rep(0,p), sigmaa)#n by p
  mydata2 <- mvrnorm(n2, rep(0,p), sigmaa2)#n by p
  for (i in 1:n1){
    # f <- matrix(rnorm(K1*1), nrow=K1)
    
    # f <- t(rmvt(1, sigma = diag(K1), df = v))/sqrt(v/(v-2))
    
    fx <- runif(K1, min=0, max=1); fy <- runif(K1, min=0, max=1)
    f <- log(fx/fy)/sqrt(2)
    X[,i] <- mu + B %*% f + mydata1[i,]########
  }
  for (i in 1:n2){
    # f2 <- matrix(rnorm(K2*1), nrow=K2)
    
    # f2 <- t(rmvt(1, sigma = diag(K2), df = v))/sqrt(v/(v-2))
    
    fx <- runif(K2, min=0, max=1); fy <- runif(K2, min=0, max=1)
    f2 <- log(fx/fy)/sqrt(2)
    Y[,i] <- mu2 + B2 %*% f2 + mydata2[i,]#########
  }
  
  X <- t(X); Y <- t(Y)
  RE2 <- pfa.test(X,Y, tval= "pval" ,mat_est="sample", plot = "linear")#poet
  summary(RE2)
  cat(jj)
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