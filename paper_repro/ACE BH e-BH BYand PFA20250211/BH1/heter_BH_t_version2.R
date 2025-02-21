rm(list = ls())
setwd("D:/R/factor model/Save2023_4_25/BH1")
t1 <- Sys.time()
p <- 2000; n <- 600 ; K <- 3; lb <- 0.1; ub <- 0.4; v <- 4
# p <- 20; n <- 10 ; K <- 3; lb <- 1; ub <- 3; v <- 4
# p <- 2000; n <- 50 ; K <- 3; lb <- 0.6; ub <- 0.9; v <- 4
PP <- 0.2 #true pi_1
gama <- 0.05 #control level
iteraton <- 10# number of repeat
library(MASS);library(mvtnorm)
ttrue_pai1 <- FFDP <- ttrue_FDP <- NNDR <- St <- RR <- rep(0,iteraton)
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
  if (length(index1) == 1) {
    sigma[index1,index1] <- 2 + mu[index1]
  } else {
    sigma[index1,index1] <- diag(c(2 + mu[index1]))
  }
  
  t_error <- t(rmvt(n, sigma = sigma, df = 4))########
  # t_error <- t(rmvt(n, sigma = sigma, df = 10))########
  for (i in 1:n){
    # f <- matrix(rnorm(K*1), nrow=K)
    
    f <- t(rmvt(1, sigma = diag(K), df = v))/sqrt(v/(v-2))
    
    # fx <- runif(K, min=0, max=1); fy <- runif(K, min=0, max=1)
    # f <- log(fx/fy)/sqrt(2)
    
    X[,i] <- mu + B %*% f + t_error[,i]# 
  }
  
  X <- t(X)
  
  Z <- sqrt(n)*colMeans(X)/sqrt(diag(cov(X)))
  pvalue <- 2 * pnorm(abs(Z), lower.tail = F)
  
  sort_p <- sort(pvalue,index.return = T)
  index_i <- (sort_p$x <= c(1:p)*gama/p)
  RR[jj] <- sum(index_i)
  
  Index_reject <- sort_p$ix[1:RR[jj]]
  false_reject <- intersect(Index_reject, index0)
  
  
  St[jj] <- RR[jj] - length(false_reject)
  NNDR[jj] <- St[jj]/length(index1)
  if (RR[jj] == 0){
    FFDP[jj] <- 0
  } else {
    FFDP[jj] <- (RR[jj]-St[jj])/RR[jj]
  }
  
  
  # RR[jj] <- sum(index_i)
  # ttrue_pai1[jj] <- true_pai1
  # if (RR[jj] == 0){
  #   St[jj] <- 0; FFDP[jj] <- 0; NNDR[jj] <- 0
  # } else {
  #   St[jj] <- RR[jj] - length(false_reject)
  #   FFDP[jj] <- (RR[jj]-St[jj])/RR[jj]
  #   NNDR[jj] <- St[jj]/length(index1)
  # }
  
  
  cat(jj)
}
(true_FDP <- mean(FFDP[1:jj], na.rm = T))
(power2 <- mean(NNDR[1:jj], na.rm = T))

t=qt((1-0.05)/2 + .5, length(FFDP[1:jj])-1)
se_FDP = sd(FFDP[1:jj], na.rm = T) / sqrt(length(FFDP[1:jj]))
(CI_FDP <- t*se_FDP)

se_power = sd(NNDR[1:jj], na.rm = T) / sqrt(length(NNDR[1:jj]))
(CI_power <- t*se_power)
t2 <- Sys.time()
(t2-t1)