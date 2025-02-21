rm(list = ls())
setwd("E:/RR/Save2022_10_14/BH1")
t1 <- Sys.time()
# p <- 2000; n1 <- 600 ; n2 <- 800; K1 <- 1; K2 <- 2; lb <- 0.1; ub <- 0.4; v <- 4
p <- 2000; n1 <- 50 ; n2 <- 100; K1 <- 1; K2 <- 2; lb <- 0.8; ub <- 1.1; v <- 4
PP <- 0.2
gama<-0.05
sigmaa2 <- read.csv("sigma_e.csv"); sigmaa2 <- sigmaa2[1:p,1:p]
library("MASS");library(mvtnorm)
iteraton <- 100
ttrue_pai1 <- FFDP <- ttrue_FDP <- NNDR <- St <- RR <- rep(0,iteraton)
for (jj in 1:iteraton){
  berlii <- rbinom(p, 1, PP)
  index0 <- which(berlii == 0); index1 <- which(berlii == 1)
  true_pai1 <- length(index1)/p
  
  mu <- matrix(rep(0, 1*p), nrow=p)
  mu[index0] <- runif(length(index0), 1, 4)##########
  mu[index1] <- runif(length(index1), 1, 4)##########
  mu2 <- mu;
  mu[index1] <- runif(length(index1), min=lb, max=ub) + mu[index1]
  
  X <- matrix(rep(0, n1*p), nrow=p)
  B <- matrix(runif(K1*p, min=-1, max=1), nrow=p)
  Y <- matrix(rep(0, n2*p), nrow=p)
  B2 <- matrix(runif(K2*p, min=-1, max=1), nrow=p)
  
  sigma1 <- diag(c(2 + mu - mu2)); sigmaa2 <- sigmaa2 - diag(rep(1,p)) + diag(c(2.5 + mu - mu2))
  diag(sigma1[index0,index0]) <- (runif(length(index0), min=1, max=1))
  diag(sigmaa2[index0,index0]) <- (runif(length(index0), min=1, max=1))
  # t_error <- t(rmvt(n1, sigma = sigma1, df = 4))
  t_error <- t(rmvt(n1, sigma = sigma1, df = 10))################
  mydata2 <- mvrnorm(n2, rep(0,p), sigmaa2)#n by p
  
  for (i in 1:n1){
    # f <- matrix(rnorm(K1*1), nrow=K1)
    
    # f <- t(rmvt(1, sigma = diag(K1), df = v))/sqrt(v/(v-2))
    
    fx <- runif(K1, min=0, max=1); fy <- runif(K1, min=0, max=1)
    f <- log(fx/fy)/sqrt(2)
    
    X[,i] <- mu + B %*% f + t_error[,i]########
  }
  for (i in 1:n2){
    # f2 <- matrix(rnorm(K2*1), nrow=K2)
    
    # f2 <- t(rmvt(1, sigma = diag(K2), df = v))/sqrt(v/(v-2))
    
    fx <- runif(K2, min=0, max=1); fy <- runif(K2, min=0, max=1)
    f2 <- log(fx/fy)/sqrt(2)
    Y[,i] <- mu2 + B2 %*% f2 + mydata2[i,]######### 
  }
  
  X <- t(X); Y <- t(Y)
  
  pvalue <- sapply(x <- 1:p, 
                   function(x, X){
                     tt <- t.test(X[,x],Y[,x])
                     return(c(tt$p.value))#tt$statistic,
                   }
                   ,X)
  
  sort_p <- sort(pvalue,index.return = T)
  index_i <- (sort_p$x <= c(1:p)*gama/p)
  
  
  correct_reject <- sapply(aa <- 1:length(sort_p$ix[1:sum(index_i)]), #
                           function(aa, bb){
                             if (sum(bb == sort_p$ix[1:sum(index_i)][aa]) == 1){
                               c_r <- 1
                             } else{
                               c_r <- 0
                             }
                             return(c_r)
                           }
                           ,bb = index1)
  
  St[jj] <- sum(correct_reject)
  RR[jj] <- sum(index_i)
  
  ttrue_pai1[jj] <- true_pai1
  NNDR[jj] <- St[jj]/length(index1)
  if (RR[jj] == 0){
    FFDP[jj] <- 0
  } else {
    FFDP[jj] <- (RR[jj]-St[jj])/RR[jj]
  }
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