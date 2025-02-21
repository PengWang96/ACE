rm(list = ls())
setwd("E:/RR/Save2022_10_14/BH1")
t1 <- Sys.time()
p <- 2000; n <- 600 ; K <- 3; lb <- 0.1; ub <- 0.4; v <- 4
# p <- 2000; n <- 50 ; K <- 3; lb <- 0.3; ub <- 0.6; v <- 4
PP <- 0.2
gama<-0.05
sigmaa <- read.csv("sigma_e.csv"); sigmaa <- sigmaa[1:p,1:p]
library(MASS);library(mvtnorm)
iteraton <- 100
ttrue_pai1 <- FFDP <- ttrue_FDP <- NNDR <- St <- RR <- rep(0,iteraton)
for (jj in 1:iteraton){
  berlii <- rbinom(p, 1, PP)
  index0 <- which(berlii == 0); index1 <- which(berlii == 1)
  true_pai1 <- length(index1)/p
  mu <- matrix(rep(0, 1*p), nrow=p)
  mu[index1] <- runif(length(index1), min=lb, max=ub)
  
  X <- matrix(rep(0, n*p), nrow=p)
  B <- matrix(runif(K*p, min=-1, max=1), nrow=p)
  # B <- matrix(runif(K*p, min=0, max=0), nrow=p)
  
  mydata <- mvrnorm(n, rep(0,p), sigmaa)#n by p
  for (i in 1:n){
    # f <- matrix(rnorm(K*1), nrow=K)
    
    f <- t(rmvt(1, sigma = diag(K), df = v))/sqrt(v/(v-2))
    
    # fx <- runif(K, min=0, max=1); fy <- runif(K, min=0, max=1)
    # f <- log(fx/fy)/sqrt(2)
    
    X[,i] <- mu + B %*% f + mydata[i,]# 
  }
  X <- t(X)
  
  pvalue <- sapply(x <- 1:p, 
                   function(x, X){
                     tt <- t.test(X[,x])
                     return(c(tt$p.value))
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