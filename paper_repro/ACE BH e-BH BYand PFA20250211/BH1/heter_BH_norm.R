rm(list = ls())
t1 <- Sys.time()
setwd("D:/R/factor model/Save2023_4_25/BH1")
p <- 2000; n <- 600 ; K <- 3; lb <- 0.1; ub <- 0.4; ita <- 1.5; v <- 4
# p <- 200; n <- 50 ; K <- 3; lb <- 0.3; ub <- 0.6; ita <- 1.5; v <- 4
PP <- 0.2 #true pi_1
gama <- 0.05 #control level
iteraton <- 100 # number of repeat
library(MASS)
meann <- rep(0,p)
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
  sigma[index1,index1] <- diag(c(2 + mu[index1]))
  mydata <- mvrnorm(n, meann, sigma)#n by p
  for (i in 1:n){
    f <- matrix(rnorm(K*1), nrow=K)
    
    # f <- t(rmvt(1, sigma = diag(K), df = v))/sqrt(v/(v-2))
    
    # fx <- runif(K, min=0, max=1); fy <- runif(K, min=0, max=1)
    # f <- log(fx/fy)/sqrt(2)
    
    X[,i] <- mu + B %*% f + mydata[i,]########
  }
  
  
  X <- t(X)
  
  pvalue <- sapply(x <- 1:p, 
                   function(x, X){
                     tt <- t.test(X[,x])
                     return(c(tt$p.value))
                   }
                   ,X)
  
  sort_p <- sort(pvalue,index.return = T)
  index_i <- (sort_p$x <= c(1:p)*gama/p)  # BH meathod
  
  
  correct_reject <- sapply(aa <- 1:length(sort_p$ix[1:sum(index_i)]), #
                           function(aa, bb){
                             if (sum(bb == sort_p$ix[1:sum(index_i)][aa]) == 1){
                               c_r <- 1
                             } else{
                               c_r <- 0
                             }
                             return(c_r)
                           }
                           ,bb = index1)  # calculate number of correct rejections
  
  St[jj] <- sum(correct_reject)
  RR[jj] <- sum(index_i)
  
  ttrue_pai1[jj] <- true_pai1
  NNDR[jj] <- St[jj]/length(index1)
  if (RR[jj] == 0){
    FFDP[jj] <- 0
  } else {
    FFDP[jj] <- (RR[jj]-St[jj])/RR[jj]
  }
}
(true_pai1 <- mean(ttrue_pai1))
(FDP <- mean(FFDP,na.rm = T))
(power <- mean(NNDR,na.rm = T))
(S <- mean(St))
(R <- mean(RR))
FDP <= (1 - true_pai1)*gama
t2 <- Sys.time()
(t2-t1)