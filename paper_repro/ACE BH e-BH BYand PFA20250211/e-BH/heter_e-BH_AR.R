rm(list = ls())
library(beepr)
library(foreach)
library(doParallel)
library(MASS)
library(mvtnorm)
library(rstudioapi)
library(rprojroot)
if (interactive()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
}
setwd(dirname(current_file))
print(getwd())
source("e-BH.R")
# 设定并行计算的线程数
numCores <- detectCores() - 1  # 获取可用的CPU核心数，减1以保留一个核心给操作系统
cl <- makeCluster(numCores)    # 创建线程
registerDoParallel(cl)         # 注册线程

# 初始化变量
t1 <- Sys.time()
# p <- 2000; n <- 600 ; K <- 3; lb <- 0.1; ub <- 0.4; v <- 4
p <- 2000; n <- 50 ; K <- 3; lb <- 0.6; ub <- 0.9; v <- 4
PP <- 0.05              #######
factor_indicator <- 0  #######
gama <- 0.05
iteraton <- 100
sqrt_factor <- sqrt(v / (v - 2))
sigmaa <- read.csv("sigma_e.csv"); sigmaa <- sigmaa[1:p,1:p]
# 并行计算
results <- foreach(jj = 1:iteraton, .combine = 'rbind', .packages = c('MASS', 'mvtnorm')) %dopar% {
  
  berlii <- rbinom(p, 1, PP)
  index0 <- which(berlii == 0)
  index1 <- which(berlii == 1)
  true_pai1 <- length(index1) / p
  mu <- matrix(rep(0, 1 * p), nrow = p)
  mu[index1] <- runif(length(index1), min = lb, max = ub)
  
  X <- matrix(rep(0, n * p), nrow = p)
  B <- factor_indicator * matrix(runif(K*p, min=-1, max=1), nrow=p)
  
  sigma <- diag(p)
  sigma[index0,index0] <- diag(runif(length(index0), min=1, max=1))
  sigma[index1,index1] <- diag(c(2 + mu[index1]))
  sigma <- sigmaa + sigma - diag(p)
  
  mydata <- mvrnorm(n, rep(0,p), sigma)#n by p
  for (i in 1:n){
    f <- t(rmvt(1, sigma = diag(K), df = v))/sqrt_factor
    X[,i] <- mu + B %*% f + mydata[i,]# 
  }
  
  X <- t(X)
  
  # p-value
  pvalue <- sapply(1:p, function(x) {
    tt <- t.test(X[, x])
    return(tt$p.value)
  })
  # e-value
  e_values <- (pvalue^(-1/2))/2
  boost <- (2/gama)^(1/2) # 8.94 #
  e_values <- e_values * boost
  
  reject_index <- e_bh(e_values, gama)
  # 计算正确拒绝的数量
  correct_reject <- sum(reject_index %in% index1)
  
  St_val <- correct_reject
  RR_val <- length(reject_index)
  
  NNDR_val <- St_val / length(index1)
  
  if (RR_val == 0) {
    FFDP_val <- 0
  } else {
    FFDP_val <- (RR_val - St_val) / RR_val
  }
  
  # 返回每个迭代的结果
  c(true_pai1, FFDP_val, St_val, RR_val, NNDR_val)
}
stopCluster(cl)

# 提取结果
ttrue_pai1 <- results[, 1]
FFDP <- results[, 2]
St <- results[, 3]
RR <- results[, 4]
NNDR <- results[, 5]

# 计算FDP和Power的均值
true_FDP <- mean(FFDP, na.rm = TRUE)
power2 <- mean(NNDR, na.rm = TRUE)

# 计算置信区间
t <- qt((1 - 0.05) / 2 + .5, length(FFDP) - 1)
se_FDP <- sd(FFDP, na.rm = TRUE) / sqrt(length(FFDP))
CI_FDP <- t * se_FDP

se_power <- sd(NNDR, na.rm = TRUE) / sqrt(length(NNDR))
CI_power <- t * se_power

t2 <- Sys.time()
time_taken <- t2 - t1
# 输出结果
list(
  true_FDP = true_FDP,
  power2 = power2,
  CI_FDP = CI_FDP,
  CI_power = CI_power,
  time_taken = time_taken
)
beep(sound = "mario")