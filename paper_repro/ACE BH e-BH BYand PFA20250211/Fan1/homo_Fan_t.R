rm(list = ls())
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
source("use.R"); source("core.R")

# 设定并行计算的线程数
numCores <- detectCores() - 1  # 获取可用的CPU核心数，减1以保留一个核心用于操作系统
cl <- makeCluster(numCores)    # 创建线程
registerDoParallel(cl)         # 注册线程

# 初始化变量
t1 <- Sys.time()
# p <- 2000; n <- 600; K <- 3; lb <- 0.1; ub <- 0.4; v <- 4; sqrt_factor <- sqrt(v / (v - 2))
p <- 2000; n <- 100 ; K <- 3; lb <- 0.3; ub <- 0.6; v <- 4; sqrt_factor <- sqrt(v / (v - 2))
# p <- 2000; n <- 50 ; K <- 3; lb <- 0.3; ub <- 0.6; v <- 4; sqrt_factor <- sqrt(v / (v - 2))
# p <- 2000; n <- 10 ; K <- 3; lb <- 1.6; ub <- 1.9; v <- 4; sqrt_factor <- sqrt(v / (v - 2))
PP <- 0.2; gama <- 0.05

iteraton <- 20

# 并行运算
results <- foreach(jj = 1:iteraton, .combine = 'rbind', .packages = c('MASS', 'mvtnorm', 'quantreg')) %dopar% {
  
  berlii <- rbinom(p, 1, PP)
  index0 <- which(berlii == 0)
  index1 <- which(berlii == 1)
  true_pai1 <- length(index1)/p
  mu <- matrix(rep(0, 1 * p), nrow = p)
  mu[index1] <- runif(length(index1), min = lb, max = ub)
  
  X <- matrix(rep(0, n * p), nrow = p)
  B <- matrix(runif(K * p, min = -1, max = 1), nrow = p)
  # B <- matrix(runif(K*p, min=0, max=0), nrow=p)
  
  t_error <- t(rmvt(n, sigma = diag(p), df = 4))########
  # t_error <- t(rmvt(n, sigma = diag(p), df = 10))########
  for (i in 1:n){
    f <- t(rmvt(1, sigma = diag(K), df = v)) / sqrt_factor
    X[,i] <- mu + B %*% f + t_error[, i] # 
  }
  
  X <- t(X)
  
  
  # 调用外部的pfa.test函数
  RE2 <- pfa.test(X, tval = "pval", mat_est = "sample", plot = "linear") # , reg = "L2"
  
  # 初始化结果
  RR_val <- St_val <- ttrue_FDP_val <- FFDP_val <- Recive_val <- NNDR_val <- 0
  ppai1_hat_val <- 1 - RE2$pi0
  
  index_FDP <- which(RE2$FDP$FDP <= gama)
  if (length(index_FDP) > 0) {
    index_FDP <- index_FDP[length(index_FDP)]
    FFDP_val <- RE2$FDP$FDP[index_FDP]
    t <- RE2$FDP$t[index_FDP]
    
    if (length(RE2$FDP$rejects[index_FDP]) > 0) {
      RR_val <- RE2$FDP$rejects[index_FDP]
      St_val <- sum(RE2$FDP$pvaluess[index1] <= t)
      ttrue_FDP_val <- (RR_val - St_val) / RR_val
      NNDR_val <- St_val / length(index1)
    }
  }
  
  Recive_val <- p - RR_val
  
  # 返回每个迭代的结果
  c(true_pai1, ppai1_hat_val, FFDP_val, ttrue_FDP_val, RR_val, St_val, Recive_val, NNDR_val)
}

stopCluster(cl)

# 提取结果
ttrue_pai1 <- results[, 1]
ppai1_hat <- results[, 2]
FFDP <- results[, 3]
ttrue_FDP <- results[, 4]
RR <- results[, 5]
St <- results[, 6]
Recive <- results[, 7]
NNDR <- results[, 8]

# 计算FDP和Power的均值
true_FDP <- mean(ttrue_FDP, na.rm = TRUE)
power2 <- mean(NNDR, na.rm = TRUE)

# 计算置信区间
t <- qt((1 - 0.05) / 2 + .5, length(ttrue_FDP) - 1)
se_FDP <- sd(ttrue_FDP, na.rm = TRUE) / sqrt(length(ttrue_FDP))
CI_FDP <- t * se_FDP

se_power <- sd(NNDR, na.rm = TRUE) / sqrt(length(NNDR))
CI_power <- t * se_power

t2 <- Sys.time()
print(t2 - t1)

# 输出结果
list(true_FDP = true_FDP, power2 = power2, CI_FDP = CI_FDP, CI_power = CI_power)
