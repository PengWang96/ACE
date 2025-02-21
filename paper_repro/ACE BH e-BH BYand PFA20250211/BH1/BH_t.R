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

# 设定并行计算的线程数
numCores <- detectCores() - 1  # 获取可用的CPU核心数，减1以保留一个核心给操作系统
cl <- makeCluster(numCores)    # 创建线程
registerDoParallel(cl)         # 注册线程


# 初始化变量
t1 <- Sys.time()
# p <- 2000; n <- 600; K <- 3; lb <- 0.1; ub <- 0.4; v <- 4; sqrt_factor <- sqrt(v / (v - 2))
p <- 2000; n <- 100 ; K <- 3; lb <- 0.3; ub <- 0.6; v <- 4; sqrt_factor <- sqrt(v / (v - 2))
# p <- 2000; n <- 50 ; K <- 3; lb <- 0.3; ub <- 0.6; v <- 4; sqrt_factor <- sqrt(v / (v - 2))
# p <- 2000; n <- 10 ; K <- 3; lb <- 1.6; ub <- 1.9; v <- 4; sqrt_factor <- sqrt(v / (v - 2))
PP <- 0.2; gama <- 0.05

iteraton <- 100

# 并行计算
results <- foreach(jj = 1:iteraton, .combine = 'rbind', .packages = c('MASS', 'mvtnorm')) %dopar% {
  
  berlii <- rbinom(p, 1, PP)
  index0 <- which(berlii == 0)
  index1 <- which(berlii == 1)
  true_pai1 <- length(index1) / p
  mu <- matrix(rep(0, 1 * p), nrow = p)
  mu[index1] <- runif(length(index1), min = lb, max = ub)
  
  X <- matrix(rep(0, n * p), nrow = p)
  B <- matrix(runif(K*p, min=-1, max=1), nrow=p)
  # B <- matrix(runif(K*p, min=0, max=0), nrow=p)
  
  t_error <- t(rmvt(n, sigma = diag(p), df = 4))########
  # t_error <- t(rmvt(n, sigma = diag(p), df = 10))########
  for (i in 1:n){
    # f <- matrix(rnorm(K*1), nrow=K)
    
    f <- t(rmvt(1, sigma = diag(K), df = v)) / sqrt_factor
    
    # fx <- runif(K, min=0, max=1); fy <- runif(K, min=0, max=1)
    # f <- log(fx/fy)/sqrt(2)
    
    X[,i] <- mu + B %*% f + t_error[,i]# 
  }
  
  X <- t(X)
  
  # 计算 p 值
  pvalue <- sapply(1:p, function(x) {
    tt <- t.test(X[, x])
    return(tt$p.value)
  })
  
  sort_p <- sort(pvalue, index.return = TRUE)
  index_i <- (sort_p$x <= c(1:p) * gama / p)
  
  # 计算正确拒绝的数量
  correct_reject <- sapply(1:sum(index_i), function(aa) {
    if (sum(index1 == sort_p$ix[1:sum(index_i)][aa]) == 1) {
      return(1)
    } else {
      return(0)
    }
  })
  
  St_val <- sum(correct_reject)
  RR_val <- sum(index_i)
  
  NNDR_val <- St_val / length(index1)
  
  if (RR_val == 0) {
    FFDP_val <- 0
  } else {
    FFDP_val <- (RR_val - St_val) / RR_val
  }
  
  # 返回每个迭代的结果
  c(true_pai1, FFDP_val, St_val, RR_val, NNDR_val)
}

# 关闭并行处理
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
