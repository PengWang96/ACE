rm(list = ls())
t1 <- Sys.time()
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

# 参数配置 ------------------------------------------------------------
param_combinations <- expand.grid(
  PP = c(0.2, 0.05),
  factor_indicator = c(1, 0),
  scenario = list(
    list(n = 600, lb = 0.1, ub = 0.4),  # 场景1
    list(n = 50, lb = 0.6, ub = 0.9)    # 场景2
  ),
  stringsAsFactors = FALSE
)

p <- 2000
K <- 3
v <- 4
gama <- 0.05
iteraton <- 100
sqrt_factor <- sqrt(v/(v-2))

# 并行初始化 ----------------------------------------------------------
numCores <- detectCores() - 5  # 保留2个核心给系统
cl <- makeCluster(numCores)
registerDoParallel(cl)
clusterExport(cl, c("K", "v", "gama", "sqrt_factor", "p"))  # 导出全局变量

# 核心计算函数 --------------------------------------------------------
simulation_task <- function(params, iter) {
  tryCatch({
    # 解包参数
    PP <- params$PP
    factor_indicator <- params$factor_indicator
    n <- params$scenario$n
    lb <- params$scenario$lb
    ub <- params$scenario$ub
    
    # 数据生成
    berlii <- rbinom(p, 1, PP)
    index0 <- which(berlii == 0)
    index1 <- which(berlii == 1)
    mu <- matrix(0, p, 1)
    mu[index1] <- runif(length(index1), lb, ub)
    
    # 矩阵化加速计算
    B <- factor_indicator * matrix(runif(K*p, -1, 1), p, K)
    f_matrix <- t(rmvt(n, sigma = diag(K), df = v)) / sqrt_factor
    
    sigma <- diag(p)
    sigma[index0,index0] <- diag(runif(length(index0), 1, 1))
    if (length(index1) == 1) {
      sigma[index1,index1] <- 2 + mu[index1]
    } else {
      sigma[index1,index1] <- diag(c(2 + mu[index1]))
    }
    
    # error_matrix <- t(rmvt(n, sigma = sigma, df = 4))
    error_matrix <- t(rmvt(n, sigma = sigma, df = 10)) #######################
    
    # 向量化计算替代循环
    X <- mu %*% rep(1, n) + B %*% f_matrix + error_matrix
    
    # 并行计算p值
    pvalue <- apply(X, 1, function(x) t.test(x)$p.value)
    p_adj <- p.adjust(pvalue, "BY")
    
    # 结果计算
    reject <- which(p_adj < gama)
    correct <- sum(reject %in% index1)
    RR <- length(reject)
    
    c(
      FFDP = ifelse(RR==0, 0, (RR - correct)/RR),
      Power = correct/max(length(index1), 1)
    )
  }, error = function(e) return(c(FFDP=NA, Power=NA)))
}

# 主执行流程 ----------------------------------------------------------
results <- foreach(i = 1:nrow(param_combinations), .combine = rbind) %:% 
  foreach(j = 1:iteraton, .combine = rbind, .packages = c("MASS", "mvtnorm")) %dopar% {
    params <- list(
      PP = param_combinations$PP[i],
      factor_indicator = param_combinations$factor_indicator[i],
      scenario = param_combinations$scenario[[i]]
    )
    res <- simulation_task(params, j)
    data.frame(
      PP = params$PP,
      factor_indicator = params$factor_indicator,
      n = params$scenario$n,
      lb = params$scenario$lb,
      ub = params$scenario$ub,
      iter = j,
      FFDP = res["FFDP"],
      Power = res["Power"]
    )
  }

stopCluster(cl)

# 结果汇总 ------------------------------------------------------------
library(dplyr)
final_results <- results %>%
  group_by(PP, factor_indicator, n, lb, ub) %>%
  summarise(
    FFDP_mean = mean(FFDP, na.rm = TRUE),
    Power_mean = mean(Power, na.rm = TRUE),
    FFDP_ci = 1.96 * sd(FFDP, na.rm = TRUE)/sqrt(n()),
    Power_ci = 1.96 * sd(Power, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )


final_results <- final_results %>%
  arrange(
    desc(PP),  # 先排列0.2再0.05
    desc(factor_indicator),  # 先排列1再0
    desc(n),  # 先排列600再50
    lb, ub
  ) %>% 
  # 添加显式的场景编号
  mutate(
    scenario_id = case_when(
      n == 600 & lb == 0.1 ~ 1,
      n == 50 & lb == 0.6 ~ 2,
      TRUE ~ NA_real_
    )
  ) %>% 
  arrange(scenario_id) %>% 
  select(-scenario_id)

print(final_results)
# write.csv(final_results, "./output/heter_BY_t4.csv", row.names = FALSE)
write.csv(final_results, "./output/heter_BY_t10.csv", row.names = FALSE)
t2 <- Sys.time()
(t2-t1)