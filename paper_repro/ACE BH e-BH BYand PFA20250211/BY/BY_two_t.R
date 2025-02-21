rm(list = ls())
t1 <- Sys.time()
library(foreach)
library(doParallel)
library(MASS)
library(mvtnorm)
library(rstudioapi)

if (interactive()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
}
setwd(dirname(current_file))

# 参数配置 ------------------------------------------------------------
param_combinations <- expand.grid(
  PP = c(0.2, 0.05),
  factor_indicator = c(1, 0),
  scenario = list(
    list(n1 = 600, n2 = 800, lb = 0.1, ub = 0.4),  # 场景1
    list(n1 = 50, n2 = 100, lb = 0.6, ub = 0.9)    # 场景2
  ),
  stringsAsFactors = FALSE
)

p <- 2000
K1 <- 1
K2 <- 2
v <- 4
gama <- 0.05
iteraton <- 100
sqrt_factor <- sqrt(v/(v-2))

# 并行初始化 ----------------------------------------------------------
numCores <- detectCores() - 2
cl <- makeCluster(numCores)
registerDoParallel(cl)
clusterExport(cl, c("p", "K1", "K2", "v", "gama", "sqrt_factor"))

# 核心计算函数 --------------------------------------------------------
simulation_task <- function(params) {
  tryCatch({
    # 解包参数
    PP <- params$PP
    factor_indicator <- params$factor_indicator
    n1 <- params$scenario$n1
    n2 <- params$scenario$n2
    lb <- params$scenario$lb
    ub <- params$scenario$ub
    
    # 数据生成
    berlii <- rbinom(p, 1, PP)
    index1 <- which(berlii == 1)
    
    # 生成mu
    mu <- matrix(runif(p, 1, 4), ncol = 1)
    mu2 <- mu
    mu[index1] <- mu[index1] + runif(length(index1), lb, ub)
    
    # 生成因子载荷矩阵
    B <- factor_indicator * matrix(runif(K1*p, -1, 1), nrow = p)
    B2 <- factor_indicator * matrix(runif(K2*p, -1, 1), nrow = p)
    
    # 生成数据矩阵
    generate_data <- function(n, K, B, mu, df) {
      f_matrix <- t(rmvt(n, sigma = diag(K), df = v)) / sqrt_factor
      error_matrix <- t(rmvt(n, sigma = diag(p), df = df))
      mu %*% rep(1, n) + B %*% f_matrix + error_matrix
    }
    
    X <- generate_data(n1, K1, B, mu, df = 4) ################################
    Y <- generate_data(n2, K2, B2, mu2, df = 4) ##############################
    
    # X <- generate_data(n1, K1, B, mu, df = 4) ################################
    # Y <- generate_data(n2, K2, B2, mu2, df = 10) #############################

    # X <- generate_data(n1, K1, B, mu, df = 10) ###############################
    # Y <- generate_data(n2, K2, B2, mu2, df = 10) #############################
    
    # 计算p值和e值
    pvalue <- sapply(1:p, function(x) t.test(X[x, ], Y[x, ])$p.value)
    p_adj <- p.adjust(pvalue, "BY")
    
    # 结果计算
    reject <- which(p_adj < gama)
    correct <- sum(reject %in% index1)
    RR <- length(reject)
    
    c(
      FFDP = ifelse(RR == 0, 0, (RR - correct)/RR),
      Power = correct/max(length(index1), 1)
    )
  }, error = function(e) return(c(FFDP = NA, Power = NA)))
}

# 主执行流程 ----------------------------------------------------------
results <- foreach(i = 1:nrow(param_combinations), .combine = rbind) %:% 
  foreach(j = 1:iteraton, .combine = rbind, .packages = c("MASS", "mvtnorm")) %dopar% {
    params <- list(
      PP = param_combinations$PP[i],
      factor_indicator = param_combinations$factor_indicator[i],
      scenario = param_combinations$scenario[[i]]
    )
    res <- simulation_task(params)
    data.frame(
      PP = params$PP,
      factor_indicator = params$factor_indicator,
      n1 = params$scenario$n1,
      n2 = params$scenario$n2,
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
  group_by(PP, factor_indicator, n1, n2, lb, ub) %>%
  summarise(
    FFDP_mean = mean(FFDP, na.rm = TRUE),
    Power_mean = mean(Power, na.rm = TRUE),
    FFDP_ci = 1.96 * sd(FFDP, na.rm = TRUE)/sqrt(n()),
    Power_ci = 1.96 * sd(Power, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  ) %>%
  arrange(
    desc(PP),
    desc(factor_indicator),
    desc(n1),
    lb, ub
  ) %>% 
  mutate(
    scenario_id = case_when(
      n1 == 600 ~ 1,
      n1 == 50 ~ 2,
      TRUE ~ NA_real_
    )
  ) %>% 
  arrange(scenario_id) %>% 
  select(-scenario_id)

# 输出结果
print(final_results)
write.csv(final_results, "./output/BY_two_t4_4.csv", row.names = FALSE)
# write.csv(final_results, "./output/BY_two_t4_10.csv", row.names = FALSE)
# write.csv(final_results, "./output/BY_two_t10_10.csv", row.names = FALSE)
t2 <- Sys.time()
(t2-t1)