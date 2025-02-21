rm(list = ls())
Z <- read.csv("./data/expr_treatment.csv", row.names = 1, stringsAsFactors = FALSE)
X <- read.csv("./data/expr_control.csv", row.names = 1, stringsAsFactors = FALSE)
gene_symbols <- read.csv("./data/gene_symbols.csv", stringsAsFactors = FALSE)
gamma <- seq(0.05, 0.15, 0.01)
gama <- gamma[1]

# hist(rowMeans(Z) - rowMeans(X))
# dat <- data.frame(sample_mean = rowMeans(Z) - rowMeans(X))
# library(ggplot2)
# ggplot(dat, aes(x=sample_mean)) + 
#   geom_histogram(bins = 20)

################# ACE #################
# 4745
t1 <- Sys.time()
library(ACE)
result_ace <- ACE(Z, X, gama = gama)
differentially_expressed_index <- which(abs(result_ace$Adjusted_statistics) > result_ace$Threshold)
num_differentially_expressed <- result_ace$Rejection
ACE_result <- list(index = differentially_expressed_index,
                   count = num_differentially_expressed)
cat("差异基因的索引：\n", differentially_expressed_index, "\n")
cat("差异基因的个数：", num_differentially_expressed, "\n")
t2 <- Sys.time()
(t2 - t1)

median((result_ace$Adjusted_statistics)^2)/qchisq(0.5, df=1)
# gene_symbols$SYMBOL[ACE_result$index]

# library(ggplot2)
# result_ace <- ACE(Z, X, gama = gama)
# statistics <- result_ace$Adjusted_statistics
# 
# qq_data <- data.frame(statistics = statistics)
# 
# p_gg <- ggplot(qq_data, aes(sample = statistics)) +
#   stat_qq() +
#   stat_qq_line(col = "red") +
#   theme_bw() +
#   theme(text = element_text(size = 10),  # 设置字体大小
#         plot.title = element_text(size = 15, hjust = 0.5),  # 加大标题字体
#         axis.title.x = element_blank(),  # 去掉 x 轴标题
#         axis.title.y = element_blank(),  # 去掉 y 轴标题
#         axis.text.x = element_text(size = 10),  # 设置坐标轴刻度字体大小
#         axis.text.y = element_text(size = 10))   # 设置坐标轴刻度字体大小
# p_gg



################# PFA #################
t1 <- Sys.time()
source("./scripts/use.R")
source("./scripts/core.R")
# RE2 <- pfa.test(t(Z), t(X), tval= "pval", mat_est="sample", plot = "linear") #
# summary(RE2) ##########################################

load("./output/PFA_all_results.rda") ##########################################

ppai1_hat <- 1 - RE2$pi0
index_FDP <- which(RE2$FDP$FDP<=gama)
differentially_expressed_index <- index_FDP
index_FDP <- index_FDP[length(index_FDP)]
FFDP <- RE2$FDP$FDP[index_FDP]
t <- RE2$FDP$t[index_FDP]
num_differentially_expressed <- sum(RE2$FDP$pvaluess<=t)
(FFDP)
(num_differentially_expressed)
(ppai1_hat)

# save(RE2, file = "./output/PFA_all_results.rda") #############################
PFA_result <- list(index = differentially_expressed_index,
                   count = num_differentially_expressed)
cat("差异基因的索引：\n", differentially_expressed_index, "\n")
cat("差异基因的个数：", num_differentially_expressed, "\n")
t2 <- Sys.time()
(t2 - t1)



################# BH #################
t1 <- Sys.time()
common_genes <- intersect(rownames(Z), rownames(X))
Z <- Z[common_genes, ]
X <- X[common_genes, ]
p_values <- numeric(length(common_genes))

# 使用 sapply 进行向量化操作
p_values <- sapply(common_genes, function(gene) {
  t.test(Z[gene, ], X[gene, ])$p.value
})

p_values_adj <- p.adjust(p_values, method = "BH")
# hist(p_values)
# hist(p_values_adj)
differentially_expressed_index <- which(p_values_adj < gama)
num_differentially_expressed <- length(differentially_expressed_index)

BH_result <- list(index = differentially_expressed_index,
                  count = num_differentially_expressed)
cat("差异基因的索引：\n", differentially_expressed_index, "\n")
cat("差异基因的个数：", num_differentially_expressed, "\n")
t2 <- Sys.time()
(t2 - t1)




################# e-BH #################
t1 <- Sys.time()
source("./scripts/e-BH.R")
# p-value
pvalue <- sapply(common_genes, function(gene) {
  t.test(Z[gene, ], X[gene, ])$p.value
})
# e-value
e_values <- (pvalue^(-1/2))/2
boost <- (2/gama)^(1/2)
e_values <- e_values * boost

differentially_expressed_index <- e_bh(e_values, gama)
num_differentially_expressed <- length(differentially_expressed_index)

e_BH_result <- list(index = differentially_expressed_index,
                    count = num_differentially_expressed)
cat("差异基因的索引：\n", differentially_expressed_index, "\n")
cat("差异基因的个数：", num_differentially_expressed, "\n")
t2 <- Sys.time()
(t2 - t1)





################# BY #################
t1 <- Sys.time()
common_genes <- intersect(rownames(Z), rownames(X))
Z <- Z[common_genes, ]
X <- X[common_genes, ]
p_values <- numeric(length(common_genes))

# 使用 sapply 进行向量化操作
p_values <- sapply(common_genes, function(gene) {
  t.test(Z[gene, ], X[gene, ])$p.value
})

p_values_adj <- p.adjust(p_values, method = "BY")
differentially_expressed_index <- which(p_values_adj < gama)
num_differentially_expressed <- length(differentially_expressed_index)

BY_result <- list(index = differentially_expressed_index,
                  count = num_differentially_expressed)
cat("差异基因的索引：\n", differentially_expressed_index, "\n")
cat("差异基因的个数：", num_differentially_expressed, "\n")
t2 <- Sys.time()
(t2 - t1)

save(ACE_result, PFA_result, BH_result, e_BH_result, BY_result, 
     file = paste0("./output/differential_results_gama_", gama, ".rda"))