rm(list = ls())
Z <- read.csv("./data/Vbench.csv")
X <- read.csv("./data/Ebench.csv")
gamma <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05)
gama <- gamma[1]

hist(rowMeans(Z))
hist(rowMeans(X))
hist(rowMeans(X) - rowMeans(Z))
################# ACE #################
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



################# PFA #################
t1 <- Sys.time()
source("./scripts/use.R")
source("./scripts/core.R")
RE2 <- pfa.test(t(Z), t(X), tval= "pval", mat_est="sample", plot = "linear") #
summary(RE2)

# load("./output/PFA_all_results.rda")

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

# save(RE2, file = "./output/PFA_all_results.rda")
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

t1 <- Sys.time()
p_values_adj <- p.adjust(p_values, method = "BH")
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

# save(ACE_result, PFA_result, BH_result, e_BH_result, BY_result, 
#      file = paste0("./output/differential_results_gama_", gama, ".rda"))