rm(list = ls())
library(limma)
library(dplyr)
library(stringr)
library(org.Hs.eg.db)
# 步骤1：读取数据 --------------------------------------------------------------
# 读取标准化矩阵
expr <- read.table(file = './data/GSE51799_Normalised_data.txt',
                   sep = '\t',
                   header = T,
                   fill = T,           # If one line in the file has less data than the other rows, a blank field is automatically added.  
                   stringsAsFactors = F, #  The string is not changed to factor
                   check.names = FALSE)
# 查看数据结构
expr[1:3, 1:5]

############# Screen probes ###############
exprSet <- avereps(expr[, -1],         # If multiple probes correspond to one gene, taking mean
                   ID = expr$Gene_ID)
exprSet <- as.data.frame(exprSet)


ensg_ids <- rownames(exprSet)
gene_symbols <- select(org.Hs.eg.db,
                       keys = ensg_ids,
                       columns = c("SYMBOL"),
                       keytype = "ENSEMBL")
head(gene_symbols)
write.csv(gene_symbols, file = "./data/gene_symbols.csv", row.names = F)
# 
# gene_symbols <- gene_symbols[!duplicated(gene_symbols$ENSEMBL), ]
# gene_symbols <- na.omit(gene_symbols)
# gene_symbols <- gene_symbols[!duplicated(gene_symbols$SYMBOL), ]
# 
# exprSet <- exprSet[gene_symbols$ENSEMBL, ]
# 
# identical(rownames(exprSet), gene_symbols$ENSEMBL)
# rownames(exprSet) <- gene_symbols$SYMBOL

# 步骤3：样本分组 -------------------------------------------------------------
expr_control <- exprSet[, 1:33]
expr_treatment <- exprSet[, -(1:33)]

write.csv(expr_control, file = "./data/expr_control.csv", row.names = TRUE)
write.csv(expr_treatment, file = "./data/expr_treatment.csv", row.names = TRUE)
