################################################################################
############  VennDiagram  ############
################################################################################
rm(list = ls())
library(ggplot2)
library(patchwork)
library(VennDiagram)
library(latex2exp)
library(ComplexUpset)
library(ggplot2)

load("./output/differential_results_gama_0.001.rda")
set1 <- ACE_result$index
set2 <- PFA_result$index
set3 <- BH_result$index
set4 <- e_BH_result$index
set5 <- BY_result$index

# union_ACE_PFA <- union(set1, set2)
# reject_only_by_PFA <- setdiff(union_ACE_PFA, set1) # unique miRNA rejected by PFA (donot consider BH!!!!)
# reject_only_by_ACE <- setdiff(union_ACE_PFA, set2) # unique miRNA rejected by ACE (donot consider BH!!!!)

# 找出唯一在 set1 中的元素
unique_in_ACE <- setdiff(set1, union(union(union(set2, set3), set4), set5))


# 首先获取所有集合的并集
all_elements <- unique(c(set1, set2, set3, set4, set5))
# 创建一个逻辑矩阵
df <- data.frame(
  ACE = all_elements %in% set1,
  PFA = all_elements %in% set2,
  BH = all_elements %in% set3,
  eBH = all_elements %in% set4,  # 使用简单的名称
  BY = all_elements %in% set5,
  row.names = all_elements
)
colnames(df) <- c("ACE", "PFA", "BH", "e-BH", "BY")  # 手动指定列名

upset(
  data = df,
  intersect = colnames(df),
  name = NULL,
  width_ratio = 0.5,
  set_sizes = FALSE,
  themes = upset_modify_themes(
    list(
      "intersections_matrix" = theme(
        axis.text.y = element_text(face = "italic", size = 14)
      )
    )
  ),
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(
        size = 6
      ),
      mapping = aes(fill = "count")
    ) + 
      scale_fill_manual(values = c("count" = "#1868B2")) +
      theme(legend.position = "none",
            axis.text.y = element_text(size = 14),  # 调整纵坐标字体大小
            axis.title.y = element_text(size = 14))  # 调整纵坐标标题字体大小
  )
)



venn.plot <- venn.diagram(
  x = list(ACE = set1, PFA = set2, BH = set3, `e-BH` = set4, BY = set5),
  filename = NULL,
  fill = c("#F94141", "#1868B2", "#018A67", "#F98F34", "#6A4C93"),
  alpha = 0.5,
  cex = 2, # 控制标签大小
  fontface = "plain",
  cat.fontface = "italic",
  margin = 0.01, # 增大或减小边距
  # cat.pos = c(-20, 20, 0), # 调整类别标签的位置
  cat.dist = 0.1, # 调整类别标签距离
  cat.cex = 2, # 调整类别标签大小
  disable.logging = T # Disable log file generation
)

grid.newpage()
grid.draw(venn.plot)








#######################################################################
#############             Volcano Plot                    #############
#######################################################################
#library(EnhancedVolcano)
library(latex2exp)
library(gcookbook)
library(ggplot2)
library(grid)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(dplyr)
library(ACE)

gene_data <- read.csv("./data/data.benchmark.csv") %>% .[, 2]
Z <- read.csv("./data/Vbench.csv")
X <- read.csv("./data/Ebench.csv")

result_ace <- ACE(Z, X, gama = 0.001)
fold_change <- result_ace$Adjusted_mean_difference
adj_pvalues <- 2 * pnorm(abs(result_ace$Adjusted_statistics), lower.tail = F)

# using adjusted p-values in ACE
dat <- data.frame(miRNA = gene_data,
                  mean_diff = fold_change,
                  log_pvalue = -log10(adj_pvalues))
# dat$log_pvalue <- pmin(-log10(adj_pvalues), 60)

## let-7a,d,e, 21, 30c,e
label_index <- c(2001, 1432, 757, 2580, 1383, 547)
label_index %in% unique_in_ACE
dat_label <- dat[label_index, ]
miRNA_name <- dat$miRNA[label_index]

sp <- ggplot(dat, aes(x=mean_diff, y=log_pvalue)) + 
  geom_point(alpha=.4) + 
  scale_x_continuous(expand = c(0,0),limits = c(-4, 4)) +
  scale_y_continuous(expand = c(0,0),limits = c(-10, 100)) +
  labs(x = "Mean expression difference: ovarian-endometrial", 
       y = TeX("$-\\log_{10}(p-value)$")) +
  geom_point(size = 2, color = "red", shape = 2, data = dat_label) +
  ggrepel::geom_label_repel(
    aes(label = miRNA_name),
    data = dat_label,
    color="black",
    size = 5,
    box.padding=unit(1.5, "lines"),
    # label.padding = unit(0.1,"lines"),
    max.overlaps = 20
  ) + 
  theme_bw() + 
  theme(text = element_text(size = 18))
sp








##############################################
#####     FDR vs number of rejection     #####
##############################################
rm(list = ls())
library(ggplot2)
FDR <- rep(c(0.001, 0.005,0.01,0.02,0.03,0.04,0.05), times = 5)
FDR
Rejection <- c(308, 417, 493, 608, 686, 725, 777,
               158, 252, 290, 395, 470, 532, 579,
               101, 123, 151, 185, 222, 251, 262,
               62, 95, 104, 112, 120, 123, 126,
               33, 53, 58, 74, 87, 95, 99)
method <- rep(c('ACE','PFA','BH', 'BY', 'e-BH'), each = 7)
method
df <- data.frame(FDR = FDR, Method = method, Rejection = Rejection)
plottt <- ggplot(df, aes(x = FDR, y = Rejection, linetype = Method, shape = Method)) + 
  geom_line(linewidth=0.7) + geom_point(size=3) + 
  theme_bw(base_size=15) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(55, "pt"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13, face="italic")
  ) + 
  ylab("Number of rejection") + xlab("FDR control level")
plottt




##############################################
#####   Hist and Genomic control factor  #####
##############################################
# median((result_ace$Adjusted_statistics)^2)/qchisq(0.5, df=1)
rm(list = ls())
Z <- read.csv("./data/Vbench.csv")
X <- read.csv("./data/Ebench.csv")
library(ACE)
library(ggplot2)
ace_result <- ACE(Z, X)
statistics <- ace_result$Adjusted_statistics
pvalues <- 2 * pnorm(abs(statistics), lower.tail = F)
pval <- data.frame(pvalues = pvalues)

p_gg <- ggplot(pval, aes(x = pvalues)) +
  geom_histogram(color="black", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
  theme_bw() +
  theme(text = element_text(size = 10),  # 设置字体大小
        plot.title = element_text(size = 15, hjust = 0.5),  # 调整标题字体大小
        axis.title.x = element_blank(),  # 去掉 x 轴标题
        axis.title.y = element_blank(),  # 去掉 y 轴标题
        axis.text.x = element_text(size = 10),  # 设置坐标轴刻度字体大小
        axis.text.y = element_text(size = 10)) +
  xlim(0, 1) +
  annotate("text", x = 0.95, y = 1000, label = "GCF = 1.99", size = 7, hjust = 1)
p_gg

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






##############################################
#####        Unique reject by ACE        #####
##############################################
rm(list = ls())
library(dplyr)
library(ACE)
gene_name <- read.csv("./data/data.benchmark.csv") %>% .[, 2]
Z <- read.csv("./data/Vbench.csv")
X <- read.csv("./data/Ebench.csv")
load("./output/differential_results_gama_0.001.rda")
set1 <- ACE_result$index
set2 <- PFA_result$index
set3 <- BH_result$index
set4 <- e_BH_result$index
set5 <- BY_result$index
# 找出唯一在 set1 中的元素
unique_in_ACE <- setdiff(set1, union(union(union(set2, set3), set4), set5))

result_ace <- ACE(Z, X, gama = 0.001)

miRNA_names <- gene_name
adj_pvalues <- 2 * pnorm(abs(result_ace$Adjusted_statistics), lower.tail = F)
mu_hat <- result_ace$Adjusted_mean_difference
statistic <- result_ace$Adjusted_statistics
sd <- sqrt(min(ncol(Z), ncol(X)))*mu_hat / statistic
regulate <- ifelse(statistic > 0, "Up", "Down")

# Combine into a data frame
result_df <- data.frame(
  miRNA_reject = miRNA_names,
  adj_pvalues = adj_pvalues,
  mu_hat = mu_hat,
  sd = sd,
  statistic = statistic,
  regulate = regulate
) %>% .[unique_in_ACE, ]

write.csv(result_df, file = "./output/ACE_unique_miRNA_results.csv", row.names = FALSE)
head(result_df)


# Create a version with unique miRNAs
result_df_unduplicated <- result_df %>% 
  distinct(miRNA_reject, .keep_all = TRUE) %>%
  arrange(desc(statistic))
write.csv(result_df_unduplicated, file = "./output/ACE_unduplicated_miRNA_results.csv", row.names = FALSE)


result_all <- data.frame(
  miRNA_reject = miRNA_names,
  adj_pvalues = adj_pvalues,
  mu_hat = mu_hat,
  sd = sd,
  statistic = statistic,
  regulate = regulate
) %>% .[set1, ]
write.csv(result_all, file = "./output/ACE_all_reject_miRNA_results.csv", row.names = FALSE)






##############################################
#####               KEGG                 #####
##############################################
rm(list = ls())
# BiocManager::install("miRBaseConverter")

library(clusterProfiler)
library(org.Hs.eg.db)
library(miRBaseConverter)
library(multiMiR)
library(dplyr)
library(ggplot2)
data <- read.csv("./output/ACE_unduplicated_miRNA_results.csv")
mirna_list <- data$miRNA_reject
# 修复miRNA名称中的let-7家族
mirna_list <- gsub("hsa-let-(7[a-z])", "hsa-miR-let-\\1", mirna_list)
mirna_list[grep("hsa-miR-let-7", mirna_list)]

# 使用multiMiR获取已验证的靶基因
targets <- get_multimir(
  org     = "hsa",         # 人类
  mirna   = mirna_list,
  table   = "validated",   # 选择已验证的相互作用
  summary = TRUE           # 汇总结果
)

# 提取靶基因Symbol并去重
target_genes <- unique(targets@data$target_symbol)


# 转换Symbol到Entrez ID
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys      = target_genes,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)

# 移除NA值并去重
entrez_ids <- na.omit(entrez_ids)
entrez_list <- unique(as.character(entrez_ids))



# 执行KEGG富集分析
kegg_result <- enrichKEGG(
  gene          = entrez_list,
  organism      = "hsa",         # 人类
  keyType       = "kegg",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

# 提取结果并保存为数据框
kegg_df <- as.data.frame(kegg_result)
head(kegg_df)



### 7. 添加miRNA对应列
# 从multiMiR结果中提取miRNA-靶基因对应关系
mirna_target_map <- targets@data %>% 
  select(mature_mirna_id, target_symbol) %>% 
  distinct()  # 去重
save(mirna_target_map, file = "./output/mirna_target_map.rda")

# 创建函数：根据基因Symbol查找对应的miRNA
find_mirnas_for_pathway <- function(gene_symbols) {
  # 分割gene_symbols字符串
  symbols <- unlist(strsplit(gene_symbols, "/"))
  
  # 查找对应的miRNA
  related_mirnas <- mirna_target_map %>%
    filter(target_symbol %in% symbols) %>%
    pull(mature_mirna_id) %>%
    unique() %>%
    sort()
  
  return(paste(related_mirnas, collapse = "/"))
}

# 添加miRNA列到KEGG结果
kegg_df <- kegg_df %>%
  mutate(
    # 将geneID拆分为单独的Entrez ID
    GeneSymbols = lapply(strsplit(geneID, "/"), function(ids) {
      # 对每个Entrez ID调用mapIds，使用multiVals="first"获取第一个Symbol
      symbols <- mapIds(org.Hs.eg.db,
                        keys = ids,
                        keytype = "ENTREZID",
                        column = "SYMBOL",
                        multiVals = "first")
      # 返回非NA的Symbol并用"/"连接
      paste(na.omit(symbols), collapse = "/")
    }),
    # 查找miRNA对应关系
    miRNAs = sapply(GeneSymbols, find_mirnas_for_pathway)
  )




# 查看添加miRNA后的结果
head(kegg_df[, c("Description", "miRNAs", "geneID")], 3)

# 保存完整结果
# write.table(as.data.frame(kegg_df), "./output/KEGG_results_with_miRNA.txt")
save(kegg_df, file = "./output/KEGG_results_with_miRNA.rda")
# 按调整后p值排序，取前10个通路

selected_pathways <- c("Signal transduction", "Cancer: specific types", 
                       "Cell growth and death", "Neurodegenerative disease", 
                       "Replication and repair", "Endocrine and metabolic disease", 
                       "Aging", "Drug resistance: antineoplastic", 
                       "Chromosome", "Global and overview maps")

# 使用 dplyr 包来筛选数据
library(dplyr)

# 筛选出包含前10个通路的行
filtered_kegg_df <- kegg_df %>% 
  filter(subcategory %in% selected_pathways)

# 查看筛选后的数据
print(filtered_kegg_df)

filtered_kegg_df <- filtered_kegg_df %>%
  group_by(subcategory) %>%  # 按照通路（子类别）分组
  slice_min(order_by = p.adjust, n = 1) %>%  # 按照 p.adjust 列的最小值选择每组的第一行
  ungroup()  # 取消分组，确保结果是一个普通的数据框

top10 <- filtered_kegg_df %>%
  arrange(p.adjust) %>%
  head(10)


top10$log_pvalue <- -log10(top10$pvalue)
top10$Description <- factor(top10$Description, 
                                levels = top10$Description[order(top10$log_pvalue, decreasing = F)])

ggplot(top10, aes(x = log_pvalue, y = Description)) +
  geom_bar(stat = "identity", fill = "#F94141") +
  labs(x = "-log10(p-value)", y = NULL, title = NULL) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),  # 增加x轴和y轴刻度字体大小
    axis.title = element_text(size = 14), # 增加轴标题字体大小
    panel.grid.major.y = element_blank(), # 去掉y轴网格线
    panel.grid.minor.y = element_blank()  # 去掉y轴次级网格线
  )





# kegg_df 中提取 log_pvalue 最大的前20个条目
load(file = "./output/KEGG_results_with_miRNA.rda")

top20 <- kegg_df %>%
  mutate(log_pvalue = -log10(pvalue)) %>%
  arrange(desc(log_pvalue)) %>%  # 按 log_pvalue 降序排列
  head(20) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))  # 确保绘图顺序正确

ggplot(top20, aes(x = log_pvalue, y = Description)) +
  geom_bar(stat = "identity", fill = "#F94141") +
  labs(x = "-log10(p-value)", y = NULL, title = NULL) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
