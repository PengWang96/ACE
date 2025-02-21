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

load("./output/differential_results_gama_0.05.rda")
set1 <- ACE_result$index
set2 <- PFA_result$index
set3 <- BH_result$index
set4 <- e_BH_result$index
set5 <- BY_result$index

union_ACE_PFA <- union(set1, set2)
reject_only_by_PFA <- setdiff(union_ACE_PFA, set1) # unique gene rejected by PFA (donot consider BH!!!!)
reject_only_by_ACE <- setdiff(union_ACE_PFA, set2) # unique gene rejected by ACE (donot consider BH!!!!)

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





################################################################################
############  Vocalno Plot  ############
################################################################################
Z <- read.csv("./data/positive_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)
X <- read.csv("./data/negative_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)
probe_gene <- read.csv("data/probe_gene_1076_in_10707.csv")
# 计算 Fold Change 和 p-value
genes <- rownames(Z)
gene_names <- probe_gene$gene
fold_change <- numeric(length(genes))
p_values <- numeric(length(genes))

for (i in seq_along(genes)) {
  # 提取每个基因在 Positive 和 Negative 样本中的表达量
  pos_expr <- as.numeric(Z[i, ])
  neg_expr <- as.numeric(X[i, ])
  
  # 计算平均表达值并求 log2 fold change
  # fold_change[i] <- log2(mean(pos_expr) / mean(neg_expr))
  fold_change[i] <- (mean(pos_expr) - mean(neg_expr))
  
  # 使用 t-test 计算 p-value
  t_test_result <- t.test(pos_expr, neg_expr)
  p_values[i] <- t_test_result$p.value
}

# 创建数据框用于绘制火山图
volcano_data <- data.frame(
  gene = genes,
  gene_name = gene_names,
  fold_change = fold_change,
  p_value = p_values,
  neg_log_pval = -log10(p_values)
)




volcano_data_PFA <- volcano_data
volcano_data_PFA$group <- ifelse(volcano_data_PFA$gene %in% reject_only_by_PFA, "Rejected only by PFA", "Other")
volcano_data_PFA$group <- factor(volcano_data_PFA$group, levels = c("Other", "Rejected only by PFA"))
volcano_data_PFA <- volcano_data_PFA[order(volcano_data_PFA$group == "Rejected only by PFA"), ]
# 绘制火山图并标记基因
ggplot(volcano_data_PFA, aes(x = fold_change, y = neg_log_pval)) +
  geom_point(aes(color = group), alpha = 1) +
  scale_color_manual(values = c("Rejected only by PFA" = "#E66101", "Other" = "black"),
                     labels = c("Rejected only by PFA" = expression("Rejected only by " * italic("PFA")),
                                "Other" = "Other")) +
  labs(title = NULL, x = "Mean expression difference: positive-negative group", y = TeX("$-\\log_{10}(p-value)$")) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.14, .95),
    legend.title = element_blank(),  # 去掉图例的标题
    text = element_text(size = 16)   # 增大整体字体大小
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.5)))  # 增大图例图标的大小




volcano_data_ACE <- volcano_data
volcano_data_ACE$group <- ifelse(volcano_data_ACE$gene %in% reject_only_by_ACE, "Rejected by ACE but not by PFA", "Other")
volcano_data_ACE$group <- factor(volcano_data_ACE$group, levels = c("Other", "Rejected by ACE but not by PFA"))
volcano_data_ACE <- volcano_data_ACE[order(volcano_data_ACE$group == "Rejected by ACE but not by PFA"), ]
# 绘制火山图并标记基因
sp <- 
ggplot(volcano_data_ACE, aes(x = fold_change, y = neg_log_pval)) +
  geom_point(aes(color = group), alpha = 1) +
  scale_color_manual(values = c("Rejected by ACE but not by PFA" = "#6186ad", "Other" = "black"),
                     labels = c("Rejected by ACE but not by PFA" = expression("Rejected by " * italic("ACE") * " but not by " * italic("PFA")),
                                "Other" = "Other")) +
  labs(title = NULL, x = "Mean expression difference: positive-negative group", y = TeX("$-\\log_{10}(p-value)$"))
sp


select_genes <- c("PTEN", "HAND1", "KIF1B", "CDKN2A", "EPHA5", "MAPK3", "STMN1", "AURKB", "CDKN1B")
indices <- match(select_genes, volcano_data_ACE$gene_name)

dat_label <- volcano_data_ACE[indices, ]
sp1 <- sp +
  geom_point(size = 2, color = "red", shape = 2, stroke = 1.2, data = dat_label) +
  # geom_text(size = 10) + 
  ggrepel::geom_label_repel(
    aes(label = select_genes),
    data = dat_label,
    color="black",
    size = 5,
    box.padding=unit(1.5, "lines"),
    # nudge_x = 0.15,  # Adjust horizontal nudge distance
    nudge_y = 0.3,  # Adjust vertical nudge distance
    # direction = "y",  # Allow labels to repel both left and right
    # force = 10,  # Increase repelling force for further label separation
    # label.padding = unit(0.1,"lines"),
    max.overlaps = 30
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.19, .95),
    legend.title = element_blank(),  # 去掉图例的标题
    text = element_text(size = 16)   # 增大整体字体大小
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.5)))  # 增大图例图标的大小
sp1






################################################################################
############  数据的分布  ############
################################################################################
rm(list = ls())
setwd("D:/R/factor model20230625/real data E-TABM-38")
Z <- read.csv("./data/positive_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)
X <- read.csv("./data/negative_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)

# p <- nrow(Z); n1 <- ncol(Z); n2 <- ncol(X)
# Y <- matrix(0, p, n1)
# for (jy in 1:n1){
#   Y[,jy] <- Z[,jy] - sqrt(n1/n2)*X[,jy] + apply(X[,1:n1], 1, sum)/sqrt(n1*n2) - apply(X,1,mean)
# }
# Y_bar <- rowMeans(Y)

Z_bar <- rowMeans(Z)
X_bar <- rowMeans(X)
ggplot(data.frame(Z_bar = Z_bar), aes(x = Z_bar)) +
  geom_histogram(color = "white", fill = "#6186ad", bins = 20, boundary = 0, alpha = 0.8) +
  # scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        text = element_text(size = 16)) +
  xlab("Data in positive group") +
  ylab("Count")

ggplot(data.frame(X_bar = X_bar), aes(x = X_bar)) +
  geom_histogram(color = "white", fill = "#6186ad", bins = 20, boundary = 0, alpha = 0.8) +
  # scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        text = element_text(size = 16)) +
  xlab("Data in negative group") +
  ylab("Count")





hist(Z, breaks = 30, col = "skyblue")
hist(X, breaks = 30, col = "skyblue")
# install.packages("moments")
library(moments)
skewness(Z)
skewness(X)

















################################################################################
############   相关性热图   ############
################################################################################
rm(list = ls())
setwd("D:/R/factor model20230625/real data E-TABM-38")
Z <- read.csv("./data/positive_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)
X <- read.csv("./data/negative_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)

p <- nrow(Z); n1 <- ncol(Z); n2 <- ncol(X)
Y <- matrix(0, p, n1)
for (jy in 1:n1){
  Y[,jy] <- Z[,jy] - sqrt(n1/n2)*X[,jy] + apply(X[,1:n1], 1, sum)/sqrt(n1*n2) - apply(X,1,mean)
}


# load("./output/differential_results.rda")
# set1 <- ACE_result$index
# set2 <- PFA_result$index
# set3 <- BH_result$index
# inter <- intersect(intersect(set1, set2), set3)
# 计算变量之间的相关性矩阵
correlation_matrix <- cor(t(Y[1:500, ]))

# 安装并加载 pheatmap 包（如果尚未安装）
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)

# 定义自定义颜色和区间
breaks <- c(-1, -0.6, -0.2, 0.2, 0.6, 1)  # 定义区间
# colors <- c("#053061", "#4393C3", "white", "#ffc0cb", "#E60000")
colors <- c("blue", "lightblue", "white", "pink", "red")  # 定义对应颜色

# 绘制相关性热图
pheatmap(correlation_matrix,
         cluster_rows = T,   # 行聚类
         cluster_cols = T,   # 列聚类
         treeheight_row = 0,     # 不显示行聚类树
         treeheight_col = 0,     # 不显示列聚类树
         show_rownames = FALSE, # 不显示行名
         show_colnames = FALSE, # 不显示列名
         color = colors,        # 使用自定义颜色
         breaks = breaks        # 使用自定义区间
)




library(ACE)
result_ace <- ACE(Z, X, gama = 0.05)
B_hat <- result_ace$Estimate_factor_loadings
f_hat <- result_ace$Estimate_factors
Y_adjust <- Y - B_hat %*% f_hat
correlation_matrix <- cor(t(Y_adjust[1:500, ]))
pheatmap(correlation_matrix,
         cluster_rows = T,   # 行聚类
         cluster_cols = T,   # 列聚类
         treeheight_row = 0,     # 不显示行聚类树
         treeheight_col = 0,     # 不显示列聚类树
         show_rownames = FALSE, # 不显示行名
         show_colnames = FALSE, # 不显示列名
         color = colors,        # 使用自定义颜色
         breaks = breaks        # 使用自定义区间
)







##############################################
#####     FDR vs number of rejection     #####
##############################################
rm(list = ls())
library(ggplot2)
FDR <- rep(seq(0.01, 0.1, 0.01), times = 5)
FDR
Rejection <- c(1147, 1446, 1620, 1785, 1932, 2078, 2196, 2304, 2431, 2537,
               851, 1094, 1310, 1472, 1611, 1743, 1877, 1995, 2089, 2197,
               299, 397, 492, 569, 649, 703, 770, 835, 875, 955,
               3, 122, 163, 183, 219, 236, 261, 272, 293, 304,
               rep(0, 10))
method <- rep(c('ACE','PFA','BH', 'BY', 'e-BH'), each = 10)
method
df <- data.frame(FDR = FDR, Method = method, Rejection = Rejection)
# plottt <- ggplot(df, aes(x = FDR, y = Rejection, linetype = Method, shape = Method)) + 
#   geom_line(linewidth=0.7) + geom_point(size=3) + 
#   theme_bw(base_size=15) + 
#   theme(legend.position = "bottom",
#         legend.key.width = unit(55, "pt"),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 15, face="italic")
#   ) + 
#   ylab("Number of rejection") + xlab("FDR control level")
# plottt

plottt <- ggplot(df, aes(x = FDR, y = Rejection, linetype = Method, shape = Method)) +
  geom_line(linewidth=0.7) + geom_point(size=3) +
  theme_bw(base_size=15) +
  theme(legend.position = "inside",
        legend.position.inside = c(.11, .87),
        legend.key.width = unit(55, "pt"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, face="italic"),
        text = element_text(size = 16)   # 增大整体字体大小
  ) +
  ylab("Number of rejection") + xlab("FDR control level")
plottt



# library(pheatmap)
# Z <- read.csv("./data/positive_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)
# X <- read.csv("./data/negative_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)
# # 假设 PFA 识别的基因索引在 set2 中
# pfa_genes <- rownames(Z)[unique_in_ACE]
# 
# # 提取 PFA 基因在两组样本中的表达矩阵
# combined_expr <- cbind(Z[pfa_genes, ], X[pfa_genes, ])
# 
# # 标注样本组
# sample_groups <- c(rep("Positive", ncol(Z)), rep("Negative", ncol(X)))
# 
# 
# annotation_col <- data.frame(Group = sample_groups)
# rownames(annotation_col) <- colnames(combined_expr)
# 
# # 绘制热图
# p2 <- 
#   pheatmap(
#     combined_expr, 
#     cluster_rows = FALSE, 
#     cluster_cols = T, 
#     scale = "row",                # 对行（基因）进行归一化
#     annotation_col = annotation_col,  # 样本的分组信息
#     show_rownames = FALSE,        # 隐藏基因名
#     show_colnames = FALSE,          # 显示样本名
#     color = colorRampPalette(c("blue", "white", "red"))(100)  # 设置更显著的颜色梯度
#   )
# print(p2)





rm(list = ls())
# 读取数据
Z <- read.csv("./data/positive_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)
X <- read.csv("./data/negative_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)

# 定义用于存储 t 统计量的向量
t_statistics <- c()

# 对每一行（每个基因）进行 t 检验
for (i in 1:nrow(Z)) {
  # 从正负样本中提取第 i 行
  Z_gene <- as.numeric(Z[i, ])
  X_gene <- as.numeric(X[i, ])
  
  # 进行双样本 t 检验
  t_test_result <- t.test(Z_gene)
  
  # 存储 t 统计量
  t_statistics <- c(t_statistics, t_test_result$statistic)
}

# 绘制 t 统计量的分布图
hist(t_statistics, breaks = 30, col = "skyblue", probability = TRUE, 
     main = "T-statistic Distribution vs Normal Distribution", 
     xlab = "T-statistic", ylab = "Density")

# 生成用于叠加正态分布的 x 轴值范围
x_values <- seq(min(t_statistics), max(t_statistics), length.out = 100)

# 计算标准正态分布的密度
normal_density <- dnorm(x_values)

# 叠加标准正态分布曲线
lines(x_values, normal_density, col = "red", lwd = 2)

# 添加图例
legend("topright", legend = c("T-statistic", "Normal Distribution"), 
       col = c("skyblue", "red"), lwd = 2)





##############################################
#####   Hist and Genomic control factor  #####
##############################################
# median((result_ace$Adjusted_statistics)^2)/qchisq(0.5, df=1)
rm(list = ls())
Z <- read.csv("./data/positive_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)
X <- read.csv("./data/negative_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)
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
  annotate("text", x = 0.95, y = 3000, label = "GCF = 3.21", size = 7, hjust = 1)
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
Z <- read.csv("./data/positive_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)
X <- read.csv("./data/negative_expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)

load("./output/differential_results_gama_0.05.rda")
set1 <- ACE_result$index
set2 <- PFA_result$index
set3 <- BH_result$index
set4 <- e_BH_result$index
set5 <- BY_result$index
# 找出唯一在 set1 中的元素
unique_in_ACE <- setdiff(set1, union(union(union(set2, set3), set4), set5))

result_ace <- ACE(Z, X, gama = 0.05)

genes <- read.csv("data/probe_gene_1076_in_10707.csv")
gene_names <- genes$gene
adj_pvalues <- 2 * pnorm(abs(result_ace$Adjusted_statistics), lower.tail = F)
mu_hat <- result_ace$Adjusted_mean_difference
statistic <- result_ace$Adjusted_statistics
sd <- sqrt(min(ncol(Z), ncol(X)))*mu_hat / statistic
regulate <- ifelse(statistic > 0, "Up", "Down")

# Combine into a data frame
result_df <- data.frame(
  gene_reject = gene_names,
  adj_pvalues = adj_pvalues,
  mu_hat = mu_hat,
  sd = sd,
  statistic = statistic,
  regulate = regulate
) %>% .[unique_in_ACE, ]

write.csv(result_df, file = "./output/ACE_unique_gene_results.csv", row.names = FALSE)
head(result_df)


# Create a version with unique genes
result_df_unduplicated <- result_df %>% 
  distinct(gene_reject, .keep_all = TRUE) %>%
  arrange(desc(statistic))
write.csv(result_df_unduplicated, file = "./output/ACE_unduplicated_gene_results.csv", row.names = FALSE)


result_all <- data.frame(
  gene_reject = gene_names,
  adj_pvalues = adj_pvalues,
  mu_hat = mu_hat,
  sd = sd,
  statistic = statistic,
  regulate = regulate
) %>% .[set1, ]
write.csv(result_all, file = "./output/ACE_all_reject_gene_results.csv", row.names = FALSE)








##############################################
#####           KEGG Top 20              #####
##############################################
################################################################################
# kegg_df 中提取 log_pvalue 最大的前20个条目
rm(list = ls())
library(ggplot2)
library(dplyr)
data <- read.table("./output/KEGG.txt", 
                   header = TRUE, sep = "\t", 
                   stringsAsFactors = FALSE)

data$log_pvalue <- -log10(data$PValue)
# 提取通路名称（去掉hsa前缀）
data$Pathway <- gsub("^hsa[0-9]+:", "", data$Term)


top20 <- data %>%
  mutate(log_pvalue = -log10(PValue)) %>%
  arrange(desc(log_pvalue)) %>%  # 按 log_pvalue 降序排列
  head(20) %>%
  mutate(Pathway = factor(Pathway, levels = rev(unique(Pathway))))  # 确保绘图顺序正确

ggplot(top20, aes(x = log_pvalue, y = Pathway)) +
  geom_bar(stat = "identity", fill = "#F94141") +
  labs(x = "-log10(p-value)", y = NULL, title = NULL) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )






##############################################
#####              KEGG                  #####
##############################################
rm(list = ls())
if (!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)
library(ggplot2)
data <- read.table("./output/KEGG.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 需要的通路名称
pathways <- c("MAPK signaling pathway", "Dopaminergic synapse", "Cholinergic synapse", 
              "Axon guidance", "Neurotrophin signaling pathway", "VEGF signaling pathway", 
              "Glutamatergic synapse", "GABAergic synapse", "Serotonergic synapse", 
              "Retrograde endocannabinoid signaling")

filtered_data <- subset(data, grepl(paste(pathways, collapse = "|"), Term))
filtered_data$log_pvalue <- -log10(filtered_data$PValue)
# 提取通路名称（去掉hsa前缀）
filtered_data$Pathway <- gsub("^hsa[0-9]+:", "", filtered_data$Term)

# 重新排列 Pathway 根据 log_pvalue
filtered_data$Pathway <- factor(filtered_data$Pathway, 
                                levels = filtered_data$Pathway[order(filtered_data$log_pvalue, decreasing = F)])

# 绘制柱状图并调整字体大小
ggplot(filtered_data, aes(x = log_pvalue, y = Pathway)) +
  geom_bar(stat = "identity", fill = "#F94141") +
  labs(x = "-log10(p-value)", y = NULL, title = NULL) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),  # 增加x轴和y轴刻度字体大小
    axis.title = element_text(size = 14), # 增加轴标题字体大小
    panel.grid.major.y = element_blank(), # 去掉y轴网格线
    panel.grid.minor.y = element_blank()  # 去掉y轴次级网格线
  )







head(filtered_data)
# Example gene list to check
genes_to_check <- c("KIF1B", "AURKB", "EPHA5", "HAND1", "CDKN2A", "CDKN1B", "PTEN", "STMN1", "MAPK3")

# Function to find pathways for a list of genes
find_pathways <- function(genes, data) {
  gene_pathways <- list()
  
  for (gene in genes) {
    pathways <- data$Pathway[sapply(data$Genes, function(x) grepl(gene, x))]
    
    # If the gene is found in any pathway, store those pathways; otherwise, store NA
    if (length(pathways) > 0) {
      gene_pathways[[gene]] <- paste(pathways, collapse = ", ")
    } else {
      gene_pathways[[gene]] <- NA
    }
  }

  return(gene_pathways)
}

# Apply the function to find pathways for the specified genes
gene_pathways_result <- find_pathways(genes_to_check, filtered_data)

# Print the results
gene_pathways_result
