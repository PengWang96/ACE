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

load("./output/differential_results_gama_0.1.rda")
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
load("./output/differential_results_gama_0.1.rda")
Z <- read.csv("./data/expr_treatment.csv", row.names = 1, stringsAsFactors = FALSE)
X <- read.csv("./data/expr_control.csv", row.names = 1, stringsAsFactors = FALSE)
ENSG_ID <- as.data.frame(rownames(Z)); colnames(ENSG_ID) <- "ENSG_ID"
ENSG_gene_symbols <- read.csv("./data/gene_symbols.csv", 
                              stringsAsFactors = FALSE)

a <- as.data.frame(table(ENSG_gene_symbols$ENSEMBL))
intersect(as.numeric(rownames(a[which(a$Freq>1),])), ACE_result$index)
# numeric(0) great
ENSG_gene_symbols <- ENSG_gene_symbols[!duplicated(ENSG_gene_symbols$ENSEMBL), ]

identical(ENSG_gene_symbols$ENSEMBL, ENSG_ID$ENSG_ID)
# TRUE, great

                                                             
# 计算 Fold Change 和 p-value
fold_change <- numeric(nrow(ENSG_ID))
p_values <- numeric(nrow(ENSG_ID))
for (i in 1:nrow(ENSG_ID)) {
  pos_expr <- as.numeric(Z[i, ])
  neg_expr <- as.numeric(X[i, ])
  fold_change[i] <- (mean(pos_expr) - mean(neg_expr))

  t_test_result <- t.test(pos_expr, neg_expr)
  p_values[i] <- t_test_result$p.value
}


# using adjusted p-values in ACE
dat <- data.frame(gene = ENSG_gene_symbols$SYMBOL,
                  mean_diff = fold_change,
                  log_pvalue = -log10(p_values))
dat[PFA_result$index, ]
10^(-dat[PFA_result$index, ]$log_pvalue)




label_index <- ACE_result$index
dat_label <- na.omit(dat[label_index, ])

selected_gene <- c("EIF5A", "ACTN3", "CNIH2", "BCL6")
dat_label <- dat_label[match(selected_gene, dat_label$gene), ]

sp <- ggplot(dat, aes(x=mean_diff, y=log_pvalue)) + 
  geom_point(alpha=.4) + 
  scale_x_continuous(expand = c(0,0),limits = c(-1.5, 1.5)) +
  scale_y_continuous(expand = c(0,0),limits = c(-0.5, 5.5)) +
  labs(x = "Mean expression difference: mutation carriers-controls", 
       y = TeX("$-\\log_{10}(p-value)$")) +
  geom_point(size = 2, color = "red", shape = 2, data = dat_label) +
  ggrepel::geom_label_repel(
    aes(label = gene),
    data = dat_label,
    color="black",
    size = 5,
    box.padding=unit(1.5, "lines"),
    # label.padding = unit(0.1,"lines"),
    max.overlaps = 5
  ) +
  theme_bw() + 
  theme(text = element_text(size = 18))
sp





################################################################################
############   相关性热图   ############
################################################################################
rm(list = ls())
X <- read.csv("./data/expr_treatment.csv", row.names = 1, stringsAsFactors = FALSE)
Z <- read.csv("./data/expr_control.csv", row.names = 1, stringsAsFactors = FALSE)

p <- nrow(Z); n1 <- ncol(Z); n2 <- ncol(X)
Y <- matrix(0, p, n1)
for (jy in 1:n1){
  Y[,jy] <- Z[,jy] - sqrt(n1/n2)*X[,jy] + apply(X[,1:n1], 1, sum)/sqrt(n1*n2) - apply(X,1,mean)
}

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
result_ace <- ACE(Z, X, gama = 0.1)
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
FDR <- rep(seq(0.05, 0.14, 0.01), times = 5)
FDR
Rejection <- c(26, 35, 37, 40, 51, 63, 74, 78, 87, 91,
               2, 3, 4, 7, 8, 8, 8, 8, 9, 13,
               c(rep(0, 10)),
               rep(0, 10),
               rep(0, 10))
method <- rep(c('ACE','PFA','BH', 'BY', 'e-BH'), each = 10)
method
df <- data.frame(FDR = FDR, Method = method, Rejection = Rejection)

plottt <- ggplot(df, aes(x = FDR, y = Rejection, linetype = Method, shape = Method)) +
  geom_line(linewidth=0.7) + geom_point(size=3) +
  theme_bw(base_size=15) +
  theme(legend.position = "inside",
        legend.position.inside = c(.11, .88),
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





##############################################
#####   Hist and Genomic control factor  #####
##############################################
# median((ace_result$Adjusted_statistics)^2)/qchisq(0.5, df=1)
rm(list = ls())
Z <- read.csv("./data/expr_treatment.csv", row.names = 1, stringsAsFactors = FALSE)
X <- read.csv("./data/expr_control.csv", row.names = 1, stringsAsFactors = FALSE)
library(ACE)
library(ggplot2)
ace_result <- ACE(Z, X)
median((ace_result$Adjusted_statistics)^2)/qchisq(0.5, df=1)

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
  ylim(0, 2000) +
  annotate("text", x = 0.95, y = 1700, label = "GCF = 1.41", size = 7, hjust = 1)
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
ENSG_gene_symbols <- read.csv("./data/gene_symbols.csv", 
                              stringsAsFactors = FALSE)
ENSG_gene_symbols <- ENSG_gene_symbols[!duplicated(ENSG_gene_symbols$ENSEMBL), ]
Z <- read.csv("./data/expr_treatment.csv", row.names = 1, stringsAsFactors = FALSE)
X <- read.csv("./data/expr_control.csv", row.names = 1, stringsAsFactors = FALSE)
load("./output/differential_results_gama_0.1.rda")
set1 <- ACE_result$index
set2 <- PFA_result$index
set3 <- BH_result$index
set4 <- e_BH_result$index
set5 <- BY_result$index
# 找出唯一在 set1 中的元素
unique_in_ACE <- setdiff(set1, union(union(union(set2, set3), set4), set5))

result_ace <- ACE(Z, X, gama = 0.1)

adj_pvalues <- 2 * pnorm(abs(result_ace$Adjusted_statistics), lower.tail = F)
mu_hat <- result_ace$Adjusted_mean_difference * (-1)
statistic <- result_ace$Adjusted_statistics * (-1)
sd <- sqrt(min(ncol(Z), ncol(X)))*mu_hat / statistic
regulate <- ifelse(statistic > 0, "Up", "Down")

# Combine into a data frame
result_df <- data.frame(
  gene_reject = ENSG_gene_symbols$SYMBOL,
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


result_df_unduplicated$gene_reject



##############################################
#####              GO                    #####
##############################################
rm(list = ls())
library(ggplot2)
library(dplyr)

# 读取数据并添加类别信息
BP <- read.table("./output/GO_BP.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
BP$Type <- "Biological Process"
MF <- read.table("./output/GO_MF.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
MF$Type <- "Molecular Function"
CC <- read.table("./output/GO_CC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
CC$Type <- "Cellular Component"
data <- rbind(BP, MF, CC)

# 处理通路名称和p值转换
data$Pathway <- gsub("^GO:\\d+~", "", data$Term)  # 提取纯通路名称
data$log_pvalue <- -log10(data$PValue)

# 全局排序（不再按类别分组排序）
data <- data %>%
  arrange(desc(log_pvalue)) %>%  # 按p值降序排列
  mutate(
    Pathway = factor(Pathway, levels = rev(unique(Pathway)))  # 反转因子确保图形顶部显示最大p值
  )

# 绘制条形图
ggplot(data, aes(x = log_pvalue, y = Pathway, fill = Type)) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(x = "-log10(p-value)", y = "") +
  scale_fill_manual(values = c("#F94144", "#577590", "#43AA8B")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.grid.major.y = element_blank()
  )
















head(filtered_data)
# Example gene list to check
genes_to_check <- c("KIF1B", "AURKB", "EPHA5", "HAND1", 
                    "CDKN2A", "CDKN1B", "PTEN", "STMN1", "MAPK3")

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
