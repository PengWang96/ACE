rm(list = ls())
library(ACE)
library(mvtnorm)
library(ggplot2)
library(gridExtra)
p <- 2000; n <- 600; K <- 3; lb <- 0.2; ub <- 0.2; v <- 4
PP <- 0.2              
factor_indicator <- 1
sqrt_factor <- sqrt(v / (v - 2))

berlii <- rbinom(p, 1, PP)
index0 <- which(berlii == 0)
index1 <- which(berlii == 1)
true_pai1 <- length(index1) / p
mu <- matrix(rep(0, 1 * p), nrow = p)
mu[index1] <- runif(length(index1), min = lb, max = ub)

X <- matrix(rep(0, n * p), nrow = p)
B <- factor_indicator * matrix(runif(K*p, min=-1, max=1), nrow=p)

sigma <- diag(p)
sigma[index0,index0] <- diag(runif(length(index0), 1, 1))
if (length(index1) == 1) {
  sigma[index1,index1] <- 2 + mu[index1]
} else {
  sigma[index1,index1] <- diag(c(2 + mu[index1]))
}

t_error <- t(rmvt(n, sigma = sigma, df = 4))
for (i in 1:n){
  f <- t(rmvt(1, sigma = diag(K), df = v))/sqrt_factor
  X[,i] <- mu + B %*% f + t_error[,i]# 
}

results <- ACE(X, H0_indicator = berlii)

X_bar <- rowMeans(X)
X_bar_factor_adjust <- results$Adjusted_mean_difference


df1 <- data.frame(
  sample_mean = as.vector(X_bar),
  group = factor(berlii)
)

df2 <- data.frame(
  sample_mean_adj = as.vector(X_bar_factor_adjust),
  group = factor(berlii)
)

x_limits <- c(min(min(df1$sample_mean), min(df2$sample_mean_adj)),
              max(max(df1$sample_mean), max(df2$sample_mean_adj)))

# 绘制 Sample mean 的直方图
p1 <- ggplot(df1, aes(x = sample_mean, fill = group)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.6, position = "identity") +
  scale_fill_manual(
    name = "Group",
    values = c("0" = "#1868B2", "1" = "#F94141"),
    labels = c("True nulls", "Non-nulls")
  ) +
  geom_vline(xintercept = lb, color = "black", linetype = "dashed", size = 1) +
  labs(title = NULL, x = "Sample mean", y = "Count") +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.8, 0.9),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 14),  # 设置横坐标名称字体大小
        axis.title.y = element_text(size = 14),  # 设置纵坐标名称字体大小
        legend.text = element_text(size = 13)) +
  coord_cartesian(xlim = x_limits, ylim = c(0, 300)) # 图例放在右上角

# 绘制 Sample mean after factor adjustment 的直方图
p2 <- ggplot(df2, aes(x = sample_mean_adj, fill = group)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.6, position = "identity") +
  scale_fill_manual(
    name = "Group",
    values = c("0" = "#1868B2", "1" = "#F94141"),
    labels = c("True nulls", "Non-nulls")
  ) +
  geom_vline(xintercept = lb, color = "black", linetype = "dashed", size = 1) +
  labs(title = NULL, 
       x = "Sample mean after factor adjustment", y = "Count",
       size = 16) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.8, 0.9),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 14),  # 设置横坐标名称字体大小
        axis.title.y = element_text(size = 14),  # 设置纵坐标名称字体大小
        legend.text = element_text(size = 13)) +
  coord_cartesian(xlim = x_limits, ylim = c(0, 300))

# 水平排列两个图
grid.arrange(p1, p2, ncol = 2)
