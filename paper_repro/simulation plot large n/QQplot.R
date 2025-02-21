rm(list = ls())
library(ggplot2)
library(mvtnorm)
library(gridExtra)
library(ACE)

p <- 2000; n1 <- 600; n2 <- 800; K1 <- 1; K2 <- 2; lb <- 0.1; ub <- 0.4; v <- 4
PP_values <- c(0.2, 0.05)
factor_indicator_values <- c(1, 0)

plot_qq <- function(PP, factor_indicator) {
  berlii <- rbinom(p, 1, PP)
  index0 <- which(berlii == 0)
  index1 <- which(berlii == 1)
  true_pai1 <- length(index1) / p

  mu <- matrix(rep(0, 1*p), nrow=p)
  mu[index0] <- runif(length(index0), 1, 4)
  mu[index1] <- runif(length(index1), 1, 4)
  mu2 <- mu
  mu[index1] <- runif(length(index1), min=lb, max=ub) + mu[index1]

  X <- matrix(rep(0, n1*p), nrow=p)
  B <- factor_indicator * matrix(runif(K1*p, min=-1, max=1), nrow=p)
  Y <- matrix(rep(0, n2*p), nrow=p)
  B2 <- factor_indicator * matrix(runif(K2*p, min=-1, max=1), nrow=p)

  sigma1 <- diag(c(2 + mu - mu2))
  sigma2 <- diag(c(2.5 + mu - mu2))
  diag(sigma1[index0,index0]) <- (runif(length(index0), min=1, max=1))
  diag(sigma2[index0,index0]) <- (runif(length(index0), min=1, max=1))
  t_error <- t(rmvt(n1, sigma = sigma1, df = 4))
  t_error2 <- t(rmvt(n2, sigma = sigma2, df = 10))

  for (i in 1:n1){
    f <- t(rmvt(1, sigma = diag(K1), df = v))/sqrt(v / (v - 2))
    X[,i] <- mu + B %*% f + t_error[,i]
  }

  for (i in 1:n2){
    f2 <- t(rmvt(1, sigma = diag(K2), df = v))/sqrt(v / (v - 2))
    Y[,i] <- mu2 + B2 %*% f2 + t_error2[,i]
  }

  ace_result <- ACE(X, Y, berlii)
  statistics <- ace_result$Adjusted_statistics


  qq_data <- data.frame(statistics = statistics)

  if (factor_indicator == 1) {
    title_expr <- bquote(pi[1] == .(PP) ~ ", B" ~ " ≠ " ~ 0)
  } else {
    title_expr <- bquote(pi[1] == .(PP) ~ ", B" ~ " = " ~ 0)
  }

  p_gg <- ggplot(qq_data, aes(sample = statistics)) +
    stat_qq() +
    stat_qq_line(col = "red") +
    theme_bw() +
    labs(title = title_expr) +
    theme(text = element_text(size = 10),  # 设置字体大小
          plot.title = element_text(size = 15, hjust = 0.5),  # 加大标题字体
          axis.title.x = element_blank(),  # 去掉 x 轴标题
          axis.title.y = element_blank(),  # 去掉 y 轴标题
          axis.text.x = element_text(size = 10),  # 设置坐标轴刻度字体大小
          axis.text.y = element_text(size = 10))   # 设置坐标轴刻度字体大小

  return(p_gg)
}


plot_list <- list()
for (PP in PP_values) {
  for (factor_indicator in factor_indicator_values) {
    plot_list[[length(plot_list) + 1]] <- plot_qq(PP, factor_indicator)
  }
}
grid.arrange(grobs = plot_list, ncol = 2, nrow = 2)
