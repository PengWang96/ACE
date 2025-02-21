e_bh <- function(e_values, alpha) {
  m <- length(e_values)
  if (m == 0) return(integer(0))
  
  # 降序排列e值并记录原始索引
  sorted_indices <- order(e_values, decreasing = TRUE)
  sorted_e <- e_values[sorted_indices]
  
  # 计算每个k对应的阈值
  k <- 1:m
  thresholds <- m / (alpha * k)
  
  # 找出满足条件的k
  satisfied <- which(sorted_e >= thresholds)
  k_hat <- ifelse(length(satisfied) > 0, max(satisfied), 0)
  
  # 返回被拒绝的假设索引
  if (k_hat == 0) {
    integer(0)
  } else {
    sorted_indices[1:k_hat]
  }
}


boostAD <- function(b){
  K <- p; a <- gama; d <- delta
  p<-rep(0,len=K)
  q<-rep(0,len=K)
  p[1]=1-pnorm(log(K/1/a/b)/d + d/2)
  q[1]=K*p[1]
  for (k in 2:K){ 
    p[k]= pnorm(log(K/(k-1)/a/b)/d + d/2) -pnorm(log(K/k/a/b)/d + d/2)
    q[k]=p[k]*K/k
  }
  return(sum(q)-a) 
}