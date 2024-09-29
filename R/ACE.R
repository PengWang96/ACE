#' @title Analysis of Correlated High-Dimensional Expression (ACE) Data
#' @description A function for estimating factor models, giving factor-adjusted statistics.
#'
#' @param Z The observed data matrix with the variables in rows and samples in columns. It is a \eqn{p}-by-\eqn{n_1} matrix.
#' @param X (Optional) The observed data matrix with the variables in rows and samples in columns. It is a \eqn{p}-by-\eqn{n_2} matrix.
#' If X is present, then perform the two-sample test; otherwise, perform one-sample test.
#' @param H0_indicator (Optional) A \eqn{p}-dimensional vector containing only 0 and 1.
#' A value of 1 means the variable/gene is non-null and a value of 0 means the gene is null.
#' @param gama FDR control level. The default is 0.05.
#' @param h_max The upper bound of the number of latent factors specified by the user. The default is 20.
#' @param h_fix (Optional) The fixed number of latent factors specified by user.
#' If `NULL`, the number of latent factors is estimated from the data. The default is `NULL`.
#' @param reg The type of regularization to use. It can be either `"L1"` for L1 regularization using quantile regression,
#' or `"L2"` for L2 regularization using least squares. The default is `"L1"`.
#'
#' @return An object with S3 class \code{ACE} containing the following items will be returned:
#' \describe{
#' \item{\code{FDP}}{If H0_indicator exists, FDP is true FDP, otherwise, it is estimated FDP.}
#' \item{\code{Power}}{If H0_indicator exists, power is output which is defined as the ratio of the number of correctly rejected to the number of non-nulls.}
#' \item{\code{Rejection}}{The number of rejections.}
#' \item{\code{Adjusted_mean_difference}}{Factor-adjusted mean difference which is a \eqn{p}-dimensional vector.}
#' \item{\code{Adjusted_statistics}}{Factor-adjusted statistics (\eqn{p}-dimensional vector).}
#' \item{\code{Threshold}}{A critical value. When absolute factor-adjusted statistics is larger than the threshold, we reject it.}
#' \item{\code{Estimated_number_factor}}{The estimated number of factors.}
#' \item{\code{pai1_hat}}{The estimated proportion of non-nulls.}
#' \item{\code{Estimate_factor_loadings}}{The estimated factor loadings.}
#' \item{\code{Estimate_factors}}{The estimated factors.}
#' }
#'
#' @importFrom stats cov pnorm quantile
#' @importFrom quantreg rq
#' @importFrom irlba irlba
#'
#' @references Cao, H., & Kosorok, M. R. (2011). Simultaneous critical values for t-tests in very high dimensions. Bernoulli, 17, 347.
#' @references Wang, P., Lyu, P., Peddada, S., Cao, H. (2024+). A powerful methodology for analyzing correlated high dimensional data using factor models. results not shown.
#'
#' @examples
#' library(mvtnorm); library(quantreg)
#' p <- 200; n <- 100; h <- 3 # the number of variables, samples and factors
#' berlii <- rbinom(p, 1, 0.2) # 1 means the variable is non-null and 0 means it is null.
#' index0 <- which(berlii == 0); index1 <- which(berlii == 1)
#'
#' mu <- matrix(rep(0, 1*p), nrow=p)
#' mu[index1] <- runif(length(index1), min=0.4, max=0.7) # expectation of data
#' B <- matrix(runif(h*p, min=-1, max=1), nrow=p) # factor loading matrix
#' t_error <- t(rmvt(n, sigma = diag(p), df = 10)) # error term followed t-distribution
#' f <- t(rmvt(n, diag(h), df = 4))/sqrt(4/(4-2)) # factor followed t-distribution
#' Y <- mu %*% matrix(rep(1, n*1), nrow=1) + B %*% f + t_error # data
#' res <- ACE(Z = Y, H0_indicator = berlii, gama = 0.05)
#' res$FDP # true FDP
#' res$Power # power
#'
#' @export
#'

ACE <- function(Z, X, H0_indicator, gama = 0.05, h_max = 20, h_fix = NULL, reg = "L1"){
  if (missing(X)) {
    Y <- Z
    n <- ncol(Y)
    p <- nrow(Y)
  } else {
    p <- nrow(Z); n1 <- ncol(Z); n2 <- ncol(X)
    if (n1 <= n2) {
      n <- n1
      Y <- matrix(0, p, n)
    } else {
      Y <- Z
      Z <- X
      X <- Y
      n1 <- ncol(Z); n2 <- ncol(X); n <- n1
    }

    for (jy in 1:n1){
      Y[,jy] <- Z[,jy] - sqrt(n1/n2)*X[,jy] + apply(X[,1:n1], 1, sum)/sqrt(n1*n2) - apply(X,1,mean)
    }
  }

  Yba <- apply(Y,1,mean); T <- (sqrt(n)* Yba); deltaa <- cov(t(Y))

  # pca <- svd(deltaa, nu=0, nv = h_max)
  pca <- irlba(deltaa, nv = h_max, symmetric = TRUE)
  lam_sort <- pca$d
  gamma_norm <- pca$v

  if (is.null(h_fix)) {
    bizhi <- lam_sort[1:(h_max-1)] / lam_sort[2:h_max]
    h_hat <- which.max(bizhi)
  } else {
    h_hat <- h_fix
  }

  if (h_hat == 1){
    B_hat <- gamma_norm[,1:h_hat] * sqrt(lam_sort[1:h_hat])
  } else {
    B_hat <- gamma_norm[,1:h_hat] %*% diag(sqrt(lam_sort[1:h_hat]))
  }


  if (reg == "L1") {
    W0_hat <- matrix(0, h_hat, n)
    for (jjj in 1:n) {
      W0_hat[, jjj] <- rq(Y[, jjj] ~ B_hat - 1, tau = 0.5)$coef
    }
  } else if (reg == "L2") {
    W0_hat <- solve(t(B_hat) %*% B_hat) %*% t(B_hat) %*% Y
  }

  W_piao <- rbind(rep(1,n), W0_hat)

  Px <- t(W_piao) %*% solve(W_piao %*% t(W_piao)) %*% W_piao
  muB_hat <- Y %*% t(W_piao) %*% solve(W_piao %*% t(W_piao))
  mu_hat <- muB_hat[,1]; B_hat <- muB_hat[,-1]

  T_k <- sqrt(n)*mu_hat
  sigma_hat <- (Y %*% (diag(n) - Px) %*% t(Y))/(n - h_hat - 1)
  bbb <- sqrt(diag(sigma_hat))
  statistics <- T_k/bbb
  abs_stat <- abs(statistics)

  if (n < 100) {
    aaaaaaa <- seq(0.3, 5, by = 0.001)
    pai <- function(c){
      (sapply(c,
              function(c, abs_stat){
                mean(pmin(abs_stat, c))
              }
              , abs_stat) / c - 2*(1-exp(-c^2/2))/(c*sqrt(2*pi)) -
         2 * pnorm(c, lower.tail = F))/(1 -2*(1-exp(-c^2/2))/(c*sqrt(2*pi)) -
                                          2 * pnorm(c, lower.tail = F))
    }
    paii <- pai(aaaaaaa)
    pai1 <- max(paii[paii >= 0 & paii <= 1])

    s0 <- 0.1*mean(quantile(bbb, seq(0, 1-pai1, by = 0.001)))
    sigma_hat <- diag((bbb + s0)^2)
    bbb <- sqrt(diag(sigma_hat))
    statistics <- T_k/bbb
    abs_stat <- abs(statistics)
  } else {pai1 <- 0}



  aaaaaa <- seq(0.01, 10, by = 0.005)
  f_t_hat <- function(x){
    2*pnorm(x, lower.tail = F) -
      gama * sapply(x, function(x, abs_stat){mean(( abs_stat >= x ))}, abs_stat)
  }
  f_t_hatt <- f_t_hat(aaaaaa)
  index <- which(f_t_hatt <= 0)[1]
  t_fdr_hat <- ifelse(length(index) <= 0, Inf, aaaaaa[index])

  R <- sum(abs_stat >= t_fdr_hat)

  if (missing(H0_indicator)) {
    FDP <- 2*p*pnorm(t_fdr_hat, lower.tail = F) / R
    return(list("FDP" = FDP, "Rejection" = R, "Adjusted_mean_difference" = mu_hat,
                "Adjusted_statistics" = statistics, "Threshold" = t_fdr_hat,
                "Estimated_number_factor" = h_hat, "pai1_hat" = pai1,
                "Estimate_factor_loadings" = B_hat, "Estimate_factors" = W0_hat))
  } else {
    index_0 <- which(H0_indicator == 0)
    false_reject <- which(abs(T_k[index_0]/bbb[index_0]) >= t_fdr_hat)
    if (R == 0){
      S <- 0; true_FDP <- 0
    } else {
      S <- R - length(false_reject)
      true_FDP <- length(false_reject)/R
    }
    power <- S/(p - length(index_0))
    return(list("FDP" = true_FDP, "Power" = power, "Rejection" = R,
                "Adjusted_mean_difference" = mu_hat, "Adjusted_statistics" = statistics,
                "Threshold" = t_fdr_hat, "Estimated_number_factor" = h_hat, "pai1_hat" = pai1,
                "Estimate_factor_loadings" = B_hat, "Estimate_factors" = W0_hat))
  }
}
