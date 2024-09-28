# ACE: Analysis of Correlated High-Dimensional Expression Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

ACE is an R package for the analysis of correlated high-dimensional expression data. It provides an implementation of the ACE method, which estimates factor models and produces factor-adjusted statistics.

## Installation

You can install the package directly from GitHub with the following R code:

```R
# install.packages("devtools")
devtools::install_github("PengWang96/ACE")
```

## Loading the ACE Package

After installing the ACE package, you can load it into your R environment using the `library()` function.

```R
# Load the ACE package
library(ACE)
```

Once the package is loaded, you can use the functions provided by the ACE package in your R scripts. Remember, if you encounter any errors during the loading process, it might be because the package is not installed. In that case, refer to the installation instructions above.

If you want to check whether the ACE package is successfully loaded, you can use the `search()` function, which shows all available packages in your current R environment. Here's how:

```R
# Check if the ACE package is loaded
search()
```

In the output, you should see "package:ACE", which indicates that the ACE package is successfully loaded and ready to use.

## Usage

The main function in this package is `ACE()`. This function takes in the following parameters:

- `Z`: The observed data matrix with the variables in rows and samples in columns. It is a p-by-n1 matrix.
- `X`: (Optional) The observed data matrix with the variables in rows and samples in columns. It is a p-by-n2 matrix. If X is present, then perform the two-sample test; otherwise, perform one-sample test.
- `H0_indicator`: (Optional) A p-dimensional vector containing only 0 and 1. A value of 1 means the variable/gene is non-null and a value of 0 means the gene is null.
- `gama`: FDR control level. The default is 0.05.
- `h_max`: The upper bound of the number of latent factors specified by the user. The default is 20.
- `h_fix`: The fixed number of latent factors specified by user. If `NULL`, the number of latent factors is estimated from the data.
- `reg`: The type of regularization to use. It can be either `"L1"` for L1 regularization using quantile regression, or `"L2"` for L2 regularization using least squares. The default is `"L1"`.


The function returns an object with S3 class `ACE` containing several important statistics and estimations.

Here is an example of how to use the `ACE()` function:

```R
library(mvtnorm); library(quantreg); library(ACE)
p <- 200; n <- 100; h <- 3 # the number of variables, samples and factors
berlii <- rbinom(p, 1, 0.2) # 1 means the variable is non-null and 0 means it is null.
index0 <- which(berlii == 0); index1 <- which(berlii == 1)

mu <- matrix(rep(0, 1*p), nrow=p)
mu[index1] <- runif(length(index1), min=0.4, max=0.7) # expectation of data
B <- matrix(runif(h*p, min=-1, max=1), nrow=p) # factor loading matrix
t_error <- t(rmvt(n, sigma = diag(p), df = 10)) # error term followed t-distribution
f <- t(rmvt(n, diag(h), df = 4))/sqrt(4/(4-2)) # factor followed t-distribution
Y <- mu %*% matrix(rep(1, n*1), nrow=1) + B %*% f + t_error # data
res <- ACE(Z = Y, H0_indicator = berlii, gama = 0.05)
res$FDP # true FDP
res$Power # power
```

## References

- Cao, H., & Kosorok, M. R. (2011). Simultaneous critical values for t-tests in very high dimensions. Bernoulli, 17, 347.
- Wang, P., Lyu, P., Peddada, S., Cao, H. (2024+). A powerful methodology for analyzing correlated high dimensional data using factor models. Results not shown.

## License

This project is licensed under the MIT License.
