
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MaximinInfer

<!-- badges: start -->
<!-- badges: end -->

The goal of MaximinInfer is to provide functionality for the paper. The
function is used to compute the bias corrected estimator of
ridge-penalized maximin effect and the point estimator of its linear
contrast. It also constructs the confidence interval for the linear
contrast.

## Installation

You can install the released version of MaximinInfer from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("MaximinInfer")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("zywang0701/MaximinInfer")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MaximinInfer)
# The number of groups
L = 2
# dimension
p = 500
# sample sizes for each group's source data
ns.source = c(500, 400)
# sample size for target data
n.target = 2000
# covariance matrix
A1gen <- function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    }
  }
  return(A1)
}
cov.source = A1gen(0.6,p)
cov.target = cov.source
for(i in 1:p) cov.target[i, i] = 1.5
for(i in 1:5){
  for(j in 1:5){
    if(i!=j) cov.target[i, j] = 0.9
  }
}
for(i in 499:500){
  for(j in 499:500){
    if(i!=j) cov.target[i, j] = 0.9
  }
}
# generate mean for source and target data
mean.source = rep(0, p)
mean.target = rep(0, p)
# coefs
Bs = matrix(0, p, L)
Bs[1:10, 1] = seq(1:10)/40
Bs[22:23,1] = 1
Bs[498,1] = 0.5
Bs[499,1] = -0.5
Bs[500,1] = -0.5
Bs[1:10, 2] = Bs[1:10, 1]
Bs[22:23, 2] = 1
Bs[500, 2] = 1
# loading
loading = rep(0, p)
loading[498:500] = 1

X.source = MASS::mvrnorm(sum(ns.source), mu=mean.source, Sigma=cov.source)
X.target = MASS::mvrnorm(n.target, mu=mean.target, Sigma=cov.target)
idx.source = rep(1:L, times=ns.source)
Y.source = rep(0, sum(ns.source))
for(l in 1:L){
  idx.l = which(idx.source==l)
  Y.source[idx.l] = X.source[idx.l, ] %*% Bs[,l] + rnorm(ns.source[l])
}
mmList <- mmInfer(X.source, Y.source, idx.source, loading, X.target, cov.target=NULL,
                  covariate.shift=TRUE, split=FALSE, delta=-1, gen.size=10)
#> [1] "delta path starts from 2.0267 which exceeds our maximum limit 2"
#> [1] "The picked delta is 2"
#> [1] "Reward Ratio is 0.9921"
#> [1] "Minimum Eigenvalue plus delta = 4.4991"
```

Weights for groups

``` r
mmList$weights
#> [1] 0.4896628 0.5103372
```

Point estimator for the linear contrast

``` r
mmList$point.est
#> [1] 0.3498103
```

Confidence Interval for point estimator

``` r
mmList$CI
#>           [,1]      [,2]
#> [1,] 0.1891654 0.4989892
```
