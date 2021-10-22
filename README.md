
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MaximinInfer

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/MaximinInfer)](https://CRAN.R-project.org/package=MaximinInfer)
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
```

``` r
set.seed(0)
## dimension
p=500

###### Source Data ######
## number of groups
L=2
## mean vector for source
mean.source = rep(0, p)
## covariance matrix for source
A1gen <- function(rho,p){
A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    }
  }
  return(A1)
}
cov.source = A1gen(0.6, p)

## 1st group's source data
n1 = 500
X1 = MASS::mvrnorm(n1, mu=mean.source, Sigma=cov.source)
b1 = rep(0, p)
b1[1:10] = seq(1:10)/40 # true coef for 1st group
Y1 = X1%*%b1 + rnorm(n1)

## 2nd group's source data
n2 = 400
X2 = MASS::mvrnorm(n2, mu=mean.source, Sigma=cov.source)
b2 = rep(0, p)
b2[1:10] = -seq(1:10)/40 # true coef for 2nd group
Y2 = X2%*%b2 + rnorm(n2)

###### Target Data ######
## covariate shift
n.target = 500
mean.target = rep(0, p)
cov.target = cov.source
for(i in 1:p) cov.target[i, i] = cov.target[i, i] + 0.1
for(i in 1:5){
  for(j in 1:5){
    if(i!=j) cov.target[i, j] = 0.9
  }
}
X.target = MASS::mvrnorm(n.target, mu=mean.target, Sigma=cov.target)

## loading
loading = rep(0, p)
loading[1:5] = 1

## call
mm <- Maximin(list(X1, X2), list(Y1, Y2), X.target, loading, covariate.shift = TRUE)
mmInfer <- infer(mm)
#> Warning in decide_delta(object$Gamma.prop, step_delta = 0.1): Fail to find a
#> suitable delta, the estimator may be not stable enough.
```

Data-dependent ridge penalty

``` r
mmInfer$delta
#> [1] 0.001883843
```

Weights for groups

``` r
mmInfer$weight
#> [1] 0.5136854 0.4863146
```

Point estimator for the linear contrast

``` r
mmInfer$point
#> [1] -0.05622201
```

Confidence Interval for point estimator

``` r
mmInfer$CI
#>           lower      upper
#> [1,] -0.2104918 0.09401328
```
