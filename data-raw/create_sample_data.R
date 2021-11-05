############## Toy Example ###############
# The proposed method is designed to high-dimensional settings
# but it can adapt to low-dimensional problems as well. For testing
# purpose, the dimension p=5 in the toy example.
##########################################
set.seed(0)

## number of groups
L=2
## dimension
p=5

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
n1 = 100
X1 = MASS::mvrnorm(n1, mu=mean.source, Sigma=cov.source)
# true coef for 1st group
b1 = rep(0, p)
b1[1] = 0.1
b1[5] = -0.5
Y1 = X1%*%b1 + rnorm(n1)

## 2nd group's source data
n2 = 100
X2 = MASS::mvrnorm(n2, mu=mean.source, Sigma=cov.source)
# true coef for 2nd group
b2 = rep(0, p)
b2[2] = 0.1
b2[5] = 0.5
Y2 = X2%*%b2 + rnorm(n2)

## Target Data, covariate shift
n.target = 100
mean.target = rep(0, p)
cov.target = cov.source
for(i in 1:p) cov.target[i, i] = 1.5
for(i in 1:5){
  for(j in 1:5){
    if(i!=j) cov.target[i, j] = 0.9
  }
}
X.target = MASS::mvrnorm(n.target, mu=mean.target, Sigma=cov.target)

sample_data = list(X1 = X1,
                   Y1 = Y1,
                   X2 = X2,
                   Y2 = Y2,
                   X.target = X.target)

usethis::use_data(sample_data, overwrite=TRUE, compress="xz")
