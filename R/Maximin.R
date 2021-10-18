#' Class Maximin
#'
#' @param X.source Design matrix for source data, of dimension \eqn{n.source} x \eqn{p}
#' @param Y.source Outcome vector for source data, of length \eqn{n.source}
#' @param idx.source Indicator vector of groups for source data, of length \eqn{n.source}
#' @param loading Loading, of length \eqn{p}
#' @param X.target Design matrix for target data, of dimension \eqn{n.target} x \eqn{p} (default = \code{NULL})
#' @param cov.target Covariance matrix for target data, of dimension \eqn{p} x \eqn{p}. If set as \code{NULL}, `cov.target` is unknown. (default = \code{NULL})
#' @param covariate.shift Covariate shifts or not between source and target data (default = \code{TRUE})
#' @param lam.value The method to be used to obtain each group's intial estimator
#' @param intercept Should intercept be fitted for the initial estimator (default = \code{TRUE})
#' @param intercept.loading Should intercept be included for the \code{loading} (default = \code{FALSE})
#'
#' @return
#' \code{Maximin} returns an object of class "Maximin". The function \code{infer} is used to do further inference.\\
#' An object of class "Maximin" is a list containing the following components.
#' \item{Gamma.prop}{The proposed debiased Weight matrix}
#' \item{Coef.est}{The initial estimators for each group}
#' \item{Point.vec}{The point estimator for each group}
#' \item{L}{The number of groups}
#' \item{gen.mu}{The mean vector for sampling method}
#' \item{gen.Cov}{The variance matrix for sampling method}
#'
#' @export
#' @importFrom stats median
#' @importFrom scalreg scalreg
#' @importFrom flare slim
#' @importFrom SIHR LF
#' @import CVXR glmnet
#'
#' @examples
#' \donttest{
#' ## number of groups
#' L=2
#' ## dimension
#' p=500
#' ## sample size for each group of source data
#' ns.source = c(500, 400)
#' ## sample size for target data
#' n.target=1000
#'
#' A1gen <- function(rho,p){
#'   A1=matrix(0,p,p)
#'   for(i in 1:p){
#'     for(j in 1:p){
#'       A1[i,j]<-rho^(abs(i-j))
#'     }
#'   }
#'   return(A1)
#' }
#' ## mean vector
#' mean.source = rep(0, p)
#' mean.target = rep(0, p)
#' ## covariate shifts
#' cov.source = A1gen(0.6, p)
#' cov.target = diag(p)
#' ## true coefficients
#' Bs = matrix(0, p, L)
#' Bs[1:10,1] = seq(1:10)/40
#' Bs[1:10,2] = -seq(1:10)/40
#' ## Data
#' X.source = MASS::mvrnorm(sum(ns.source), mu=mean.source, Sigma=cov.source)
#' X.target = MASS::mvrnorm(n.target, mu=mean.target, Sigma=cov.target)
#' idx.source = rep(1:L, times=ns.source)
#' Y.source = rep(0, sum(ns.source))
#' for(l in 1:L){
#'   idx.l = which(idx.source==l)
#'   Y.source[idx.l] = X.source[idx.l, ] %*% Bs[,l] + rnorm(ns.source[l])
#' }
#' loading = rep(0, p)
#' loading[1:5] = 1
#' mm <- Maximin(X.source, Y.source, idx.source, loading, X.target, covariate.shift = TRUE)
#' mmInfer <- infer(mm, gen.size=100, delta=0)
#' }
Maximin <- function(X.source, Y.source, idx.source, loading, X.target=NULL, cov.target=NULL,
                    covariate.shift=TRUE, lam.value=c("CV","CV.min","scalreg","slim"),
                    intercept=TRUE, intercept.loading=FALSE){
  if(!intercept) intercept.loading=FALSE
  lam.value = match.arg(lam.value)

  ##################################
  ########### Check Input ##########
  ##################################
  if((dim(X.source)[1]!=length(Y.source))||(dim(X.source)[1]!=length(idx.source))){
    stop("Error: check dimensions of Source Data")
  }
  if(dim(X.source)[2]!=length(loading)){
    stop("Error: check length of loading")
  }
  if((!is.null(X.target))&&(dim(X.source)[2]!=dim(X.target)[2])){
    stop("Error: check dimensions of Target Data")
  }
  if((!is.null(cov.target))&&((dim(X.source)[2]!=dim(cov.target)[2])||(dim(cov.target)[1]!=dim(cov.target)[2]))){
    stop("Error: check dimensions of cov.target")
  }

  if(is.null(X.target)) X.target = X.source  # the code adapts to No target setting
  n.source = nrow(X.source)
  n.target = nrow(X.target)
  p = ncol(X.source) + as.integer(intercept)
  L = length(unique(idx.source))
  uni_groups = sort(unique(idx.source))

  Coef.est = matrix(0, p, L)  # estimators of groups
  Pred.vec = rep(0, n.source)  # predicted outcome of source data
  Pred.mat.target = matrix(0, n.target, L)  # predicted outcomes of target data
  Point.vec = rep(0, L)  # debiased point estimators
  Var.vec = rep(0, L)  # variance of residual
  SE.vec = rep(0, L)  # SE of loading
  for(l in 1:L){
    ## obtain estimators of group l using Lasso
    index.set = which(idx.source==uni_groups[l])
    X = X.source[index.set, ]
    Y = Y.source[index.set]
    col.norm = 1/sqrt(1/nrow(X)*diag(t(X)%*%X))
    X.norm = X %*% diag(col.norm)
    Coef.est[, l] = Lasso(X.norm, Y, lambda=lam.value, intercept=intercept)
    if(intercept){
      Coef.est[-1, l] = Coef.est[-1, l]*col.norm
      Pred.vec[index.set] = X%*%Coef.est[-1, l] + Coef.est[1, l]
      Pred.mat.target[, l] = X.target%*%Coef.est[-1, l] + Coef.est[1, l]
    }else{
      Coef.est[, l] = Coef.est*col.norm
      Pred.vec[index.set] = X%*%Coef.est[, l]
      Pred.mat.target[, l] = X.target%*%Coef.est[, l]
    }
    ## obtain variance of residual for group l
    supp.l = which(abs(Coef.est[, l])>0.01)
    n.eff = max(0.9*nrow(X), nrow(X)-length(supp.l))
    Var.vec[l] = sum((Y - Pred.vec[index.set])^2) / n.eff
    ## bias correction for <loading, b^{(l)}>, leveraged by package SIHR
    est <- LF(X, Y, loading, intercept.loading=intercept.loading, intercept=intercept, init.coef=Coef.est[,l], verbose=FALSE)
    SE.vec[l] = est$se
    Point.vec[l] = est$prop.est
  }
  #########################################
  ######### compte Gamma.plugin ###########
  #########################################
  if(!is.null(cov.target)){
    Sigma.target.est = matrix(0, p, p)
    if(intercept){
      Sigma.target.est[1, 1] = 1
      Sigma.target.est[-1, -1] = cov.target
    }else{
      Sigma.target.est = cov.target
    }
  }else{
    if(covariate.shift){
      if(intercept){
        X.target.b = cbind(1, X.target)
        Sigma.target.est = (1/n.target)*(t(X.target.b)%*%X.target.b)
      }else{
        Sigma.target.est = (1/n.target)*(t(X.target)%*%X.target)
      }
    }else{
      if(intercept){
        X.source.b = cbind(1, X.source)
        X.target.b = cbind(1, X.target)
        Sigma.target.est = (t(X.source.b)%*%X.source.b + t(X.target.b)%*%X.target.b)/(n.source+n.target)
      }else{
        Sigma.target.est = (t(X.source)%*%X.source + t(X.target)%*%X.target)/(n.source + n.target)
      }
    }
  }
  Gamma.plugin = t(Coef.est)%*%Sigma.target.est%*%Coef.est
  Omega.est = Sigma.target.est%*%Coef.est

  ####################################################
  ##### conduct bias correction for Gamma.plugin #####
  ####################################################
  Gamma.prop = Gamma.plugin
  Proj.array = array(NA, dim=c(L, L, p))
  for(l in 1:L){
    for(k in l:L){
      index.set.l = which(idx.source==uni_groups[l])
      index.set.k = which(idx.source==uni_groups[k])
      X.l = X.source[index.set.l, ]
      X.k = X.source[index.set.k, ]
      Y.l = Y.source[index.set.l]
      Y.k = Y.source[index.set.k]
      Pred.l = Pred.vec[index.set.l]
      Pred.k = Pred.vec[index.set.k]

      if(intercept){
        X.l = cbind(1, X.l)
        X.k = cbind(1, X.k)
      }
      if(covariate.shift){
        output <- Gamma.shift(Gamma.plugin[l, k], X.l, X.k, Omega.est[, l], Omega.est[, k],
                              Y.l, Y.k, Pred.l, Pred.k)
        Gamma.prop[l, k] = output$est
        Proj.array[l, k, ] = output$proj.lk
        Proj.array[k, l, ] = output$proj.kl
      }else{
        Gamma.prop[l, k] = Gamma.plugin[l, k]+t(Coef.est[,l])%*%t(X.k)%*%(Y.k-Pred.k)/nrow(X.k)+
          t(Coef.est[,k])%*%t(X.l)%*%(Y.l-Pred.l)/nrow(X.l)
        Proj.array[l, k, ] = Coef.est[,k]
        Proj.array[k, l, ] = Coef.est[,l]
      }
    }
  }
  for(l in 2:L){
    for(k in 1:(l-1)){
      Gamma.prop[l, k] = Gamma.prop[k, l]
    }
  }
  ######################################################################
  ################## to obtain sampling materials ######################
  ## compute mean and covariance matrix for the sampling distribution ##
  ######################################################################
  gen.mu = Gamma.prop[lower.tri(Gamma.prop, diag=TRUE)]
  gen.dim = L*(L+1)/2
  gen.Cov = matrix(NA, nrow=gen.dim, ncol=gen.dim)
  for(k1 in 1:L){
    for(l1 in k1:L){
      index1 = index.map(L, l1, k1)
      for(k2 in 1:L){
        for(l2 in k2:L){
          index2 = index.map(L, l2, k2)
          index.set.l1 = which(idx.source==uni_groups[l1])
          index.set.k1 = which(idx.source==uni_groups[k1])
          index.set.l2 = which(idx.source==uni_groups[l2])
          index.set.k2 = which(idx.source==uni_groups[k2])
          X.l1 = X.source[index.set.l1, ]
          X.k1 = X.source[index.set.k1, ]
          X.l2 = X.source[index.set.l2, ]
          X.k2 = X.source[index.set.k2, ]
          if(intercept){
            X.l1 = cbind(1, X.l1)
            X.k1 = cbind(1, X.k1)
            X.l2 = cbind(1, X.l2)
            X.k2 = cbind(1, X.k2)
          }
          if(!is.null(cov.target)){
            gen.Cov[index1, index2] <- cov.inner.shift.known(Var.vec, l1, k1, l2, k2,
                                                             X.l1, X.k1, X.l2, X.k2,
                                                             Proj.array)
          }else{
            gen.Cov[index1, index2] <- cov.inner.shift(Var.vec, l1, k1, l2, k2,
                                                       X.l1, X.k1, X.l2, X.k2,
                                                       Pred.mat.target, Proj.array)
          }
        }
      }
    }
  }
  tau = 0.2
  gen.Cov = gen.Cov + diag(max(tau*diag(gen.Cov), 1/floor(n.source/L)), dim(gen.Cov)[2])

  out = list(Gamma.prop = Gamma.prop,
             Coef.est = Coef.est,
             Point.vec = Point.vec,
             SE.vec = SE.vec,
             L = L,
             gen.mu = gen.mu,
             gen.Cov = gen.Cov)
  structure(out, class = "Maximin")
}



