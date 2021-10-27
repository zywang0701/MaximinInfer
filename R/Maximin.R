#' Class Maximin
#'
#' @description
#' `Maximin` returns the class "Maximin", which provides materials for later inference method.
#'
#' @details
#' The algorithm implemented scenarios with or without covariate shift. If `cov.target` is specified,
#' the `X.target` will be ignored; if not, while `X.target` is specified, `cov.target` will be estimated
#' by `X.target`. If both are not specified, the algorithm will automatically set `covariate.shift` as
#' `FALSE`.
#'
#' @param Xlist list of design matrix for source data, of length \eqn{L}
#' @param Ylist list of outcome vector for source data, of length \eqn{L}
#' @param loading Loading, of length \eqn{p}
#' @param X.target Design matrix for target data, of dimension \eqn{n.target} x \eqn{p} (default = `NULL`)
#' @param cov.target Covariance matrix for target data, of dimension \eqn{p} x \eqn{p} (default = `NULL`)
#' @param covariate.shift Covariate shifts or not between source and target data (default = `TRUE`)
#' @param lam.value The method to be used to obtain Lasso estimator of high-dimensional regression vector for each group
#' @param intercept Should intercept be fitted for the initial estimator (default = `TRUE`)
#' @param intercept.loading Should intercept be included for the loading (default = `FALSE`)
#'
#' @return
#' `Maximin` returns an object of class "Maximin". The function `infer` is used to do further inference.
#' An object of class "Maximin" is a list containing the following components.
#' \item{Gamma.prop}{The proposed debiased regression covariance matrix}
#' \item{Coef.est}{matrix, of dimension \eqn{p(+1)} x \eqn{L} where each column corresponds to the Lasso estimator of the high-dimensional regression vector for a given group}
#' \item{Point.vec}{vector, of length \eqn{L} with the l-th entry as the debiased estimator of the linear combination of the l-th high-dimensional regression vector}
#' \item{L}{The number of groups}
#' \item{gen.mu}{The mean vector for sampling the regression covariance matrix}
#' \item{gen.Cov}{The variance matrix for sampling the regression covariance matrix}
#'
#' @export
#' @importFrom stats median
#' @importFrom SIHR LF
#' @import CVXR glmnet
Maximin <- function(Xlist, Ylist, loading, X.target=NULL, cov.target=NULL,
                    covariate.shift=TRUE, lam.value=c("CV","CV.min"),
                    intercept=TRUE, intercept.loading=FALSE){

  if(!intercept) intercept.loading=FALSE
  lam.value = match.arg(lam.value)

  ############################################
  ########### Transfer Source Data ###########
  ############################################
  if((is.list(Xlist)==FALSE)||(is.list(Ylist)==FALSE)) stop("Error: check the type of Xlist and Ylist, they must be list")
  if(length(Xlist)!=length(Ylist)) stop("Error: check the length of Xlist and Ylist")
  L = length(Xlist)
  n.x.vec = rep(0, L)
  n.y.vec = rep(0, L)
  p.vec = rep(0, L)
  for(l in 1:L){
    n.x.vec[l] = dim(Xlist[[l]])[1]
    p.vec[l] = dim(Xlist[[l]])[2]
    n.y.vec[l] = length(Ylist[[l]])
    if(n.x.vec[l]!=n.y.vec[l]) stop(paste("Error: X and Y to the group",l,"has different number of samples"))
  }
  if(!all(p.vec==p.vec[1])) stop("Error: check the dimension p of each X, they must be the same")
  X.source = do.call(rbind, Xlist)
  Y.source = do.call(c, Ylist)
  idx.source = rep(1:L, times=n.y.vec)

  ##################################
  ########### Check Input ##########
  ##################################
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



