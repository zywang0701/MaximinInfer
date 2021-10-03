#' Inference for Ridge-Penalized Maximin Effect
#'
#' @description
#' This function is used to compute the bias corrected estimator of ridge-penalized maximin effect and
#' the point estimator of its linear contrast. It also constructs the confidence interval for the linear contrast.
#'
#' @param X.source Design matrix of source data, consisting of multi-groups.
#' @param Y.source Outcome vector of source data, consisting of multi-groups.
#' @param idx.source Index vector of source data, each element indicates the group from which the sample is drawn
#' @param loading Loading, the linear contrast of interest
#' @param X.target Design matrix of target data, the function investigates source data only if \code{NULL} (default = \code{NULL})
#' @param cov.target Covariance matrix for target data, the covariance is unkown if \code{NULL} (default = \code{NULL})
#' @param covariate.shift Covariate is shifted or not between target and source data (default = \code{TRUE})
#' @param split Samples splitted or not to fit algorithm (default = \code{FALSE})
#' @param lam.value The tuning parameter in the construction of LASSO estimator (default = "CV.min")
#' @param intercept Should intercept(s) be fitted or not (default = \code{TRUE})
#' @param delta The ridge penalty level, -1 indicates a data-dependent way to decide the penalty (default = -1)
#' @param gen.size The number of generated samples for inference (default = 500)
#' @param threshold version of selecting generated samples
#' @param alpha confidence value
#'
#' @return
#' \item{mm.est}{The estimator of maximin effects}
#' \item{weights}{The weights of groups, sum up to 1}
#' \item{point.est}{The point estimator of the linear contrast}
#' \item{CI}{Confidence Interval for the linear contrast}
#' \item{CI.length}{Length of the confidence interval}
#'
#' @importFrom stats na.omit
#' @importFrom intervals Intervals interval_union
#' @import CVXR glmnet
#' @export
#'
#' @examples
#' L=2
#' p=500
#' n=500
#' n.target=2000
#' A1gen <- function(rho,p){
#'   A1=matrix(0,p,p)
#'   for(i in 1:p){
#'     for(j in 1:p){
#'       A1[i,j]<-rho^(abs(i-j))
#'     }
#'   }
#'   return(A1)
#' }
#' Bs = matrix(0, p, L)
#' Bs[1:10,1] = seq(1:10)/40
#' Bs[1:10,2] = -seq(1:10)/40
#' cov.source = A1gen(0.6, p)
#' cov.target = cov.source
#' diag(cov.target) = diag(cov.target) + 0.2
#' for(i in 1:5){
#'   for(j in 1:5){
#'     if(i!=j) cov.target[i, j] = 0.9
#'   }
#' }
#' X.source = MASS::mvrnorm(n*L, mu=rep(0,p), Sigma=cov.source)
#' X.target = MASS::mvrnorm(n.target, mu=rep(0,p), Sigma=cov.target)
#' idx.source = rep(1:L, times=rep(n,L))
#' Y.source = rep(0, n*L)
#' for(l in 1:L) Y.source[seq((l-1)*n+1, l*n)] = X.source[seq((l-1)*n+1, l*n), ] %*% Bs[,l] + rnorm(n)
#' loading = rep(0, p)
#' loading[1:5] = 1
#' mmList <- mmInfer(X.source, Y.source, idx.source, loading, X.target, cov.target=NULL,
#'                   covariate.shift=TRUE, split=FALSE, delta=-1)
mmInfer <- function(X.source, Y.source, idx.source, loading, X.target=NULL,
                    cov.target=NULL, covariate.shift=TRUE, split=FALSE,
                    lam.value="CV.min", intercept=TRUE, delta=-1,
                    gen.size=500, threshold=1, alpha=0.01){
  if(split){
    if((covariate.shift==FALSE)||(is.null(cov.target)==FALSE)){
      split = FALSE
      print("In the No Covariate Shift Setting or when Cov.target is known, the sample splitting has no effect.")
    }
  }

  if(split) s1 <- mm.s1.split(X.source, Y.source, idx.source, X.target, loading, NULL,
                              TRUE, TRUE, lam.value, intercept)
  if(!split) s1 <- mm.s1.nosplit(X.source, Y.source, idx.source, X.target, loading, cov.target,
                                 covariate.shift, FALSE, lam.value, intercept)

  s2 <- mm.s2(s1$Gamma.prop, s1$Coef.est, s1$Point.vec, delta)
  s3 <- mm.s3(s1$gen.mu, s1$gen.Cov, s1$gen.dim, gen.size, threshold, alpha)
  s4 <- mm.s4(s1$Point.vec, s1$SE.vec, s1$L, s3$gen.samples, s3$gen.size, s2$delta)

  returnList <- list("mm.est"=s2$mm.est,
                     "weights"=s2$weight.prop,
                     "point.est"=s2$point.est,
                     "CI"=s4$CI.union,
                     "CI.length"=s4$CI.length
                     )
  return(returnList)
}

mm.s1.nosplit <- function(X.source, Y.source, idx.source, X.target=NULL, loading=NULL, cov.target=NULL,
                                   covariate.shift=TRUE, lam.value="CV.min", intercept=TRUE){
  ####################### 1st Part: to obtain Gamma.prop #######################
  if(is.null(X.target)) X.target = X.source  # the code can adapt to No target setting
  n.source = nrow(X.source)
  n.target = nrow(X.target)
  p = ncol(X.source)
  L = length(unique(idx.source))
  uni_groups = sort(unique(idx.source))

  ## check dimensions
  if(ncol(X.source)!=ncol(X.target)) stop("X.source and X.target have different dimensions")
  if(n.source!=length(Y.source)) stop("X.source and Y.source have different number of samples")
  if(n.source!=length(idx.source)) stop("length of idx.source differs from X.source")

  if(intercept) p = p+1

  Coef.est = matrix(0, p, L)  # estimators of groups
  Pred.vec = rep(0, n.source)  # predicted outcome via corresponding group's estimator
  Pred.mat.target = matrix(0, n.target, L)  # predicted target outcome via L groups' estimators
  Point.vec = rep(0, L)  # hat of <loading, b^{(l)}>
  Var.vec = rep(0, L)  # hat of sigma^2_l
  SE.vec = rep(0, L)  # later used in SE(loading)
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
    ## bias correction for <loading, b^{(l)}>
    est <- LF(X, Y, loading, intercept=intercept, init.Lasso=Coef.est[, l])
    SE.vec[l] = est$se
    Point.vec[l] = est$prop.est
  }
  ## compte Gamma.plugin
  if(!is.null(cov.target)){
    Sigma.target.est = matrix(0, p, p)
    if(intercept){
      if(dim(cov.target)[1] != (p-1)) stop("cov.target has wrong dimensions")
      Sigma.target.est[1, 1] = 1
      Sigma.target.est[-1, -1] = cov.target
    }else{
      if(dim(cov.target)[1] != p) stop("cov.target has wrong dimensions")
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

  ## conduct bias correction for Gamma.plugin
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

  ################## 2nd Part: to obtain sampling materials ####################
  ## compute mean and covariance matrix for the sampling distribution
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
  ## conduct sample procedure
  tau = 0.2
  gen.Cov = gen.Cov + diag(max(tau*diag(gen.Cov), 1/floor(n.source/L)), dim(gen.Cov)[2])

  returnList = list("Gamma.prop"=Gamma.prop,
                    "Coef.est"=Coef.est,
                    "Point.vec"=Point.vec,
                    "SE.vec"=SE.vec,
                    "L"=L,
                    "gen.mu"=gen.mu,
                    "gen.Cov"=gen.Cov,
                    "gen.dim"=gen.dim)
  return(returnList)
}

mm.s1.split <- function(X.source, Y.source, idx.source, X.target=NULL, loading=NULL, cov.target=NULL,
                        covariate.shift=TRUE, split=TRUE, lam.value="CV.min", intercept=TRUE){

  # In the split version, it must be covariate shift setting and cov.target is unknown.
  ####################### 1st Part: to obtain Gamma.prop #######################
  if(is.null(X.target)) X.target = X.source
  n.source = nrow(X.source)
  n.target = nrow(X.target)
  p = ncol(X.source)
  L = length(unique(idx.source))
  uni_groups = sort(unique(idx.source))

  ## check dimensions
  if(ncol(X.source)!=ncol(X.target)) stop("X.source and X.target have different dimensions")
  if(n.source!=length(Y.source)) stop("X.source and Y.source have different number of samples")
  if(n.source!=length(idx.source)) stop("length of idx.source differs from X.source")
  ## intercept or not
  if(intercept) p = p+1

  ## sample splitting
  # part: source data
  index.A.List = list()
  index.B.List = list()
  for(l in 1:L){
    index.set = which(idx.source==uni_groups[l])
    n.l = length(index.set)
    if(n.l < 2) stop("Some group has samples fewer than 2, which cannot be splitted")
    n.l.A = floor(n.l/2)
    index.A = sample(index.set, size=n.l.A, replace=FALSE)
    index.B = setdiff(index.set, index.A)
    index.A.List[[l]] = index.A
    index.B.List[[l]] = index.B
  }
  index.A = do.call(c, index.A.List)
  index.B = do.call(c, index.B.List)
  n.source.A = length(index.A)
  n.source.B = length(index.B)
  X.source.A = X.source[index.A,]
  X.source.B = X.source[index.B,]
  Y.source.A = Y.source[index.A]
  Y.source.B = Y.source[index.B]
  idx.source.A = idx.source[index.A]
  idx.source.B = idx.source[index.B]
  # part: target data
  n.target.A = floor(n.target/2)
  n.target.B = n.target - n.target.A
  index.A = sample(seq(1,n.target), size = n.target.A, replace = FALSE)
  index.B = setdiff(seq(1,n.target), index.A)
  X.target.A = X.target[index.A,]
  X.target.B = X.target[index.B,]

  ### MAIN PROGRAM ###
  Coef.est = matrix(0, p, L)  # estimators of groups
  Pred.vec = rep(0, n.source.B)  # predicted outcome via corresponding group's estimator
  Pred.mat.target = matrix(0, n.target, L)  # predicted target outcome via L groups' estimators
  Point.vec = rep(0, L)  # hat of <loading, b^{(l)}>
  Var.vec = rep(0, L)  # hat of sigma^2_l
  SE.vec = rep(0, L)  # later used in SE(loading)

  for(l in 1:L){
    index.set = which(idx.source==uni_groups[l])
    index.set.A = which(idx.source.A==uni_groups[l])
    index.set.B = which(idx.source.B==uni_groups[l])
    X = X.source[index.set,]
    Y = Y.source[index.set]
    X.A = X.source.A[index.set.A, ]
    Y.A = Y.source.A[index.set.A]
    X.B = X.source.B[index.set.B, ]
    Y.B = Y.source.B[index.set.B]
    ## obtain estimators of group l by applying Lasso to the sub-sample source.A
    col.norm = 1/sqrt(1/nrow(X.A)*diag(t(X.A)%*%X.A))
    X.A.norm = X.A %*% diag(col.norm)
    Coef.est[, l] = Lasso(X.A.norm, Y.A, lambda=lam.value, intercept=intercept)
    if(intercept){
      Coef.est[-1, l] = Coef.est[-1, l]*col.norm
      # only Source.B is relevant in eq(19)
      Pred.vec[index.set] = X.B%*%Coef.est[-1, l] + Coef.est[1, l]
      # use all target in eq(21)
      Pred.mat.target[, l] = X.target%*%Coef.est[-1, l] + Coef.est[1, l]
    }else{
      Coef.est[, l] = Coef.est*col.norm
      Pred.vec[index.set] = X.B%*%Coef.est[, l]
      Pred.mat.target[, l] = X.target%*%Coef.est[, l]
    }
    ## obtain variance of residual for group l
    supp.l = which(abs(Coef.est[, l])>0.01)
    n.eff = max(0.9*nrow(X), nrow(X)-length(supp.l))
    # use all source in eq(21)
    if(intercept){
      Var.vec[l] = sum((Y - (X%*%Coef.est[-1,l]+Coef.est[1,l]))^2) / n.eff
    }else{
      Var.vec[l] = sum((Y - (X%*%Coef.est[,l]))^2) / n.eff
    }
    ## bias correction for <loading, b^{(l)}>
    # use all source in eq(25)
    est <- LF(X, Y, loading, intercept=intercept, init.Lasso=Coef.est[, l])
    SE.vec[l] = est$se
    Point.vec[l] = est$prop.est
  }
  ## compute Gamma.plugin
  # covariate shift, shift unknown
  if(intercept){
    X.target.B.b = cbind(1, X.target.B)
    Sigma.target.hat = t(X.target.B.b)%*%X.target.B.b/nrow(X.target.B.b)
    X.target.A.b = cbind(1, X.target.A)
    Sigma.target.tilde = t(X.target.A.b)%*%X.target.A.b/nrow(X.target.A.b)
  }else{
    Sigma.target.hat = t(X.target.B)%*%X.target.B/nrow(X.target.B)
    Sigma.target.tilde = t(X.target.A)%*%X.target.A/nrow(X.target.A)
  }
  Gamma.plugin = t(Coef.est)%*%Sigma.target.hat%*%Coef.est
  Omega.est = Sigma.target.tilde%*%Coef.est

  ## conduct bias correction for Gamma.plugin, where Source.B is used only.
  Gamma.prop = Gamma.plugin
  Proj.array = array(NA, dim=c(L, L, p))
  for(l in 1:L){
    for(k in l:L){
      index.set.l = which(idx.source.B==uni_groups[l])
      index.set.k = which(idx.source.B==uni_groups[k])
      X.l = X.source.B[index.set.l, ]
      X.k = X.source.B[index.set.k, ]
      Y.l = Y.source.B[index.set.l]
      Y.k = Y.source.B[index.set.k]
      Pred.l = Pred.vec[index.set.l]
      Pred.k = Pred.vec[index.set.k]

      if(intercept){
        X.l = cbind(1, X.l)
        X.k = cbind(1, X.k)
      }
      # covariate.shift
      output <- Gamma.shift(Gamma.plugin[l, k], X.l, X.k, Omega.est[, l], Omega.est[, k],
                            Y.l, Y.k, Pred.l, Pred.k)
      Gamma.prop[l, k] = output$est
      Proj.array[l, k, ] = output$proj.lk
      Proj.array[k, l, ] = output$proj.kl
    }
  }
  for(l in 2:L){
    for(k in 1:(l-1)){
      Gamma.prop[l, k] = Gamma.prop[k, l]
    }
  }

  ################## 2nd Part: to obtain sampling materials ####################
  ## compute mean and covariance matrix for the sampling distribution
  gen.mu = Gamma.prop[lower.tri(Gamma.prop, diag=TRUE)]
  gen.dim = L*(L+1)/2
  gen.Cov = matrix(NA, nrow=gen.dim, ncol=gen.dim)
  for(k1 in 1:L){
    for(l1 in k1:L){
      index1 = index.map(L, l1, k1)
      for(k2 in 1:L){
        for(l2 in k2:L){
          index2 = index.map(L, l2, k2)
          index.set.l1 = which(idx.source.B==uni_groups[l1])
          index.set.k1 = which(idx.source.B==uni_groups[k1])
          index.set.l2 = which(idx.source.B==uni_groups[l2])
          index.set.k2 = which(idx.source.B==uni_groups[k2])
          X.l1 = X.source.B[index.set.l1, ]
          X.k1 = X.source.B[index.set.k1, ]
          X.l2 = X.source.B[index.set.l2, ]
          X.k2 = X.source.B[index.set.k2, ]
          if(intercept){
            X.l1 = cbind(1, X.l1)
            X.k1 = cbind(1, X.k1)
            X.l2 = cbind(1, X.l2)
            X.k2 = cbind(1, X.k2)
          }

          gen.Cov[index1, index2] <- cov.inner.shift(Var.vec, l1, k1, l2, k2,
                                                     X.l1, X.k1, X.l2, X.k2,
                                                     Pred.mat.target, Proj.array)

        }
      }
    }
  }
  ## conduct sample procedure
  tau = 0.2
  gen.Cov = gen.Cov + diag(max(tau*diag(gen.Cov), 1/floor(n.source/L)), dim(gen.Cov)[2])


  returnList = list("Gamma.prop"=Gamma.prop,
                    "Coef.est"=Coef.est,
                    "Point.vec"=Point.vec,
                    "SE.vec"=SE.vec,
                    "L"=L,
                    "gen.mu"=gen.mu,
                    "gen.Cov"=gen.Cov,
                    "gen.dim"=gen.dim)
  return(returnList)
}

mm.s2 <- function(Gamma.prop, Coef.est, Point.vec, delta=-1){
  if(delta==-1) delta = decide_delta(Gamma.prop, step_delta=0.1)
  ####################### Maximin Effects ##############################
  # aggregated weights
  solution = opt.weight(Gamma.prop, delta, report.reward=FALSE)
  weight.vector = solution$weight
  # point estimation of <loading, \beta_\delta^*>
  point.est = sum(Point.vec * weight.vector)
  # maximin effect
  mm.est = Coef.est %*% weight.vector
  returnList = list("mm.est"=mm.est,
                    "weight.prop"=weight.vector,
                    "point"=point.est,
                    "delta"=delta)
  return(returnList)
}

mm.s3 <- function(gen.mu, gen.Cov, gen.dim, gen.size=500, threshold=0, alpha=0.01){
  ## Varying Approaches to generate samples
  if(threshold==0){
    # 0 stands for the original approach
    gen.samples = matrix(mvrnorm(gen.size, mu=gen.mu, Sigma=gen.Cov), nrow=gen.size, ncol=gen.dim)
  }
  if(threshold==1){
    # 1 stands for chi square threshold
    UDV_list = svd(gen.Cov)
    U = UDV_list$u
    D = UDV_list$d
    V = UDV_list$v
    D.sqrt = sqrt(D)
    gen.Cov.sqrt = U %*% diag(D.sqrt) %*% t(V)
    thres = 1.0 * qchisq((1-alpha), df=gen.dim) #1.05
    gen.samples = matrix(0, nrow=gen.size, ncol=gen.dim)
    n.picked = 0
    for(i.gen in 1:gen.size){
      Z = mvrnorm(1, mu=rep(0, gen.dim), Sigma=diag(gen.dim))
      Z.normsq = sum(Z^2)
      if(Z.normsq <= thres){
        n.picked = n.picked + 1
        gen.samples[n.picked, ] = gen.mu + gen.Cov.sqrt %*% Z
      }
    }
    gen.samples = gen.samples[1:n.picked, ]
    gen.size = n.picked
  }
  if(threshold == 2){
    thres = qnorm(1-alpha/(gen.dim*2))
    gen.samples = matrix(0, nrow=gen.size, ncol=gen.dim)
    n.picked = 0 # records how many samples are picked

    for(i.gen in 1:gen.size){
      S = mvrnorm(1, mu=rep(0, gen.dim), Sigma=gen.Cov)
      if(max(abs(S / sqrt(diag(gen.Cov)))) <= thres){
        n.picked = n.picked + 1
        gen.samples[n.picked, ] = gen.mu + S
      }
    }
    gen.samples = gen.samples[1:n.picked, ]
    gen.size = n.picked
  }

  returnList <- list("gen.samples"=gen.samples,
                     "gen.size"=gen.size)
  return(returnList)
}

mm.s4 <- function(Point.vec, SE.vec, L, gen.samples, gen.size, delta){
  ####################### Sampling #####################################
  gen.weight.mat = matrix(NA, nrow=gen.size, ncol=L)
  gen.est = matrix(NA, nrow=gen.size, ncol=2)
  # construct CI
  for(g in 1:gen.size){
    gen.matrix = matrix(NA, nrow=L, ncol=L)
    gen.matrix[lower.tri(gen.matrix, diag=TRUE)] = gen.samples[g, ]
    for(l in 1:L){
      for(k in 2:L){
        gen.matrix[l, k] = gen.matrix[k, l]
      }
    }
    gen.solution = opt.weight(gen.matrix, delta, report.reward=FALSE)
    gen.weight.vector = gen.solution$weight
    gen.weight.mat[g, ] = gen.weight.vector
    gen.point = sum(Point.vec * gen.weight.vector)
    gen.se = sqrt(sum(gen.weight.vector^2 * SE.vec^2))
    gen.est[g, 1] = gen.point
    gen.est[g, 2] = gen.se
  }
  CI.original = cbind(gen.est[, 1]-1.96*gen.est[, 2], gen.est[,1]+1.96*gen.est[,2])
  CI = na.omit(CI.original)
  uni = Intervals(CI)
  CI.union = as.matrix(interval_union(uni))
  CI.length = sum(CI.union[,2] - CI.union[,1])

  returnList = list("gen.est"=gen.est,
                    "CI.union"=CI.union,
                    "CI.length"=CI.length)
  return(returnList)
}
