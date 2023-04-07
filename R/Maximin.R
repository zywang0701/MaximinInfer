#' Returns a list that provides materials for later inference method.
#' @description Given list of observations, compute the bias-corrected initial estimators and do bias-correction to the regressopm covariance matrix.
#' @details
#' The algorithm implemented scenarios with or without covariate shift. If \code{cov0} is specified,
#' the \code{X0} will be ignored; if not, while \code{X0} is specified, \code{cov0} will be estimated
#' by \code{X0}. If both are not specified, the algorithm will automatically set \code{cov.shift} as
#' \code{FALSE}.
#' @param Xlist list of design matrix for source data, of length \eqn{L}
#' @param Ylist list of outcome vector for source data, of length \eqn{L}
#' @param loading.mat Loading matrix, of dimension \eqn{n.loading} x \eqn{p}, each column corresponds to a
#'   loading of interest
#' @param X0 design matrix for target data, of dimension \eqn{n0} x \eqn{p} (default =
#'   \code{NULL})
#' @param cov.shift Covariate shifts or not between source and target data (default = \code{TRUE})
#' @param cov0 Covariance matrix for target data, of dimension \eqn{p} x \eqn{p} (default = \code{NULL})
#' @param intercept Should intercept be fitted for the initial estimator
#'   (default = \code{TRUE})
#' @param intercept.loading Should intercept term be included for the loading
#'   (default = \code{FALSE})
#' @param lambda The tuning parameter in fitting initial model. If \code{NULL},
#'   it will be picked by cross-validation. (default = \code{NULL})
#' @param verbose Should intermediate message(s) be printed. (default = \code{FALSE})
#'
#' @return The returned list contains the following components:
#' \item{Gamma.plugin}{The plugin regression covariance matrix}
#' \item{Gamma.debias}{The proposed debiased regression covariance matrix}
#' \item{Var.Gamma}{The variance matrix for sampling the regression covariance matrix}
#' \item{fits.info}{The list of length \eqn{L}, that contains the initial coefficient estimators and variance of fitted residuals.}
#' \item{Points.info}{The list of length \eqn{L}, that contains the initial debiased estimator for linear combinations and its corresponding standard error.}
#' @export
#'
#' @importFrom stats coef
#' @importFrom SIHR LF
#' @import CVXR glmnet
#'
#' @examples
#' L = 2
#' n1 = n2 = 100; p = 4
#' X1 = MASS::mvrnorm(n1, rep(0,p), Sigma=diag(p))
#' X2 = MASS::mvrnorm(n2, rep(0,p), Sigma=0.5*diag(p))
#' b1 = seq(1,4)/10; b2 = rep(0.2, p)
#' y1 = as.vector(X1%*%b1+rnorm(n1)); y2 = as.vector(X2%*%b2+rnorm(n2))
#' loading1 = rep(0.4, p)
#' loading2 = c(-0.5, -0.5, rep(0,p-2))
#' loading.mat = cbind(loading1, loading2)
#' cov0 = diag(p)
#' mm = Maximin(list(X1,X2),list(y1,y2),loading.mat,cov0=cov0)
#'
#' # inference
#' out = Infer(mm, gen.size=10)
Maximin <- function(Xlist, Ylist, loading.mat, X0=NULL, cov.shift=TRUE, cov0=NULL, intercept=TRUE, intercept.loading=FALSE, lambda=NULL, verbose=FALSE){

  ### Basic Preparation ###
  Xlist = lapply(Xlist, FUN=as.matrix)
  Ylist = lapply(Ylist, FUN=as.vector)
  loading.mat = as.matrix(loading.mat)
  L = length(Xlist)
  p = ncol(Xlist[[1]]) + as.integer(intercept)
  ns = sapply(Xlist, nrow)
  n.loading = ncol(loading.mat)

  ### Check arguments ###
  if(!is.logical(verbose)) verbose=TRUE
  if(intercept==FALSE && intercept.loading==TRUE){
    intercept.loading = FALSE
    cat("Argument 'intercept.loading' is set to FALSE, because 'intercept' is FALSE \n")
  }
  check.args(Xlist, Ylist, loading.mat, X0, cov.shift, cov0, intercept, intercept.loading, lambda, verbose)
  if(!is.null(X0)) X0 = as.matrix(X0)
  if(is.null(X0)&&is.null(cov0)){
    cov.shift=FALSE
    X.source = do.call(rbind, Xlist)
    X0 = X.source
  }

  ### Specify relevant functions ###
  funs.all = relevant.funs(intercept=intercept)
  train.fun = funs.all$train.fun
  pred.fun = funs.all$pred.fun
  dev.fun = funs.all$dev.fun

  ####################################################################
  #################### Part - Initial Estimators #####################
  ####################################################################
  if(verbose) cat("======> Bias Correction for initial estimators.... \n")
  fits.info = rep(list(NA), L)
  Points.info = rep(list(NA), L)
  for(l in 1:L){
    Xlist[[l]] = scale(Xlist[[l]], center=TRUE, scale=FALSE)
    y = Ylist[[l]]
    X = Xlist[[l]]
    beta.init = as.vector(train.fun(X, y, lambda=lambda)$lasso.est)
    sparsity = sum(abs(beta.init)>1e-4)
    pred = pred.fun(X, beta.init)
    dev = dev.fun(pred, y, sparsity)
    Est = LF(X, y, loading.mat, model='linear', intercept=intercept, intercept.loading=intercept.loading,
             beta.init=beta.init, verbose=verbose)
    fits.info[[l]] = list(beta.init = beta.init,
                          dev = dev)
    Points.info[[l]] = list(est.debias.vec = Est$est.debias.vec,
                            se.vec = Est$se.vec)
  }

  ##############################################################################
  #################### Bias-corrected estimators for Gamma #####################
  ##############################################################################
  if(!is.null(X0)){
    ### centralize X0 ###
    X0 = scale(X, center=TRUE, scale=FALSE)
    ### pred0.mat ###
    pred0.mat = matrix(NA, nrow=nrow(X0), ncol=L)
    for(l in 1:L){
      pred0.mat[,l] = pred.fun(X0, fits.info[[l]]$beta.init)
    }
  }

  ### compute Sigma0 ###
  if(is.null(cov0)){
    if(intercept) X0 = cbind(1, X0)
    Sigma0 = t(X0)%*%X0/nrow(X0)
  }
  if(!is.null(cov0)){
    Sigma0 = matrix(0, p, p)
    if(intercept){
      Sigma0[1,1] = 1; Sigma0[-1,-1] = cov0
    }else{
      Sigma0 = cov0
    }
  }

  ### Gamma.plugin ###
  Gamma.plugin = matrix(0, L, L)
  for(l in 1:L) for(k in l:L) Gamma.plugin[l,k] = as.numeric(t(fits.info[[l]]$beta.init)%*%Sigma0%*%(fits.info[[k]]$beta.init))
  for(l in 2:L) for(k in 1:(l-1)) Gamma.plugin[l,k] = Gamma.plugin[k,l]

  ### Bias-corrected estimators: Gamma.debias ###
  if(verbose) cat("======> Bias Correction for matrix Gamma.... \n")
  correct.mat = matrix(0, L, L)
  Proj.array = array(NA, dim=c(L,L,p))
  for(l in 1:L){
    for(k in 1:L){
      loading = as.vector(Sigma0 %*% fits.info[[k]]$beta.init)
      Est.lk = LF(Xlist[[l]], Ylist[[l]], loading, intercept=intercept, beta.init=fits.info[[l]]$beta.init, verbose=verbose)
      correct.mat[l,k] = Est.lk$est.debias.vec - Est.lk$est.plugin.vec
      Proj.array[l,k,] = as.vector(Est.lk$proj.mat)
    }
  }
  Gamma.debias = matrix(0, L, L)
  for(l in 1:L) for(k in l:L) Gamma.debias[l,k] = Gamma.plugin[l,k] + correct.mat[l,k] + correct.mat[k,l]
  for(l in 2:L) for(k in 1:(l-1)) Gamma.debias[l,k] = Gamma.debias[k,l]

  ############################################################
  #################### Variance matrix V #####################
  ############################################################
  # gen.mu = Gamma.debias[lower.tri(Gamma.debias, diag=TRUE)]
  gen.dim = L*(L+1)/2
  Var.Gamma = matrix(NA, nrow=gen.dim, ncol=gen.dim)
  for(k1 in 1:L){
    for(l1 in k1:L){
      for(k2 in 1:L){
        for(l2 in k2:L){
          ind1 = index.map(L,l1,k1)
          ind2 = index.map(L,l2,k2)

          X.l1 = Xlist[[l1]]; X.k1 = Xlist[[k1]]
          if(intercept){
            X.l1 = cbind(1, X.l1); X.k1 = cbind(1, X.k1)
          }
          Sigma.l1 = t(X.l1)%*%X.l1/nrow(X.l1); Sigma.k1 = t(X.k1)%*%X.k1/nrow(X.k1)

          val1 = fits.info[[l1]]$dev/nrow(X.l1)*Proj.array[l1,k1,]%*%Sigma.l1%*%(Proj.array[l2,k2,]*(l2==l1) + Proj.array[k2,l2,]*(k2==l1))
          val2 = fits.info[[k1]]$dev/nrow(X.k1)*Proj.array[k1,l1,]%*%Sigma.k1%*%(Proj.array[l2,k2,]*(l2==k1) + Proj.array[k2,l2,]*(k2==k1))
          if(is.null(cov0)){
            val3 = mean((diag(pred0.mat[,k1]%*%t(pred0.mat[,l1]))-mean(pred0.mat[,k1]*pred0.mat[,l1]))*(diag(pred0.mat[,k2]%*%t(pred0.mat[,l2]))-mean(pred0.mat[,k2]*pred0.mat[,l2])))
            val3 = val3/nrow(X0)
          }else{
            val3 = 0
          }
          val = val1+val2+val3

          Var.Gamma[ind1, ind2] = val
        }
      }
    }
  }
  tau = 0.2
  Var.Gamma = Var.Gamma + diag(max(tau*diag(Var.Gamma), 1/min(ns)), gen.dim)

  obj = list(Gamma.plugin = Gamma.plugin,
             Gamma.debias = Gamma.debias,
             Var.Gamma = Var.Gamma,
             fits.info = fits.info,
             Points.info = Points.info)
  obj
}

#' Inference method
#' @description Given the returned list of Maximin, compute the Point estimator and Confidence interval.
#'
#' @param obj returned list of Maximin
#' @param delta The ridge penalty (Default = 0)
#' @param gen.size The generating sample size (Default = 500)
#' @param threshold Should generated samples be filtered or not?
#' if 0, use normal threshold to filter;
#' if 1, use chi-square threshold to filter;
#' if 2, do not filter (Default = 0)
#' @param alpha confidence value to construct confidence interval (Default = 0.05)
#' @param alpha.thres confidence value to select generated samples (Default = 0.01)
#'
#' @return
#' \item{weight}{The weight vector for groups, of length \eqn{L}}
#' \item{mm.effect}{The aggregated maximin effect (coefficients), of length \eqn{p} or \eqn{p+1}}
#' \item{mminfer}{The list of length \eqn{n.loading}, each contains the point estimator and confidence interval}
#' @export
#'
#' @importFrom stats na.omit qchisq qnorm
#' @importFrom intervals Intervals interval_union
#' @importFrom MASS mvrnorm
#' @import CVXR
Infer <- function(obj, delta=0, gen.size=500, threshold=0, alpha=0.05, alpha.thres=0.01){
  n.loading = length(obj$Points.info[[1]]$est.debias.vec)
  L = nrow(obj$Gamma.debias)
  ## points.mat records each loading each group's debiased point estimator
  points.mat = do.call(cbind, lapply(obj$Points.info, FUN=function(x) x$est.debias.vec)) # dim of (n.loading, L)
  ## se.mat records each loading each group's standard error
  se.mat = do.call(cbind, lapply(obj$Points.info, FUN=function(x) x$se.vec)) # dim of (n.loading, L)

  ####################### Bias-corrected Point Estimation #######################
  ### solve weights using Gamma.debias ###
  sol = opt.weight(obj$Gamma.debias, delta, report.reward=FALSE)
  weight = sol$weight

  ### Point Estimation of <loading, \beta_\delta^*>
  Point.debias = as.vector(points.mat%*%weight)

  ### Maximin Effect ###
  Coefs = do.call(cbind, lapply(obj$fits.info, FUN=function(x) x$beta.init)) # (p, L)
  mm.effect = as.vector(Coefs %*% weight)

  ####################### Sampling #######################
  gen.mu = obj$Gamma.debias[lower.tri(obj$Gamma.debias, diag=TRUE)]
  gen.Cov = obj$Var.Gamma
  gen.samples = gensamples(gen.mu, gen.Cov, gen.size, threshold, alpha.thres)
  ### weights for each generated sample ###
  gen.weight.mat = matrix(NA, nrow=gen.size, ncol=L)
  for(g in 1:gen.size){
    gen.matrix = matrix(NA, nrow=L, ncol=L)
    gen.matrix[lower.tri(gen.matrix, diag=TRUE)] = gen.samples[g, ]
    for(l in 1:L){
      for(k in l:L){
        gen.matrix[l, k] = gen.matrix[k, l]
      }
    }
    gen.sol = opt.weight(gen.matrix, delta, report.reward=FALSE)
    gen.weight.mat[g, ] = gen.sol$weight
  }

  #################### Construct CI ######################
  CIs = rep(list(NA), n.loading)
  for(i.loading in 1:n.loading){
    points = points.mat[i.loading,]
    ses = se.mat[i.loading,]
    point.gen.vec = as.vector(gen.weight.mat %*% points)
    se.gen.vec = sqrt(as.vector(gen.weight.mat^2 %*% ses^2))
    CI.ori = cbind(point.gen.vec - qnorm(1-alpha/2)*se.gen.vec,
                   point.gen.vec + qnorm(1-alpha/2)*se.gen.vec)
    CI = na.omit(CI.ori)
    uni = Intervals(CI)
    CI.union = as.matrix(interval_union(uni))
    colnames(CI.union) <- c('lower', 'upper')
    CIs[[i.loading]] = CI.union
  }

  #################### Output #####################
  mminfer = rep(list(NA), n.loading)
  for(i.loading in 1:n.loading){
    mminfer[[i.loading]] = list(point = Point.debias[i.loading],
                                CI = CIs[[i.loading]])
  }
  out = list(weight = weight,
             mm.effect = mm.effect,
             mminfer = mminfer)
  return(out)
}

#' measurement of instability
#' @description compute the instability measurement given a specific ridge penalty
#'
#' @param obj The returned list of Maximin
#' @param delta The ridge penalty (Default = 0)
#' @param gen.size The generating sample size (Default = 500)
#' @param threshold Should generated samples be filtered or not?
#' if 0, use normal threshold to filter;
#' if 1, use chi-square threshold to filter;
#' if 2, do not filter. (Default = 0)
#' @param alpha.thres The confidence value to select generated samples (Default = 0.01)
#'
#' @return The measurement of instability
#' @export
measure_instability <- function(obj, delta=0, gen.size=500, threshold=0, alpha.thres=0.01){
  L = nrow(obj$Gamma.debias)
  gen.dim = L*(L+1)/2
  gen.mu = obj$Gamma.debias[lower.tri(obj$Gamma.debias, diag=TRUE)]
  gen.Cov = obj$Var.Gamma

  spnorm.Gamma.diff = rep(0, gen.size)
  l2norm.gamma.diff = rep(0, gen.size)

  ################################################
  ###### Compute Gamma.diff and gamma.diff #######
  ################################################
  gen.samples = gensamples(gen.mu, gen.Cov, gen.size, threshold, alpha.thres)

  gen.Gamma.array = array(NA, dim=c(gen.size, L, L))
  gen.weight.mat = matrix(NA, nrow=gen.size, L)

  sol = opt.weight(obj$Gamma.debias, delta, report.reward=FALSE)
  weight.prop = sol$weight

  for(g in 1:gen.size){
    gen.matrix = matrix(NA, nrow=L, ncol=L)
    gen.matrix[lower.tri(gen.matrix, diag=TRUE)] = gen.samples[g,]
    for(l in 1:L){
      for(k in l:L){
        gen.matrix[l, k] = gen.matrix[k, l]
      }
    }
    gen.Gamma.array[g,,] = gen.matrix
    gen.sol = opt.weight(gen.matrix, delta, report.reward=FALSE)
    gen.weight.mat[g,] = gen.sol$weight
  }
  for(g in 1:gen.size){
    Gamma.diff = gen.Gamma.array[g,,] - obj$Gamma.debias
    spnorm.Gamma.diff[g] = sqrt(max((eigen(Gamma.diff)$values)^2))
    gamma.diff = gen.weight.mat[g,] - weight.prop
    l2norm.gamma.diff[g] = sqrt(sum(gamma.diff^2))
  }
  measure = mean(l2norm.gamma.diff^2 / spnorm.Gamma.diff^2)
  return(measure)
}


#' Decide ridge penalty data-dependently
#' @description  To tell if the estimator is stable or not without ridge penalty
#' at first. If instable, it picks a ridge penalty data-dependently.
#' @param obj The returned list of Maximin
#' @param gen.size The generating sample size (Default = 500)
#' @param step_delta The step size of searching delta (Default = 0.1)
#' @param MAX_iter Maximum of iterations for searching (Default = 100)
#' @param verbose Print information about delta and reward (Default = \code{FALSE})
#'
#' @return
#' \item{delta}{The data-dependent ridge penalty}
#' \item{reward.ratio}{The ratio of penalized reward over non-penalized reward}
#' @export
decide_delta <- function(obj, gen.size=500, step_delta=0.1, MAX_iter=100, verbose=FALSE){
  #######################################
  ############### STEP-1 ################
  ## Compute Measure square on delta=0 ##
  #######################################
  measure = measure_instability(obj, delta=0, gen.size=gen.size, threshold=0, alpha.thres=0.01)
  if(measure <=0.5){
    print(paste("ridge penalty 0 suffices to yield a stable estimator"))
    out <- list(delta = 0,
                reward.ratio = 1)
  }
  ####################################################
  ##################### STEP-2 #######################
  ## If rejected above, pick delta data-dependently ##
  ####################################################
  if(measure > 0.5){

    sol0 = opt.weight(obj$Gamma.debias, 0)
    reward0 = sol0$reward

    delta = 0
    reward = reward0
    for(i in 1:MAX_iter){
      delta_new = delta + step_delta
      sol = opt.weight(obj$Gamma.debias, delta_new)
      reward = sol$reward
      if(reward <= 0.95 * reward0) break
      delta = delta_new
      if(delta > 2){
        cat("The picked delta reaches the maximum limit 2.\n")
        break
      }
    }

    sol = opt.weight(obj$Gamma.debias, delta)
    reward = sol$reward
    if(verbose){
      print(paste0("The picked delta is ", round(delta,4)))
      print(paste0("Reward Ratio is ", round(reward / reward0, 4)))
    }
    out <- list(delta = delta,
                reward.ratio = reward/reward0)
  }

  return(out)
}
