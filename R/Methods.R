#' Inference method for class "Maximin"
#' @description Point estimator and Confidence interval based on Maximin object
#'
#' @param object Object of class inheriting from "Maximin"
#' @param delta The ridge penalty (Default = 0)
#' @param gen.size The generating sample size (Default = 500)
#' @param threshold Should generated samples be filtered or not?
#' If 0, use normal threshold to filter;
#' if 1, use chi-square threshold to filter;
#' if 2, do not filter (Default = 0)
#' @param alpha confidence value to construct confidence interval (Default = 0.05)
#' @param alpha.thres confidence value to select generated samples (Default = 0.01)
#'
#' @return
#' \item{weight}{The weight vector for groups, of length \eqn{L}}
#' \item{point}{The point estimator of the linear combination}
#' \item{mm.effect}{The aggregated maximin effect (coefficients), of length \eqn{p} or \eqn{p+1}}
#' \item{CI}{Confidence interval for the linear combination}
#'
#' @export
#' @importFrom stats na.omit coef qchisq qnorm
#' @importFrom intervals Intervals interval_union
#' @importFrom MASS mvrnorm
#' @import CVXR
#'
#' @examples
#' \donttest{
#' ## number of groups
#' L=2
#' ## dimension
#' p=500
#'
#' ## mean vector for source
#' mean.source = rep(0, p)
#' ## covariance matrix for source
#' A1gen <- function(rho,p){
#'   A1=matrix(0,p,p)
#'   for(i in 1:p){
#'     for(j in 1:p){
#'       A1[i,j]<-rho^(abs(i-j))
#'     }
#'   }
#'   return(A1)
#' }
#' cov.source = A1gen(0.6, p)
#'
#' ## 1st group's source data
#' n1 = 500
#' X1 = MASS::mvrnorm(n1, mu=mean.source, Sigma=cov.source)
#' b1 = rep(0, p)
#' b1[1:10] = seq(1:10)/40 # true coef for 1st group
#' Y1 = X1%*%b1 + rnorm(n1)
#'
#' ## 2nd group's source data
#' n2 = 400
#' X2 = MASS::mvrnorm(n2, mu=mean.source, Sigma=cov.source)
#' b2 = rep(0, p)
#' b2[1:10] = -seq(1:10)/40 # true coef for 2nd group
#' Y2 = X2%*%b2 + rnorm(n2)
#'
#' ## Target Data, covariate shift
#' n.target = 500
#' mean.target = rep(0, p)
#' cov.target = cov.source
#' for(i in 1:p) cov.target[i, i] = cov.target[i, i] + 0.1
#' for(i in 1:5){
#'   for(j in 1:5){
#'     if(i!=j) cov.target[i, j] = 0.9
#'   }
#' }
#' X.target = MASS::mvrnorm(n.target, mu=mean.target, Sigma=cov.target)
#'
#' ## loading
#' loading = rep(0, p)
#' loading[1:5] = 1
#'
#' ## call
#' mm <- Maximin(list(X1, X2), list(Y1, Y2), loading, X.target, covariate.shift = TRUE)
#' mmInfer <- infer(mm)
#' }
infer <- function(object, delta=0, gen.size=500, threshold=c(0,1,2), alpha=0.05, alpha.thres=0.01){
  ##############################
  ####### Maximin Effects ######
  ##############################
  # aggregated weights
  solution = opt.weight(object$Gamma.prop, delta, report.reward=FALSE)
  weight = solution$weight
  # point estimation of <loading, \beta_\delta^*>
  point = sum(object$Point.vec * weight)
  # maximin effect
  mm.effect = object$Coef.est %*% weight

  #############################
  ##### Generate Samples ######
  #############################
  gen.dim = object$L*(object$L+1)/2
  threshold = match.arg(threshold)
  if(threshold==0){
    thres = qnorm(1-alpha.thres/(gen.dim*2))
    gen.samples = matrix(0, nrow=gen.size, ncol=gen.dim)
    n.picked = 0 # records how many samples are picked
    while(n.picked < gen.size){
      S = mvrnorm(1, mu=rep(0, gen.dim), Sigma=object$gen.Cov)
      if(max(abs(S / sqrt(diag(object$gen.Cov)))) <= thres){
        n.picked = n.picked + 1
        gen.samples[n.picked, ] = object$gen.mu + S
      }
    }
  }
  if(threshold==1){
    # 1 stands for chi square threshold
    UDV_list = svd(object$gen.Cov)
    U = UDV_list$u
    D = UDV_list$d
    V = UDV_list$v
    D.sqrt = sqrt(D)
    gen.Cov.sqrt = U %*% diag(D.sqrt) %*% t(V)
    thres = 1.0 * qchisq((1-alpha.thres), df=gen.dim)
    gen.samples = matrix(0, nrow=gen.size, ncol=gen.dim)
    n.picked = 0
    while(n.picked < gen.size){
      Z = mvrnorm(1, mu=rep(0, gen.dim), Sigma=diag(gen.dim))
      Z.normsq = sum(Z^2)
      if(Z.normsq <= thres){
        n.picked = n.picked + 1
        gen.samples[n.picked, ] = object$gen.mu + gen.Cov.sqrt %*% Z
      }
    }
  }
  if(threshold==2){
    gen.samples = matrix(mvrnorm(gen.size, mu=object$gen.mu, Sigma=object$gen.Cov),
                         nrow=gen.size, ncol=gen.dim)
  }

  ##################################
  ######### Construct CI ###########
  ##################################
  gen.weight.mat = matrix(NA, nrow=gen.size, ncol=object$L)
  gen.est = matrix(NA, nrow=gen.size, ncol=2)
  # construct CI
  for(g in 1:gen.size){
    gen.matrix = matrix(NA, nrow=object$L, ncol=object$L)
    gen.matrix[lower.tri(gen.matrix, diag=TRUE)] = gen.samples[g, ]
    for(l in 1:object$L){
      for(k in 2:object$L){
        gen.matrix[l, k] = gen.matrix[k, l]
      }
    }
    gen.solution = opt.weight(gen.matrix, delta, report.reward=FALSE)
    gen.weight.vector = gen.solution$weight
    gen.weight.mat[g, ] = gen.weight.vector
    gen.point = sum(object$Point.vec * gen.weight.vector)
    gen.se = sqrt(sum(gen.weight.vector^2 * object$SE.vec^2))
    gen.est[g, 1] = gen.point
    gen.est[g, 2] = gen.se
  }

  CI.original = cbind(gen.est[, 1]-qnorm(1-alpha/2)*gen.est[, 2], gen.est[,1]+qnorm(1-alpha/2)*gen.est[,2])
  CI = na.omit(CI.original)
  uni = Intervals(CI)
  CI.union = as.matrix(interval_union(uni))
  colnames(CI.union) <- c("lower", "upper")

  out = list(weight = weight,
             point = point,
             mm.effect = mm.effect,
             CI = CI.union)
  return(out)
}

#' Wrapper function for Maximin inference
#' @description \code{MaximinInfer} is a wrapper for class \code{Maximin} and the method \code{infer}.
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
#' @param delta The ridge penalty (Default = 0)
#' @param gen.size The generating sample size (Default = 500)
#' @param threshold Should generated samples be filtered or not?
#' If 0, use normal threshold to filter;
#' if 1, use chi-square threshold to filter;
#' if 2, do not filter. (Default = 0)
#' @param alpha confidence value to construct confidence interval (Default = 0.05)
#' @param alpha.thres confidence value to select generated samples (Default = 0.01)
#'
#' @return
#' \item{weight}{The weight vector for groups, of length \eqn{L}}
#' \item{point}{The point estimator of the linear combination}
#' \item{mm.effect}{The aggregated maximin effect (coefficients), of length \eqn{p} or \eqn{p+1}}
#' \item{CI}{Confidence interval for the linear combination}
#' @export
#'
#' @examples
#' \donttest{
#' ## number of groups
#' L=2
#' ## dimension
#' p=500
#'
#' ## mean vector for source
#' mean.source = rep(0, p)
#' ## covariance matrix for source
#' A1gen <- function(rho,p){
#'   A1=matrix(0,p,p)
#'   for(i in 1:p){
#'     for(j in 1:p){
#'       A1[i,j]<-rho^(abs(i-j))
#'     }
#'   }
#'   return(A1)
#' }
#' cov.source = A1gen(0.6, p)
#'
#' ## 1st group's source data
#' n1 = 500
#' X1 = MASS::mvrnorm(n1, mu=mean.source, Sigma=cov.source)
#' b1 = rep(0, p)
#' b1[1:10] = seq(1:10)/40 # true coef for 1st group
#' Y1 = X1%*%b1 + rnorm(n1)
#'
#' ## 2nd group's source data
#' n2 = 400
#' X2 = MASS::mvrnorm(n2, mu=mean.source, Sigma=cov.source)
#' b2 = rep(0, p)
#' b2[1:10] = -seq(1:10)/40 # true coef for 2nd group
#' Y2 = X2%*%b2 + rnorm(n2)
#'
#' ## Target Data, covariate shift
#' n.target = 500
#' mean.target = rep(0, p)
#' cov.target = cov.source
#' for(i in 1:p) cov.target[i, i] = cov.target[i, i] + 0.1
#' for(i in 1:5){
#'   for(j in 1:5){
#'     if(i!=j) cov.target[i, j] = 0.9
#'   }
#' }
#' X.target = MASS::mvrnorm(n.target, mu=mean.target, Sigma=cov.target)
#'
#' ## loading
#' loading = rep(0, p)
#' loading[1:5] = 1
#'
#' ## call
#' out <- MaximinInfer(list(X1, X2), list(Y1, Y2), loading, X.target, covariate.shift = TRUE)
#' out$CI
#' }
MaximinInfer <- function(Xlist, Ylist, loading, X.target=NULL, cov.target=NULL,
                         covariate.shift=TRUE, lam.value=c("CV","CV.min"),
                         intercept=TRUE, intercept.loading=FALSE, gen.size=500,
                         delta=0, threshold=c(0,1,2), alpha=0.05, alpha.thres=0.01){
  lam.value = match.arg(lam.value)
  threshold = match.arg(threshold)
  mm <- Maximin(Xlist, Ylist, loading, X.target, cov.target, covariate.shift, lam.value, intercept, intercept.loading)
  out <- infer(mm, gen.size, delta, threshold, alpha, alpha.thres)
  return(out)
}

#' instability measurement
#'
#' @description compute the instability measurement given a specific ridge penalty
#' @param object Object of class inheriting from "Maximin"
#' @param delta The ridge penalty (Default = 0)
#' @param gen.size The generating sample size (Default = 500)
#' @param threshold Should generated samples be filtered or not?
#' If 0, use normal threshold to filter;
#' if 1, use chi-square threshold to filter;
#' if 2, do not filter. (Default = 0)
#' @param alpha.thres confidence value to select generated samples (Default = 0.01)
#'
#' @return
#' \item{measure}{The measurement of instability}
#' @export
measure_instability <- function(object, delta=0, gen.size=500, threshold=c(0,1,2), alpha.thres=0.01){
  spnorm.Gamma.diff = rep(0, gen.size)
  l2norm.gamma.diff = rep(0, gen.size)

  #############################
  ##### Generate Samples ######
  #############################
  gen.dim = object$L*(object$L+1)/2
  threshold = match.arg(threshold)
  if(threshold==0){
    thres = qnorm(1-alpha.thres/(gen.dim*2))
    gen.samples = matrix(0, nrow=gen.size, ncol=gen.dim)
    n.picked = 0 # records how many samples are picked
    while(n.picked < gen.size){
      S = mvrnorm(1, mu=rep(0, gen.dim), Sigma=object$gen.Cov)
      if(max(abs(S / sqrt(diag(object$gen.Cov)))) <= thres){
        n.picked = n.picked + 1
        gen.samples[n.picked, ] = object$gen.mu + S
      }
    }
  }
  if(threshold==1){
    # 1 stands for chi square threshold
    UDV_list = svd(object$gen.Cov)
    U = UDV_list$u
    D = UDV_list$d
    V = UDV_list$v
    D.sqrt = sqrt(D)
    gen.Cov.sqrt = U %*% diag(D.sqrt) %*% t(V)
    thres = 1.0 * qchisq((1-alpha.thres), df=gen.dim)
    gen.samples = matrix(0, nrow=gen.size, ncol=gen.dim)
    n.picked = 0
    while(n.picked < gen.size){
      Z = mvrnorm(1, mu=rep(0, gen.dim), Sigma=diag(gen.dim))
      Z.normsq = sum(Z^2)
      if(Z.normsq <= thres){
        n.picked = n.picked + 1
        gen.samples[n.picked, ] = object$gen.mu + gen.Cov.sqrt %*% Z
      }
    }
  }
  if(threshold==2){
    gen.samples = matrix(mvrnorm(gen.size, mu=object$gen.mu, Sigma=object$gen.Cov),
                         nrow=gen.size, ncol=gen.dim)
  }
  ################################################
  ###### Compute Gamma.diff and gamma.diff #######
  ################################################
  gen.Gamma.array = array(NA, dim=c(gen.size, object$L, object$L))
  gen.weight.mat = matrix(NA, nrow=gen.size, object$L)

  solution = opt.weight(object$Gamma.prop, delta, report.reward=FALSE)
  weight.prop = solution$weight

  for(g in 1:gen.size){
    gen.matrix = matrix(NA, nrow=L, ncol=L)
    gen.matrix[lower.tri(gen.matrix, diag=TRUE)] = gen.samples[g,]
    for(l in 1:L){
      for(k in 2:L){
        gen.matrix[l, k] = gen.matrix[k, l]
      }
    }
    gen.Gamma.array[g,,] = gen.matrix
    gen.solution = opt.weight(gen.matrix, delta, report.reward=FALSE)
    gen.weight.mat[g,] = gen.solution$weight
  }
  for(g in 1:gen.size){
    Gamma.diff = gen.Gamma.array[g,,] - object$Gamma.prop
    spnorm.Gamma.diff[g] = sqrt(max((eigen(Gamma.diff)$values)^2))
    gamma.diff = gen.weight.mat[g,] - weight.prop
    l2norm.gamma.diff[g] = sqrt(sum(gamma.diff^2))
  }
  measure = mean(l2norm.gamma.diff^2 / spnorm.Gamma.diff^2)
  out <- list(measure = measure)
  return(out)
}

#' decide delta data-dependently
#'
#' @description \code{decide_delta} will tell if the estimator is stable or not without ridge penalty
#' at first. If instable, it picks a ridge penalty data-dependently.
#'
#' @param object Object of class inheriting from "Maximin"
#' @param step_delta The step size of searching delta (Default = 0.1)
#' @param MAX_iter Maximum of iterations for searching (Default = 100)
#' @param verbose Print information about delta and reward (Default = `FALSE`)
#'
#' @return
#' \item{delta}{The data-dependent ridge penalty}
#' \item{reward.ratio}{The ratio of penalized reward over non-penalized reward}
#' @export
decide_delta <- function(object, step_delta=0.1, MAX_iter=100, verbose=FALSE){
  #######################################
  ############### STEP-1 ################
  ## Compute Measure square on delta=0 ##
  #######################################
  measure = measure_instability(object, delta=0, gen.size=500, threshold=0, alpha.thres=0.01)$measure
  if(threshold <=0.5){
    print(paste("No ridge penalty suffices to yield a stable estimator"))
    out <- list(delta = 0,
                reward.ratio = 1)
  }
  ####################################################
  ##################### STEP-2 #######################
  ## If rejected above, pick delta data-dependently ##
  ####################################################
  if(threshold > 0.5){

    solution0 = opt.weight(object$Gamma.prop, 0)
    reward0 = solution0$reward

    delta = step_delta
    solution = opt.weight(object$Gamma.prop, delta)
    reward = solution$reward
    i = 1
    while(reward >= 0.95 * reward0){
      delta_new = delta + step_delta
      solution = opt.weight(object$Gamma.prop, delta_new)
      reward = solution$reward
      if(reward <= 0.95 * reward0) break
      if(delta_new > 2){
        warning(paste("The picked delta is over our maximum limit 2. Early Stopping at iteration", i))
        break
      }
      delta = delta_new
      i = i+1
      if((i %% 10 == 0)&verbose){
        print(paste("Iteration ", i, "delta =", round(delta, 4), "reward ratio = ", round(reward/reward0, 4)))
      }
      if(i >= MAX_iter){
        warning("Delta searching stops, because it reach the Max iterations.")
        break
      }
    }
    solution = opt.weight(Gamma, delta)
    reward = solution$reward
    if(verbose){
      print(paste0("The picked delta is ", round(delta,4)))
      print(paste0("Reward Ratio is ", round(reward / reward0, 4)))
    }
    out <- list(delta = delta,
                reward.ratio = reward/reward.ratio)
  }
  return(out)
}
