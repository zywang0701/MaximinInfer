#' Inference method for class "Maximin"
#' @description Point estimator and Confidence interval based on Maximin object
#'
#' @param object Object of class inheriting from "Maximin"
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
#' ## The problem is low-dimensional and we do sampling only 5 times instead of 500 for testings
#' ## heterogenous data and covariates shift
#' X1 = sample_data$X1
#' X2 = sample_data$X2
#' Y1 = sample_data$Y1
#' Y2 = sample_data$Y2
#' X.target = sample_data$X.target
#'
#' ## loading
#' loading = rep(0, 5) # dimension p=5
#' loading[5] = 1
#'
#' ## call
#' mm <- Maximin(list(X1, X2), list(Y1, Y2), loading, X.target, covariate.shift = TRUE)
#' mmInfer <- infer(mm, gen.size=5)
infer <- function(object, delta=0, gen.size=500, threshold=0, alpha=0.05, alpha.thres=0.01){
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

#' The Wrapper function for Maximin inference
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
#' ## The problem is low-dimensional and we do sampling only 5 times instead of 500 for testings
#' ## heterogenous data and covariates shift
#' X1 = sample_data$X1
#' X2 = sample_data$X2
#' Y1 = sample_data$Y1
#' Y2 = sample_data$Y2
#' X.target = sample_data$X.target
#'
#' ## loading
#' loading = rep(0, 5) # dimension p=5
#' loading[5] = 1
#'
#' ## call
#' mmInfer <- MaximinInfer(list(X1, X2), list(Y1, Y2), loading, X.target, gen.size=5)
MaximinInfer <- function(Xlist, Ylist, loading, X.target=NULL, cov.target=NULL,
                         covariate.shift=TRUE, lam.value=c("CV","CV.min"),
                         intercept=TRUE, intercept.loading=FALSE, delta=0,
                         gen.size=500, threshold=0, alpha=0.05, alpha.thres=0.01){
  lam.value = match.arg(lam.value)
  mm <- Maximin(Xlist, Ylist, loading, X.target, cov.target, covariate.shift, lam.value, intercept, intercept.loading)
  out <- infer(mm, delta, gen.size, threshold, alpha, alpha.thres)
  return(out)
}

#' measurement of instability
#' @description compute the instability measurement given a specific ridge penalty
#'
#' @param object Object of class inheriting from "Maximin"
#' @param delta The ridge penalty (Default = 0)
#' @param gen.size The generating sample size (Default = 500)
#' @param threshold Should generated samples be filtered or not?
#' if 0, use normal threshold to filter;
#' if 1, use chi-square threshold to filter;
#' if 2, do not filter. (Default = 0)
#' @param alpha.thres confidence value to select generated samples (Default = 0.01)
#'
#' @return
#' \item{measure}{The measurement of instability}
#' @export
#' @examples
#' ## The problem is low-dimensional and we do sampling only 5 times instead of 500 for testings
#' ## heterogenous data and covariates shift
#' X1 = sample_data$X1
#' X2 = sample_data$X2
#' Y1 = sample_data$Y1
#' Y2 = sample_data$Y2
#' X.target = sample_data$X.target
#'
#' ## loading
#' loading = rep(0, 5)
#' loading[5] = 1
#'
#' ## call
#' mm <- Maximin(list(X1, X2), list(Y1, Y2), loading, X.target, covariate.shift = TRUE)
#' out <- measure_instability(mm, gen.size=5)
#' out$measure
measure_instability <- function(object, delta=0, gen.size=500, threshold=0, alpha.thres=0.01){
  spnorm.Gamma.diff = rep(0, gen.size)
  l2norm.gamma.diff = rep(0, gen.size)

  #############################
  ##### Generate Samples ######
  #############################
  gen.dim = object$L*(object$L+1)/2
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
    gen.matrix = matrix(NA, nrow=object$L, ncol=object$L)
    gen.matrix[lower.tri(gen.matrix, diag=TRUE)] = gen.samples[g,]
    for(l in 1:object$L){
      for(k in 2:object$L){
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
#' @param gen.size The generating sample size (Default = 500)
#' @param step_delta The step size of searching delta (Default = 0.1)
#' @param MAX_iter Maximum of iterations for searching (Default = 100)
#' @param verbose Print information about delta and reward (Default = `FALSE`)
#'
#' @return
#' \item{delta}{The data-dependent ridge penalty}
#' \item{reward.ratio}{The ratio of penalized reward over non-penalized reward}
#' @export
#' @examples
#' ## The problem is low-dimensional for testings
#' ## heterogenous data and covariates shift
#' X1 = sample_data$X1
#' X2 = sample_data$X2
#' Y1 = sample_data$Y1
#' Y2 = sample_data$Y2
#' X.target = sample_data$X.target
#'
#' ## loading
#' loading = rep(0, 5) # dimension p=5
#' loading[5] = 1
#'
#' ## call
#' mm <- Maximin(list(X1, X2), list(Y1, Y2), loading, X.target)
#' out <- decide_delta(mm, gen.size=5)
#' out$delta
#' out$reward.ratio
decide_delta <- function(object, gen.size=500, step_delta=0.1, MAX_iter=100, verbose=FALSE){
  #######################################
  ############### STEP-1 ################
  ## Compute Measure square on delta=0 ##
  #######################################
  measure = measure_instability(object, delta=0, gen.size=gen.size, threshold=0, alpha.thres=0.01)$measure
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
    solution = opt.weight(object$Gamma.prop, delta)
    reward = solution$reward
    if(verbose){
      print(paste0("The picked delta is ", round(delta,4)))
      print(paste0("Reward Ratio is ", round(reward / reward0, 4)))
    }
    out <- list(delta = delta,
                reward.ratio = reward/reward0)
  }
  return(out)
}
