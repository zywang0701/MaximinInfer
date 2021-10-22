#' Maximin Inference
#' @description
#' `infer` is a generic function for inference for Maximin model.
#' @param object a "Maximin" object
#' @param ... additional arguments affecting inference
#'
#' @export
infer <- function(object, ...){
  UseMethod("infer", object)
}

#' Inference method for Maximin
#' @description Point estimator and Confidence interval based on Maximin object
#'
#' @param object Object of class inheriting from "Maximin"
#' @param gen.size The generated sample size
#' @param delta The ridge penalty. If set as negative value, the penalty is decided data-dependently. (Default = -1)
#' @param threshold Should generated samples be filter or not? If 0, do not filter;
#' if 1, use chi-square threshold to filter;
#' if 2, use normal threshold to filter. (Default = 2)
#' @param alpha confidence value to select generated samples
#' @param ... further arguments passed
#'
#' @return
#' \item{delta}{The ridge penalty used}
#' \item{weight}{The weight vector for groups, of length \eqn{L}}
#' \item{point}{The point estimator of the linear contrast}
#' \item{mm.effect}{The aggregated maximin effect (coefficients), of length \eqn{p} or \eqn{p+1}}
#' \item{CI}{Confidence interval for the linear contrast}
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
#' mm <- Maximin(list(X1, X2), list(Y1, Y2), X.target, loading, covariate.shift = TRUE)
#' mmInfer <- infer(mm)
#' }
infer.Maximin <- function(object, gen.size=500, delta=-1, threshold=2, alpha=0.01,...){
  if(delta<0) delta = decide_delta(object$Gamma.prop, step_delta=0.1)
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
    gen.samples = matrix(mvrnorm(gen.size, mu=object$gen.mu, Sigma=object$gen.Cov),
                         nrow=gen.size, ncol=gen.dim)
  }
  if(threshold==1){
    # 1 stands for chi square threshold
    UDV_list = svd(object$gen.Cov)
    U = UDV_list$u
    D = UDV_list$d
    V = UDV_list$v
    D.sqrt = sqrt(D)
    gen.Cov.sqrt = U %*% diag(D.sqrt) %*% t(V)
    thres = 1.0 * qchisq((1-alpha), df=gen.dim)
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
  if(threshold == 2){
    thres = qnorm(1-alpha/(gen.dim*2))
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
  CI.original = cbind(gen.est[, 1]-1.96*gen.est[, 2], gen.est[,1]+1.96*gen.est[,2])
  CI = na.omit(CI.original)
  uni = Intervals(CI)
  CI.union = as.matrix(interval_union(uni))
  colnames(CI.union) <- c("lower", "upper")

  out = list(delta = delta,
             weight = weight,
             point = point,
             mm.effect = mm.effect,
             CI = CI.union)
  return(out)
}
