#' Inference for a linear combination of regression coefficients in high dimensional linear regression
#'
#' @description
#' Computes the bias corrected estimator of the linear combination of regression coefficients and
#' the corresponding standard error. It also constructs the confidence interval for the linear combination.
#'
#' @param X Design matrix, of dimension \eqn{n} x \eqn{p}
#' @param y Outcome vector, of length \eqn{n}
#' @param loading Loading, of length \eqn{p}
#' @param intercept Should intercept(s) be fitted or not (default = \code{TRUE})
#' @param init.Lasso Initial LASSO estimator of the regression vector (default = \code{TRUE})
#' @param lambda The tuning parameter in the construction of LASSO estimator of the regression vector (default = \code{NULL})
#' @param mu The dual tuning parameter used in the construction of the projection direction (default = \code{NULL})
#' @param step The step size used to compute \code{mu}; if set to \code{NULL} it is computed to be
#' the number of steps (\code{maxiter}) to obtain the smallest \code{mu} such that the dual optimization problem
#' for constructing the projection direction converges (default = \code{NULL})
#' @param resol The factor by which \code{mu} is increased/decreased to obtain the smallest \code{mu}
#' such that the dual optimization problem for constructing the projection direction converges (default = 1.5)
#' @param maxiter Maximum number of steps along which \code{mu} is increased/decreased to obtain the smallest \code{mu}
#' such that the dual optimization problem for constructing the projection direction converges (default = 6)
#'
#' @return
#' \item{prop.est}{The bias-corrected estimator for the linear combination of regression coefficients}
#' \item{se}{The standard error of the bias-corrected estimator}
#' \item{proj}{The projection direction, of length \eqn{p}}
#' \item{plug.in}{The plug-in LASSO estimator for the linear combination}
#'
#' @examples
#' \dontrun{
#' n <- 90
#' p <- 200
#' A1gen <- function(rho,p){
#'  A1=matrix(0,p,p)
#'   for(i in 1:p){
#'    for(j in 1:p){
#'      A1[i,j] <- rho^(abs(i-j))
#'    }
#'   }
#'  A1
#' }
#' mu <- rep(0,p)
#' rho <- 0.5
#' Cov <- (A1gen(rho,p))/2
#' beta <- rep(0,p)
#' beta[1:10] <- c(1:10)/5
#' X <- MASS::mvrnorm(n,mu,Cov)
#' y <- X%*%beta + rnorm(n)
#' loading <- c(1,rep(0,(p-1)))
#' Est <- LF(X = X, y = y, loading = loading, intercept = TRUE)
#' }
LF<-function(X,y,loading,intercept=TRUE,init.Lasso=NULL,lambda=NULL,mu=NULL,step=NULL,resol = 1.5,maxiter=10){
  ### Option 1: search tuning parameter with steps determined by the ill conditioned case (n=p/2)
  ### Option 2: search tuning parameter with maximum 10 steps.
  ### Option 3: fixed tuning parameter and this is not recommended without exploring the tuning parameter selection
  xnew<-loading
  p <- ncol(X);
  n <- nrow(X);
  n_y <- length(y)

  if(n_y!=n)
  {
    stop("Error: Check dimensions of X and y")
  }
  else
  {
    data = na.omit(data.frame(y,X))
    X <- as.matrix(data[,-1])
    y <- as.vector(data[,1])
    p <- ncol(X)
    n <- nrow(X)

    # implement the correction of the initial estimator
    # set up the randomization step
    if (intercept==TRUE){
      Xc <- cbind(rep(1,n),X);
      pp <- (p+1);
    } else {
      Xc <- X;
      pp <- p
    }
    if(is.null(init.Lasso)){
      # implement a lasso algorithm to get beta and sigma
      init.Lasso<-Initialization.step(X,y,lambda,intercept)
      htheta<-init.Lasso$lasso.est
    }else{
      htheta<-init.Lasso
    }
    sparsity <- sum(abs(htheta) > 0.001)
    sd.est <- sqrt(sum((y - Xc %*% htheta)^2) / n)

    # compute the initial estimator
    if(intercept==TRUE){
      loading=rep(0,pp)
      loading[1]=1
      loading[-1]=xnew
    }else{
      loading=xnew
    }
    loading.norm<-sqrt(sum(loading^2))
    lasso.plugin<-sum(loading*htheta)

    #############################################
    ############# Correction step ###############

    if ((n>=6*p)){
      sigma.hat <- (1/n)*(t(Xc)%*%Xc);
      tmp <- eigen(sigma.hat)
      tmp <- min(tmp$values)/max(tmp$values)
    }else{
      tmp <- 0
    }
    sigma.hat <- (1/n)*(t(Xc)%*%Xc);
    if ((n>=6*p)&&(tmp>=1e-4)){
      direction <- solve(sigma.hat)%*%loading/loading.norm
    }else{
      if(is.null(step)){
        step.vec<-rep(NA,3)
        for(t in 1:3){
          index.sel<-sample(1:n,size=ceiling(0.5*min(n,p)), replace=FALSE)
          Direction.Est.temp<-Direction_searchtuning_lin(Xc[index.sel,],loading,mu=NULL, resol, maxiter)
          step.vec[t]<-Direction.Est.temp$step
        }
        step<-getmode(step.vec)
      }
      Direction.Est<-Direction_fixedtuning_lin(Xc,loading,mu=sqrt(2.01*log(pp)/n)*resol^{-(step-1)})
      while(is.na(Direction.Est) || length(Direction.Est$proj)==0){
        step<-step-1
        Direction.Est <- Direction_fixedtuning_lin(Xc, loading, mu = sqrt(2.01 * log(pp) / n) * resol^{-(step - 1)})
      }
      direction<-Direction.Est$proj
    }
    correction = t(Xc%*%direction)%*%(y - Xc%*%htheta)/n;
    debias.est=lasso.plugin+correction*loading.norm
    se<-sd.est*sqrt(sum((Xc%*%direction)^2)/(n)^2)*loading.norm

    returnList <- list("prop.est" = debias.est,
                       "se" = se,
                       "proj"=direction,
                       "plug.in"=lasso.plugin
    )
    return(returnList)
  }
}


#' Compute Ridge-type weight vector
#'
#' @description \eqn{\argmin_{\gamma \in \Delta} \gamma^T (\Gamma + \delta * I) \gamma}.
#'
#' @param Gamma regression covariance matrix, of dimension \eqn{L} x \eqn{L}.
#' @param delta the ridge penalty level, non-positive.
#' @param report.reward the reward is computed or not (Default = \code{TRUE})
#'
#' @return
#' \item{weight}{the minimizer \eqn{\gamma}}
#' \item{reward}{the value of penalized reward}
#'
opt.weight<-function(Gamma,delta,report.reward=TRUE){
  L<-dim(Gamma)[2]
  opt.weight<-rep(NA, L)
  opt.reward<-NA
  # Problem definition
  v<-Variable(L)
  Diag.matrix<-diag(eigen(Gamma)$values)
  for(ind in 1:L){
    Diag.matrix[ind,ind]<-max(Diag.matrix[ind,ind],0.001)
  }
  Gamma.positive<-eigen(Gamma)$vectors%*%Diag.matrix%*%t(eigen(Gamma)$vectors)
  objective <- Minimize(quad_form(v,Gamma.positive+diag(delta,L)))
  constraints <- list(v >= 0, sum(v)== 1)
  prob.weight<- Problem(objective, constraints)
  if(is_dcp(prob.weight)){
    result<- solve(prob.weight)
    opt.status<-result$status
    opt.sol<-result$getValue(v)
    for(l in 1:L){
      opt.weight[l]<-opt.sol[l]*(abs(opt.sol[l])>10^{-8})
    }
  }
  if(report.reward){
    v<-Variable(L)
    objective<-Minimize(2*t(v)%*%Gamma.positive%*%opt.weight-t(opt.weight)%*%Gamma.positive%*%opt.weight)
    constraints<-list(v >= 0, sum(v)== 1)
    delta.optim<-Problem(objective, constraints)
    result<- solve(delta.optim)
    opt.reward<-result$value
    returnList <- list("weight" = opt.weight,
                       "reward" = opt.reward)
  }else{
    returnList <- list("weight" = opt.weight)
  }
  return(returnList)
}

#' Bias correction for initial estimator of Gamma target
#'
#' @param plug.in Initial estimator of Gamma target, of dimension \eqn{L} x \eqn{L}
#' @param X.l Design matrix of label \eqn{l} in training data, of dimension \eqn{n_l} x \eqn{p}
#' @param X.k Design matrix of label \eqn{k} in training data, of dimension \eqn{n_k} x \eqn{p}
#' @param omega.l The l-th column of Omega matrix
#' @param omega.k The k-th column of Omega matrix
#' @param Y.l Outcome vector of label \eqn{l} in training data, of length \eqn{n_l}
#' @param Y.k Outcome vector of label \eqn{k} in training data, of length \eqn{n_k}
#' @param Pred.l Predicted outcome vector of label \eqn{l} in training data
#' @param Pred.k Predicted outcome vector of label \eqn{k} in training data
#'
#' @return
#' \item{est}{The proposed bias-corrected estimator of Gamma target}
#' \item{proj.lk}{The projection direction of index (l, k)}
#' \item{proj.kl}{The projection direction of index (k, l)}
Gamma.shift<-function(plug.in,X.l,X.k,omega.l,omega.k,Y.l,Y.k,Pred.l,Pred.k){
  u.lk<-proj.direction(X.l,omega.k)
  u.kl<-proj.direction(X.k,omega.l)
  n.k<-nrow(X.k)
  n.l<-nrow(X.l)
  prop.est<-plug.in+t(u.kl)%*%t(X.k)%*%(Y.k-Pred.k)/n.k+t(u.lk)%*%t(X.l)%*%(Y.l-Pred.l)/n.l
  returnList <- list("est" = prop.est,
                     "proj.lk" = u.lk,
                     "proj.kl" = u.kl
                     )
  return(returnList)
}

#' Estimate the covariance between the pi(l1,k1) entry and pi(l2,k2) entry
#' when the covariance for the targeted distribution is unknown
#'
#' @param Var.vec Variance of residuals in groups, of length \eqn{L}
#' @param l1 Index l1
#' @param k1 Index k1
#' @param l2 Index l2
#' @param k2 Index k2
#' @param X.l1 Design matrix of group \eqn{l1} in training data
#' @param X.k1 Design matrix of group \eqn{k1} in training data
#' @param X.l2 Design matrix of group \eqn{l2} in training data
#' @param X.k2 Design matrix of group \eqn{k2} in training data
#' @param Pred.mat.target Predicted outcome matrix for target design matrix
#' using L fitted coefficients, of dimension \eqn{n.target} x \eqn{L}
#' @param Proj.array Projection directions, of dimension \eqn{L} x \eqn{L} x \eqn{p}
#'
#' @return covariance between the pi(l1, k1) entry and pi(l2, k2) entry
cov.inner.shift<-function(Var.vec,l1,k1,l2,k2,X.l1,X.k1,X.l2,X.k2,Pred.mat.target,Proj.array){
  N<-dim(Pred.mat.target)[1]
  Sigma.est.l1<-(1/dim(X.l1)[1])*(t(X.l1)%*%X.l1)
  Sigma.est.k1<-(1/dim(X.k1)[1])*(t(X.k1)%*%X.k1)
  var1<-0
  if(l2==l1){
    var1<-var1+Var.vec[l1]*Proj.array[l1,k1,]%*%Sigma.est.l1%*%Proj.array[l2,k2,]/dim(X.l1)[1]
  }
  if(k2==l1){
    var1<-var1+Var.vec[l1]*Proj.array[l1,k1,]%*%Sigma.est.l1%*%Proj.array[k2,l2,]/dim(X.l1)[1]
  }
  if(l2==k1){
    var1<-var1+Var.vec[k1]*Proj.array[k1,l1,]%*%Sigma.est.k1%*%Proj.array[l2,k2,]/dim(X.k1)[1]
  }
  if(k2==k1){
    var1<-var1+Var.vec[k1]*Proj.array[k1,l1,]%*%Sigma.est.k1%*%Proj.array[k2,l2,]/dim(X.k1)[1]
  }
  var2<-mean((diag(Pred.mat.target[,k1]%*%t(Pred.mat.target[,l1]))-mean(Pred.mat.target[,k1]*Pred.mat.target[,l1]))*(diag(Pred.mat.target[,k2]%*%t(Pred.mat.target[,l2]))-mean(Pred.mat.target[,k2]*Pred.mat.target[,l2])))
  var<-var1+var2/N
  return((var))
}


#' Estimate the covariance between the pi(l1,k1) entry and pi(l2,k2) entry
#' when the covariance for the targeted distribution is known
#'
#' @param Var.vec Variance of residuals in groups, of length \eqn{L}
#' @param l1 Index l1
#' @param k1 Index k1
#' @param l2 Index l2
#' @param k2 Index k2
#' @param X.l1 Design matrix of group \eqn{l1} in training data
#' @param X.k1 Design matrix of group \eqn{k1} in training data
#' @param X.l2 Design matrix of group \eqn{l2} in training data
#' @param X.k2 Design matrix of group \eqn{k2} in training data
#' @param Proj.array Projection directions, of dimension \eqn{L} x \eqn{L} x \eqn{p}
#'
#' @return covariance between the pi(l1, k1) entry and pi(l2, k2) entry
cov.inner.shift.known<-function(Var.vec,l1,k1,l2,k2,X.l1,X.k1,X.l2,X.k2,Proj.array){
  Sigma.est.l1<-(1/dim(X.l1)[1])*(t(X.l1)%*%X.l1)
  Sigma.est.k1<-(1/dim(X.k1)[1])*(t(X.k1)%*%X.k1)
  var1<-0
  if(l2==l1){
    var1<-var1+Var.vec[l1]*Proj.array[l1,k1,]%*%Sigma.est.l1%*%Proj.array[l2,k2,]/dim(X.l1)[1]
  }
  if(k2==l1){
    var1<-var1+Var.vec[l1]*Proj.array[l1,k1,]%*%Sigma.est.l1%*%Proj.array[k2,l2,]/dim(X.l1)[1]
  }
  if(l2==k1){
    var1<-var1+Var.vec[k1]*Proj.array[k1,l1,]%*%Sigma.est.k1%*%Proj.array[l2,k2,]/dim(X.k1)[1]
  }
  if(k2==k1){
    var1<-var1+Var.vec[k1]*Proj.array[k1,l1,]%*%Sigma.est.k1%*%Proj.array[k2,l2,]/dim(X.k1)[1]
  }
  var<-var1
  return((var))
}
