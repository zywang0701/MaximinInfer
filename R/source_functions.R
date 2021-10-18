opt.weight<-function(Gamma,delta,report.reward=TRUE){
  ## Purpose: Compute Ridge-type weight vector
  ## Returns: weight:  the minimizer \eqn{\gamma}
  ##          reward:  the value of penalized reward
  ## ----------------------------------------------------
  ## Arguments: Gamma: regression covariance matrix, of dimension \eqn{L} x \eqn{L}
  ##            delta the ridge penalty level, non-positive.
  ##            report.reward the reward is computed or not (Default = `TRUE`)
  ## ----------------------------------------------------

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

Gamma.shift<-function(plug.in,X.l,X.k,omega.l,omega.k,Y.l,Y.k,Pred.l,Pred.k){
  ## Purpose: Bias correction for initial estimator of Gamma when covariate shifts
  ## Returns: est:       The proposed bias-corrected estimator of Gamma target
  ##          proj.lk:   The projection direction of index (l, k)
  ##          proj.kl:   The projection direction of index (k, l)
  ## ----------------------------------------------------------------------
  ## Arguments: plug.in: Initial estimator of Gamma target, of dimension \eqn{L} x \eqn{L}
  ##            X.l:     Design matrix of label \eqn{l} in training data, of dimension \eqn{n_l} x \eqn{p}
  ##            X.k:     Design matrix of label \eqn{k} in training data, of dimension \eqn{n_k} x \eqn{p}
  ##            omega.l: The l-th column of Omega matrix
  ##            omega.k: The k-th column of Omega matrix
  ##            Y.l:     Outcome vector of label \eqn{l} in training data, of length \eqn{n_l}
  ##            Y.k:     Outcome vector of label \eqn{k} in training data, of length \eqn{n_k}
  ##            Pred.l:  Predicted outcome vector of label \eqn{l} in training data
  ##            Pred.k:  Predicted outcome vector of label \eqn{k} in training data
  ## ----------------------------------------------------------------------
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

cov.inner.shift<-function(Var.vec,l1,k1,l2,k2,X.l1,X.k1,X.l2,X.k2,Pred.mat.target,Proj.array){
  ## Purpose: Estimate the covariance between the pi(l1,k1) entry and pi(l2,k2) entry
  ##          when the covariance for the target distribution is unknown
  ## Returns: covariance between the pi(l1, k1) entry and pi(l2, k2) entry
  ## ----------------------------------------------------------------------
  ## Arguments: Var.vec:         Variance of residuals in groups, of length \eqn{L}
  ##            l1:              Index l1
  ##            k1:              Index k1
  ##            l2:              Index l2
  ##            k2:              Index k2
  ##            X.l1:            Design matrix of group \eqn{l1} in training data
  ##            X.k1:            Design matrix of group \eqn{k1} in training data
  ##            X.l2:            Design matrix of group \eqn{l2} in training data
  ##            X.k2:            Design matrix of group \eqn{k2} in training data
  ##            Pred.mat.target: Predicted outcome matrix for target design matrix
  ##                             L fitted coefficients, of dimension \eqn{n.target} x \eqn{L}
  ##            Proj.array:      Projection directions, of dimension \eqn{L} x \eqn{L} x \eqn{p}
  ## ----------------------------------------------------------------------
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

cov.inner.shift.known<-function(Var.vec,l1,k1,l2,k2,X.l1,X.k1,X.l2,X.k2,Proj.array){
  ## Purpose: Estimate the covariance between the pi(l1,k1) entry and pi(l2,k2) entry
  ##          when the covariance for the target distribution is known
  ## Returns: covariance between the pi(l1, k1) entry and pi(l2, k2) entry
  ## ----------------------------------------------------------------------
  ## Arguments: Var.vec:         Variance of residuals in groups, of length \eqn{L}
  ##            l1:              Index l1
  ##            k1:              Index k1
  ##            l2:              Index l2
  ##            k2:              Index k2
  ##            X.l1:            Design matrix of group \eqn{l1} in training data
  ##            X.k1:            Design matrix of group \eqn{k1} in training data
  ##            X.l2:            Design matrix of group \eqn{l2} in training data
  ##            X.k2:            Design matrix of group \eqn{k2} in training data
  ##            Proj.array:      Projection directions, of dimension \eqn{L} x \eqn{L} x \eqn{p}
  ## ----------------------------------------------------------------------
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

Lasso <- function(X, y, lambda = NULL, intercept = TRUE) {
  p <- ncol(X)
  n <- nrow(X)

  htheta <- if (is.null(lambda)) {
    lambda <- sqrt(qnorm(1 - (0.1 / p)) / n)
    outLas <- slim(X, y, lambda = lambda, method = "lq", q = 2,
                   verbose = FALSE)
    # Objective : sqrt(RSS/n) + lambda * penalty
    c(as.vector(outLas$intercept), as.vector(outLas$beta))
  } else if (lambda == "CV") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else if (lambda == "CV.min") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.min))
  } else if (lambda == "scalreg") {
    Xc <- if (intercept) {
      cbind(rep(1, n), X)
    } else {
      X
    }
    outLas <- scalreg(Xc, y)
    # return object
    if (intercept) {
      outLas$coefficients
    } else {
      # add a coefficient for the (not estimated) intercept b/c of implementation
      c(0, outLas$coefficients)
    }
  } else {
    outLas <- glmnet(X, y, family = "gaussian", alpha = 1,
                     intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = lambda))
  }

  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}

getmode <- function(v) {
  tbl <- table(v)
  if (all(tbl == 1)) {
    median(v)
  } else {
    as.numeric(names(which.max(tbl)))
  }
}

Direction_fixedtuning_lin<-function(X,loading,mu=NULL){
  pp<-ncol(X)
  n<-nrow(X)
  if(is.null(mu)){
    mu<-sqrt(2.01*log(pp)/n)
  }
  loading.norm<-sqrt(sum(loading^2))
  if (loading.norm==0){
    H <- cbind(loading, diag(1, pp))
  }else{
    H <- cbind(loading / loading.norm, diag(1, pp))
  }
  v<-Variable(pp+1)
  obj<-1/4*sum((X%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
  prob<-Problem(Minimize(obj))
  result<-solve(prob)
  opt.sol<-result$getValue(v)
  cvxr_status<-result$status
  direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  returnList <- list("proj"=direction)
  return(returnList)
}

Direction_searchtuning_lin<-function(X,loading,mu=NULL, resol = 1.5, maxiter = 10){
  pp<-ncol(X)
  n<-nrow(X)
  tryno = 1;
  opt.sol = rep(0,pp+1);
  lamstop = 0;
  cvxr_status = "optimal";
  mu = sqrt(2.01*log(pp)/n);
  #mu.initial= mu;
  while (lamstop == 0 && tryno < maxiter){
    ###### This iteration is to find a good tuning parameter
    lastv = opt.sol;
    lastresp = cvxr_status;
    loading.norm<-sqrt(sum(loading^2))
    if (loading.norm==0){
      H <- cbind(loading, diag(1, pp))
    }else{
      H <- cbind(loading / loading.norm, diag(1, pp))
    }
    v<-Variable(pp+1)
    obj<-1/4*sum((X%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
    prob<-Problem(Minimize(obj))
    result<-solve(prob)
    cvxr_status<-result$status

    if(tryno==1){
      if(cvxr_status=="optimal"){
        incr = 0;
        mu=mu/resol;
        opt.sol<-result$getValue(v) ### we should move this line from above to here
        temp.vec<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
        initial.sd<-sqrt(sum((X%*% temp.vec)^2)/(n)^2)*loading.norm ##what's this?
        temp.sd<-initial.sd
      }else{
        incr = 1;
        mu=mu*resol;
      }
    }else{
      if(incr == 1){ ### if the tuning parameter is increased in the last step
        if(cvxr_status=="optimal"){
          lamstop = 1;
          opt.sol<-result$getValue(v)
        }else{
          mu=mu*resol;
        }
      }else{
        if(cvxr_status=="optimal"&&temp.sd<3*initial.sd){ ##Why this condition on sd?
          mu = mu/resol;
          opt.sol<-result$getValue(v)
          temp.vec<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
          temp.sd<-sqrt(sum((X%*% temp.vec)^2)/(n)^2)*loading.norm
          #print(temp.sd)
        }else{
          mu=mu*resol;
          opt.sol=lastv;
          lamstop=1;
          tryno=tryno-1
        }
      }
    }
    tryno = tryno + 1;
  }
  direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  step<-tryno-1
  returnList <- list("proj"=direction,
                     "step"=step)
  return(returnList)
}

index.map<-function(L,l,k){
  return((2*L-k)*(k-1)/2+l)
}

proj.direction<-function(Xc,loading,maxiter=6,resol=1.25){
  n<-dim(Xc)[1]
  p<-dim(Xc)[2]
  loading.norm<-sqrt(sum(loading^2))
  sigma.hat <- (1/n)*(t(Xc)%*%Xc);
  if ((n>=6*p)){
    tmp <- eigen(sigma.hat)
    tmp <- min(tmp$values)/max(tmp$values)
  }else{
    tmp <- 0
  }
  if ((n>=6*p)&&(tmp>=1e-4)){
    direction <- solve(sigma.hat)%*%loading
  }else{
    step.vec<-rep(NA,3)
    for(t in 1:3){
      index.sel<-sample(1:n,size=ceiling(0.5*min(n,p)), replace=FALSE)
      Direction.Est.temp<-Direction_searchtuning_lin(Xc[index.sel,],loading,mu=NULL, resol, maxiter)
      step.vec[t]<-Direction.Est.temp$step
    }
    step<-getmode(step.vec)
    Direction.Est<-Direction_fixedtuning_lin(Xc,loading,mu=sqrt(2.01*log(p)/n)*resol^{-(step-1)})
    while(is.na(Direction.Est) || length(Direction.Est$proj)==0){
      step<-step-1
      Direction.Est <- Direction_fixedtuning_lin(Xc, loading, mu = sqrt(2.01 * log(p) / n) * resol^{-(step - 1)})
    }

    direction<-loading.norm*Direction.Est$proj
  }
  return(direction)
}

#'
#'
#' @param
#' @param
#' @param
#' @param
#'
#' @return
#' \item{delta}{}
decide_delta <- function(Gamma, step_delta=0.1, MAX_iter=100, verbose=FALSE){
  ## Purpose: Decide ridge penalty data-dependently
  ## Returns: The data-dependent ridge penalty
  ## ----------------------------------------------------------------------
  ## Arguments: Gamma Weight Matrix
  ##            step_delta The step size of searching delta (default = 0.1)
  ##            MAX_iter Maximum of iterations for searching
  ##            verbose Print information about delta and reward (default = `FALSE`)
  ## ----------------------------------------------------------------------
  L = dim(Gamma)[1]
  min_eigen = min(eigen(Gamma)$values)
  max_eigen = max(eigen(Gamma)$values)
  solution0 = opt.weight(Gamma, 0)
  reward0 = solution0$reward
  delta_min = 0.1 * reward0 * (2 * L) / (L - 1)
  if(delta_min > 2){
    if(verbose) print(paste("delta path starts from", round(delta_min, 4), "which exceeds our maximum limit 2"))
    delta = 2
  }else{
    delta = delta_min
    solution = opt.weight(Gamma, delta)
    reward = solution$reward
    i = 1
    while(reward >= 0.9 * reward0){
      delta_new = delta + step_delta
      solution = opt.weight(Gamma, delta_new)
      reward = solution$reward
      if(reward <= 0.9 * reward0) break
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
  }

  if((min_eigen + delta) < 0.5) warning("Fail to find a suitable delta, the estimator may be not stable enough.")

  solution = opt.weight(Gamma, delta)
  reward = solution$reward
  if(verbose){
    print(paste0("The picked delta is ", round(delta,4)))
    print(paste0("Reward Ratio is ", round(reward / reward0, 4)))
    print(paste0("Minimum Eigenvalue plus delta = ", round(min_eigen + delta, 4)))
  }
  return(delta)
}
