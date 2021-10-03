getmode <- function(v) {
  tbl <- table(v)
  if (all(tbl == 1)) {
    median(v)
  } else {
    as.numeric(names(which.max(tbl)))
  }
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

Initialization.step <- function(X, y, lambda = NULL, intercept = FALSE) {
  n <- nrow(X)
  col.norm <- 1 / sqrt((1 / n) * diag(t(X) %*% X))
  Xnor <- X %*% diag(col.norm)

  ### Call Lasso
  htheta <- Lasso(Xnor, y, lambda = lambda, intercept = intercept)

  ### Calculate return quantities
  if (intercept == TRUE) {
    Xb <- cbind(rep(1, n), Xnor)
    col.norm <- c(1, col.norm)
  } else {
    Xb <- Xnor
  }
  sparsity <- sum(abs(htheta) > 0.001)
  sd.est <- sqrt(sum((y - Xb %*% htheta)^2) / n)
  htheta <- htheta * col.norm
  returnList <- list("lasso.est" = htheta,
                     "sigma" = sd.est,
                     "sparsity" = sparsity)
  return(returnList)
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

decide_delta <- function(Gamma, step_delta=0.1, MAX_iter=100, verbose=TRUE){
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