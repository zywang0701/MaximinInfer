index.map <- function(L, l, k){
  return((2*L-k)*(k-1)/2+l)
}

relevant.funs <- function(intercept=TRUE){
  ## high dimensional linear regression ##
  train.fun <- function(X, y, lambda=NULL){
    if(is.null(lambda)) lambda = 'CV.min'
    p = ncol(X)
    htheta <- if (lambda == "CV.min") {
      outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                          intercept = intercept, standardize = T)
      as.vector(coef(outLas, s = outLas$lambda.min))
    } else if (lambda == "CV") {
      outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                          intercept = intercept, standardize = T)
      as.vector(coef(outLas, s = outLas$lambda.1se))
    } else {
      outLas <- glmnet(X, y, family = "gaussian", alpha = 1,
                       intercept = intercept, standardize = T)
      as.vector(coef(outLas, s = lambda))
    }
    if(intercept==FALSE) htheta = htheta[2:(p+1)]

    return(list(lasso.est = htheta))
  }
  ## predict ##
  pred.fun <- function(X, htheta){
    X = as.matrix(X)
    if((length(htheta) - ncol(X)) == 1) X = cbind(1,X)
    pred = as.vector(X%*%htheta)
    return(pred)
  }
  ## variance of residual ##
  dev.fun <- function(pred, y, sparsity=0){
    pred = as.vector(pred); y = as.vector(y)
    n = length(y)
    sigmasq.hat = sum((y - pred)^2) / max(0.7*n, n-sparsity)
    return(sigmasq.hat)
  }
  return(list(train.fun = train.fun,
              pred.fun = pred.fun,
              dev.fun = dev.fun))
}


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

gensamples <- function(gen.mu, gen.Cov, gen.size=500, threshold=0, alpha.thres=0.01){
  gen.dim = length(gen.mu)
  if(threshold==0){
    thres = qnorm(1-alpha.thres/(gen.dim*2))
    gen.samples = matrix(0, nrow=gen.size, ncol=gen.dim)
    n.picked = 0 # records how many samples are picked
    while(n.picked < gen.size){
      S = mvrnorm(1, mu=rep(0, gen.dim), Sigma=gen.Cov)
      if(max(abs(S / sqrt(diag(gen.Cov)))) <= thres){
        n.picked = n.picked + 1
        gen.samples[n.picked, ] = gen.mu + S
      }
    }
  }
  if(threshold==1){
    # 1 stands for chi square threshold
    UDV_list = svd(gen.Cov)
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
        gen.samples[n.picked, ] = gen.mu + gen.Cov.sqrt %*% Z
      }
    }
  }
  if(threshold==2){
    gen.samples = matrix(mvrnorm(gen.size, gen.mu, Sigma=gen.Cov),
                         nrow=gen.size, ncol=gen.dim)
  }
  return(gen.samples)
}



