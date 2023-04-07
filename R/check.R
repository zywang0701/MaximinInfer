check.args <- function(Xlist=NULL, Ylist=NULL, loading.mat=NULL, X0=NULL, cov.shift=NULL,
                       cov0=NULL, intercept=NULL, intercept.loading=NULL, lambda=NULL, verbose=NULL){
  if(is.null(Xlist) || (!is.list(Xlist))) stop("Xlist must be a list")
  if(is.null(Ylist) || (!is.list(Ylist))) stop("Ylist must be a list")
  if(length(Xlist) != length(Ylist)) stop("The number of groups in Xlist and Ylist must match")
  if(length(unique(sapply(Xlist, FUN=ncol))) != 1) stop("each group should have the same dimension of covariates")
  p = unique(sapply(Xlist, FUN=ncol))
  if(any(sapply(Xlist, nrow)!= sapply(Ylist, length))) stop("The group should match in Xlist and Ylist:
                                                       they should have the same number of observations for each group")
  if(is.null(loading.mat) || !is.numeric(loading.mat)) stop("loading must be a numeric matrix")
  if(p != nrow(loading.mat)) stop("ncol(X) and nrow(loading) must match")
  if(!is.null(X0)){
    if(!is.numeric(X0)|| ncol(X0)!=p) stop("X0 must be a numeric matrix with the same dimension as source covariates.")
  }
  if(is.null(cov.shift) || !is.logical(cov.shift) || length(cov.shift)!=1 ) stop("cov.shift must be a Boolean")
  if(!is.null(cov0)){
    if(!isSymmetric.matrix(cov0)) stop('cov0 must be a symmetric matrix.')
    if(any(eigen(cov0)$values < 0)) stop('cov0 must be positive semi-definite.')
    if(nrow(cov0)!=p) stop('The dimension of cov0 should match with the dimension of covariates')
  }
  if(is.null(intercept) || length(intercept)!=1 || !is.logical(intercept)){
    stop("intercept must be a Boolean")
  }
  if(is.null(intercept.loading) || length(intercept.loading)!=1 || !is.logical(intercept.loading)){
    stop("intercept.loading must be a Boolean")
  }
  if(!is.null(lambda)){
    if(length(lambda)!=1) stop("lambda must be length of 1")
    if(!is.numeric(lambda)) if(!(lambda %in% c("CV.min","CV"))) stop("lambda must be a number, except from 'CV.min','CV'.")
  }
  if(is.null(verbose) || !is.logical(verbose) || length(verbose)!=1 ) stop("verbose must be a Boolean")
}
