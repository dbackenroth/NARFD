library(nloptr)

# objective for estimating FPC spline coefficients
Objective3 <- function(x, X, y, pmat, lambda){
  beta <- x
  mu <- as.numeric(X %*% beta)
  log.lik <- sum(y * log(mu+1e-16) - mu)
  pen.term <- lambda * t(beta) %*% pmat %*% beta
  objective <- log.lik - pen.term
  return(-objective)
}

# objective for estimating scores
Objective4 <- function(x, X, y){
  beta <- x
  mu <- as.numeric(X %*% beta)
  log.lik <- sum(y * log(mu+1e-16) - mu)
  return(-log.lik)
}

Gradient3 <- function(x, X, y, pmat, lambda){
  beta <- x
  mu <- as.numeric(X %*% beta)
  grad <- as.matrix(colSums(X * (y / (mu+1e-16) - 1))) - 2 * lambda * pmat %*% beta
  return(-grad)
}

Gradient4 <- function(x, X, y){
  beta <- x
  mu <- as.numeric(X %*% beta)
  grad <- as.matrix(colSums(X * (y / (mu+1e-16) - 1)))
  return(-grad)
}

# non-negative regression without regularization
FitNL.noR <- function(X, y, init){
  LB <- rep(0, length(init))
  UB <- rep(Inf, length(init))
  l <- list()
  obj <- list()
  # gradient descent with supplied initial values
  nlopt.fit <- nloptr(x0=init, eval_f=Objective4, 
                      eval_grad=Gradient4, lb=LB, 
                      ub=UB, opts=list("algorithm"="NLOPT_LD_LBFGS", 
                                       "xtol_rel"=1e-04), 
                      X=X, y=y)
  l[[1]] <- nlopt.fit$solution
  obj[[1]] <- Objective4(l[[1]], X, y)
  # nnlm without initial values
  l[[2]] <- nnlm(X,
                   as.matrix(y), 
                   loss="mkl",
                   check.x=F)$coefficients[, 1]
  obj[[2]] <- Objective4(l[[2]], X, y)
  # nnlm with initial values
  l[[3]] <- nnlm(X, 
                    as.matrix(y), loss="mkl", check.x=F, 
                 init=init)$coefficients[, 1]
  obj[[3]] <- Objective4(l[[3]], X, y)
  # nnlm with nnlm.mse initial values
  nnlm.mse.init <- nnlm(X, as.matrix(y), check.x=F)$coefficients[, 1]
  l[[4]] <- nnlm(X, as.matrix(y), loss="mkl", check.x=F, 
                 init=as.matrix(nnlm.mse.init))$coefficients[, 1]
  obj[[4]] <- Objective4(l[[4]], X, y)
  # pick the best out of the 4 methods; we use this method since
  # each method has its own characteristic numerical instabilities
  selected <- l[[which.min(unlist(obj))]]
  if (any(y>0 & X %*% matrix(selected)==0)) stop("ZERO PREDICTION!!!!\n")
  selected
}

# non-negative regression with regularization
FitNL <- function(X, y, init, lambda, pmat, type="pcs"){
  LB <- rep(0, length(init))
  UB <- rep(Inf, length(init))
  res <- tryCatch({
    nloptr(x0=init, eval_f=Objective3, 
           eval_grad=Gradient3, lb=LB, 
           ub=UB, opts=list("algorithm"="NLOPT_LD_LBFGS", 
                            "xtol_rel"=1e-04), 
           X=X, y=y, lambda=lambda, pmat=pmat)
  }, error=function(e) {
    t.file <- tempfile(tmpdir=getwd())
    cat(t.file, "\n")
    save(init, LB, UB, X, y, lambda, pmat, 
         file=t.file)

    stop("Numerical error in BFGS\n")
  })
  if (type=="scores"){
    if (any(res$solution==0) & any(y>0 & (X %*% matrix(res$solution))==0)){
      nnlm.res <- nnlm(X,
                     as.matrix(y), 
                     loss="mkl",
                     check.x=F)$coefficients[, 1]
      res$solution <- nnlm.res
    }
  }
  res$solution
}