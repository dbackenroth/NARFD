library(splines)
library(pbs)
library(refund)
library(NNLM)   # provides mse.mkl function
library(dplyr)
library(tidyr)
library(data.table)
library(glmnet)
library(mgcv)
library(tibble)   # provides rownames_to_column
source("nloptr.R")
options(stringsAsFactors=F)

# inputs:
#   long--data.frame with columns .index (time value), .id (curve label), 
#         .observed (observation of curve at that time value)
#   npc--how many FPCs/functional prototypes to estimate
#   nbasis--how many spline basis functions to use
#   type--NARFD of GFPCA
#   D--how many knots to use in spline basis
#   seed--integer, seed to use for random initialization
#   folds--integer, how many folds to use in cross-validation
#   iter--integer, maximum number of iterations to be used in alternating minimization
#   periodic--T if the spline basis used should be periodic, F otherwise
#   penalty.seq--numeric vector of penalty values to try

# returns a list with elements
#   pred: data.frame with columns 
#         .index--time value
#         .id--curve label
#         .observed--count data
#         .pred--predicted mean (for GFPCA, after exponentiation)
#          FPC1, FPC2, ...--contributions to predicted mean from each FPC
#             (for GFPCA, prior to exponentiation)  [CHECK THIS]
#   scores: data.frame with columns FPC, .id and Value 
#           (Value is the estimated score for that FPC and curve)
#   pcs: data.frame with columns FPCs, .index and Value
#           (Value is the value of the FPC at the corresponding .index)
#   penalty: the value of the selected penalty
#   mean: for NARFD, NA, 
#   for GFPCA, a data.frame with columns .index and Value
#      
RegPath <- function(long, npc, nbasis, type, D, 
                    seed, folds, iter, periodic, 
                    penalty.seq){
  subjs <- unique(long$.id)
  n.curves <- length(subjs)
  pred.error <- list()
  
  # get spline basis functions
  if (periodic) bs.func <- pbs else bs.func <- bs
  Theta <- bs.func(1:D, df=nbasis, intercept=T)
  # this is for fast access to the spline basis
  Theta.dt <- GetThetaDt(Theta, unique.points=unique(long$.index))
  # get initial values
  set.seed(seed)
  inits <- GetInits(type=type, long=long, npc=npc, Theta=Theta)
  # set up folds
  set.seed(seed)
  fold.df <- data.table(Curve=subjs, 
                        Which.fold=sample(rep(1:folds, ceiling(n.curves / folds))[1:n.curves], 
                                          size=n.curves, 
                                          replace=F))
  long <- as.data.table(long)
  
  all.errors <- list()
  for (f in 1:folds){
    cat("Doing fold", f, "\n")
    for (i in 1:length(penalty.seq)){
      p <- penalty.seq[i]
      cat("penalty", p, "\n")
      r <- CrossValidate(fold.df=fold.df, f=f, p=p, long=long, 
                         npc=npc, inits=inits, Theta=Theta, 
                         type=type, iter=iter, Theta.dt=Theta.dt)
      this.error <- r$error
      if (i > 1){
        fold.errors <- rbind(fold.errors, this.error)
      } else {
        fold.errors <- this.error
      }
    }
    all.errors[[f]] <- fold.errors
  }
  errors <- bind_rows(all.errors)
  
  error.summ <- group_by(errors, Penalty) %>% summarize(Mean.MKL=mean(mkl)) %>% arrange(Mean.MKL)
  best.penalty <- filter(error.summ, Mean.MKL==min(Mean.MKL))$Penalty
  # once the best penalty is fixed, redo with all the data and 
  # the selected value
  res.l <- AlternatingFit(dat=long, npc=npc, 
                          inits=inits,
                          Theta=Theta, 
                          type=type, 
                          iter=iter,
                          lambda=best.penalty, 
                          decompose=T)
  pred <- res.l$pred %>% as.data.frame()
  pcs <- ProcessPCs(Theta %*% t(res.l$pc.coefs))
  scores <- ProcessScores(res.l$scores)
  l <- ReorderAndScale(pred, scores, pcs, type=type)
  if (type=="GFPCA"){
    mean <- data.frame(.index=1:nrow(Theta), 
                       Value=Theta %*% t(res.l$mean.coefs))
  } else {
    mean <- NA
  }
  # do PCs, do scores, do mean, return pred, return pen.summ and pen
  list(pred=l$pred, scores=l$scores, pcs=l$pcs, 
       penalty=best.penalty, 
       mean=mean)
}

EstimateScores <- function(curr.scores, curr.pc.coefs, curr.mean.coefs, 
                           dat=dat, type=type, 
                           loss, Theta.dt){
  # update scores
  curve.names <- unique(dat$.id)
  new.scores <- matrix(0, nrow=length(curve.names), ncol=nrow(curr.pc.coefs))
  rownames(new.scores) <- curve.names
  setkey(dat, .id)
  for (cn in curve.names){
    n.curr.pcs <- as.matrix(Theta.dt[dat[cn, .index, ], !"gridval"]) %*% 
      t(curr.pc.coefs)
    observations <- dat[cn, .observed]
    if (type=="NARFD"){
      PenMat = diag(1, nrow=ncol(n.curr.pcs))
      init.vals <- as.matrix(curr.scores[cn, ])
      nl.res <- FitNL.noR(X=n.curr.pcs, y=observations, init=init.vals)
      new.scores[cn, ] <- nl.res
    } else if (type=="GFPCA"){
      n.offset <- as.matrix(Theta.dt[dat[cn, .index, ], !"gridval"]) %*% t(curr.mean.coefs)
      npc <- ncol(n.curr.pcs)
      # sometimes glm will converge without initial values
      # but won't with initial values
      poisson.mod <- tryCatch({
        glm(observations ~ n.curr.pcs + 0, 
            family=poisson(link=log), offset=n.offset,
            start=curr.scores[cn, ])
      }, error=function(e) {
        poisson.mod <- glm(observations ~ n.curr.pcs + 0, 
                           family=poisson(link=log), offset=n.offset) 
      })
      # if it doesn't converge, try with random starting values anyway
      if (!poisson.mod$converged) {
        poisson.mod <- glm(observations ~ n.curr.pcs + 0, 
                           family=poisson(link=log), offset=n.offset)
      }
      new.coefs <- coef(poisson.mod)
      new.scores[cn, ] <- new.coefs
    }
  }
  new.scores
}

EstimatePCs <- function(curr.scores, dat, 
                        curr.pc.coefs, curr.mean.coefs, type, 
                        npc, 
                        Theta.dt, lambda,
                        Theta){
  new.mean.coefs <- NA
  curve.names <- rownames(curr.scores)
  Theta.mats.kp <- vector('list', length(curve.names))
  Y.l <- vector('list')
  id.l <- vector('list')
  for (cn in curve.names){
    cs <- curr.scores[cn, , drop=F]
    if (type=="GFPCA"){
      cs <- cbind(matrix(1), cs)
    }
    Theta.mats.kp[[cn]] <- cs %x% 
      as.matrix(Theta.dt[dat[cn, .index,], !"gridval"])
    Y.l[[cn]] <- dat[cn, .observed]
    id.l[[cn]] <- rep(cn, length(Y.l[[cn]]))
  }
  A <- do.call('rbind', Theta.mats.kp)
  Y <- do.call('c', Y.l)
  ids <- data.table(ID=do.call('c', id.l), Index=1:length(Y))
  setkey(ids, ID)
  diff <- diff(diag(nrow(Theta)), diff = 2)
  
  if (type=="GFPCA"){
    init.vals <- as.numeric(t(rbind(curr.mean.coefs, curr.pc.coefs)))
    PenMat = diag(1, nrow=npc + (type=="GFPCA")) %x% 
        (t(Theta) %*% t(diff) %*% diff %*% Theta)
    
    nl.fit <- gam(Y ~ A - 1, 
                   paraPen=list(A=list(PenMat, sp=lambda)), 
                   drop.intercept=T, 
                   family="poisson", start=init.vals)
    nl.res <- coef(nl.fit)
    new.coefs <- matrix(nl.res, nrow=npc+(type=="GFPCA"), byrow=T)
    new.pc.coefs <- new.coefs[-1, , drop=F]
    new.mean.coefs <- new.coefs[1, , drop=F]
  } else {
    PenMat = diag(1, nrow=npc) %x% (t(Theta) %*% t(diff) %*% diff %*% Theta)
    init.vals <- as.matrix(as.numeric(t(curr.pc.coefs)))
    nl.res <- FitNL(X=A, y=Y, init=init.vals, 
                    lambda=lambda, pmat=PenMat)
    new.pc.coefs <- matrix(nl.res, nrow=npc, byrow=T)
    new.mean.coefs <- NULL
  }
  #print(Theta %*% t(new.pc.coefs) %>% ProcessPCs() %>% ggplot(aes(x=.index, y=Value, col=FPC)) + geom_line()+theme_bw() + ggtitle(p))
  list(new.pc.coefs=new.pc.coefs, new.mean.coefs=new.mean.coefs)
}

# predict function
# grid is either common grid (numeric vector) or grid for each individual
# also creates decomposition for each subject
# pc.coefs has the coefficients for each FPC in a different row
# scores has the scores for each curve, the row names are the curve names
# grid is either a numeric vector or a data.frame with columns .id and .index
# with different points for prediction for each subject
Predict <- function(Theta.dt, pc.coefs, 
                    mean.coefs=NULL, scores,
                    grid, type, decompose=F){
  if (is.numeric(grid)){
    predict.grid <- expand.grid(.id=rownames(scores), 
                                .index=grid, 
                                stringsAsFactors=F) %>% as.data.table()
  } else {
    predict.grid <- grid
  }
  setkey(predict.grid, .id)
  subjects <- unique(predict.grid$.id)
  # if we're estimating a mean, add mean.coefs to the fpc coefficients
  # and a column of 1 to the scores
  if (!type=="NARFD" & !is.null(mean.coefs)){
    use.coefs <- rbind(mean.coefs, pc.coefs)
    use.scores <- cbind(rep(1, nrow(scores)), scores)
    cns <- c("MEAN", paste0("FPC", 1:ncol(scores)))
  } else {
    use.coefs <- pc.coefs
    use.scores <- scores
    cns <- paste0("FPC", 1:ncol(scores))
  }
  rownames(use.coefs) <- cns
  colnames(use.scores) <- cns
  for (s in subjects){
    predict.grid[s, `:=`(.pred=as.matrix(Theta.dt[.index, !"gridval"]) %*% 
                           t(use.coefs) %*% 
                           t(use.scores[s, , drop=F]))]
    if (decompose){
      for (cn in cns){
        predict.grid[s, (cn):=as.matrix(Theta.dt[.index, !"gridval"]) %*% 
                     t(use.coefs[cn, , drop=F]) * 
                     use.scores[s, cn]]
      }
    }
  }
  if (type=="GFPCA"){
    predict.grid[, `:=`(.pred=exp(.pred))]
  }
  predict.grid
}

GetThetaDt <- function(Theta, unique.points){
  Theta.dt <- predict(Theta, unique.points) 
  attr(Theta.dt, "class") <- "matrix"
  Theta.dt <- as.data.table(Theta.dt)[, `:=`(gridval=unique.points)]
  setkey(Theta.dt, gridval)
  Theta.dt
}

Objective <- function(scores, pc.coefs, pred, lambda, PenMat, mean.coefs=NULL){
  s.pc.coefs <-  as.matrix(c(mean.coefs, as.numeric(t(pc.coefs))))
  
  obj <- -sum(dpois(x=pred$.observed, lambda=pred$.pred+1e-16, log=T)) + lambda * t(s.pc.coefs) %*% PenMat %*% s.pc.coefs
  obj
}

AlternatingFit <- function(dat, npc=2, 
                           inits,
                           Theta, 
                           iter=20, 
                           type="NARFD", 
                           tolerance=1e-03, 
                           lambda, decompose=F){
  error.iters <- rep(NA, iter+1)
  curr.scores <- inits$scores
  curr.pc.coefs <- inits$coefs.pcs
  curr.mean.coefs <- inits$coefs.mean
  
  curve.names <- intersect(unique(dat$.id), rownames(curr.scores))
  curr.scores <- curr.scores[curve.names, , drop=F]
  ccc <- suppressWarnings(curve.names==as.numeric(curve.names))
  if (any(!is.na(ccc) & ccc)) stop("No numeric curve names allowed")
  
  dat <- filter(dat, !is.na(.observed), .id %in% curve.names) %>% 
    as.data.table()
  Theta.dt <- GetThetaDt(Theta, unique.points=unique(dat$.index))
  
  curr.pred <- Predict(Theta.dt=Theta.dt, pc.coefs=curr.pc.coefs, 
                       mean.coefs=curr.mean.coefs, 
                       scores=curr.scores,
                       grid=dat, type=type, decompose=F)
  diff <- diff(diag(nrow(Theta)), diff = 2)
  PenMat = diag(1, nrow=npc + (type=="GFPCA")) %x% 
      (t(Theta) %*% t(diff) %*% diff %*% Theta)
  
  old.err <- GetStats(observed=curr.pred$.observed, 
             predicted=curr.pred$.pred, print=F)
  error.iters[1] <- old.err
  for (i in 1:iter){
    obj <- Objective(scores=curr.scores, pc.coefs=curr.pc.coefs, 
                     pred=curr.pred, lambda=lambda, PenMat=PenMat, 
                     mean.coefs=curr.mean.coefs)
    if (iter==1) cat("\n")
    cat(i, ":", obj, "..")
    pc.l <- EstimatePCs(curr.scores=curr.scores, dat=dat, 
                        curr.pc.coefs=curr.pc.coefs, 
                        curr.mean.coefs=curr.mean.coefs, 
                        type=type, 
                        lambda=lambda,
                        npc=npc, Theta.dt=Theta.dt,  
                        Theta=Theta)
    curr.pc.coefs <- pc.l$new.pc.coefs
    curr.mean.coefs <- pc.l$new.mean.coefs
    curr.pred <- Predict(Theta.dt=Theta.dt, pc.coefs=curr.pc.coefs, 
                         mean.coefs=curr.mean.coefs, 
                         scores=curr.scores,
                         grid=dat, type=type, decompose=F)
    if (type=="GFPCA"){
      curr.pcs <- Theta %*% t(curr.pc.coefs)
      orth.curr.pcs <- tryCatch({
         svd(curr.scores %*% t(curr.pcs))$v[, 1:npc]
      }, error = function(e) {
         curr.pcs
      })
      proj.matrix <- solve(t(Theta) %*% (Theta)) %*% t(Theta)
      orth.curr.pc.coefs <- proj.matrix %*% orth.curr.pcs
      curr.pc.coefs <- t(orth.curr.pc.coefs)
    }
    ptm <- proc.time()
    curr.scores <- EstimateScores(curr.scores=curr.scores, 
            curr.pc.coefs=curr.pc.coefs, 
            curr.mean.coefs=curr.mean.coefs,
            dat=dat, type=type, 
            Theta.dt=Theta.dt)
   if (type=="GFPCA"){   # center scores for GFPCA
      for (j in 1:npc){
        mean.score <- mean(curr.scores[, j])
        curr.scores[, j] <- curr.scores[, j] - mean.score
        curr.mean.coefs <- curr.mean.coefs + mean.score * curr.pc.coefs[j, ]
      }
    }
    curr.pred <- Predict(Theta.dt=Theta.dt, pc.coefs=curr.pc.coefs, 
                         mean.coefs=curr.mean.coefs, 
                         scores=curr.scores,
                         grid=dat, type=type, decompose=F)
    new.err <- GetStats(observed=curr.pred$.observed, 
                          predicted=curr.pred$.pred, print=F)
    err.rel.change <- abs((old.err - new.err) / old.err)
    error.iters[i+1] = new.err
    if (err.rel.change < tolerance & i > 5) break
    old.err <- new.err
  }
  cat("\n")
  pred <- Predict(Theta.dt=Theta.dt, pc.coefs=curr.pc.coefs, 
                  mean.coefs=curr.mean.coefs, 
                  scores=curr.scores,
                  grid=dat, type=type, decompose=decompose)
  list(scores=curr.scores, pc.coefs=curr.pc.coefs, 
       mean.coefs=curr.mean.coefs, pred=pred, pen.summ=pc.l$pen.summ, 
       penalty=pc.l$penalty, error.iters=error.iters, Theta.dt=Theta.dt)
}


ReorderAndScale <- function(pred, scores, pcs, type){
  norm <- group_by(pcs, FPC) %>%
    summarize(Norm=sqrt(sum(Value^2)))
  pcs <- merge(pcs, norm) %>% mutate(Value=Value / Norm) %>% dplyr::select(-Norm)
  scores <- merge(scores, norm) %>% mutate(Value=Value * Norm) %>% dplyr::select(-Norm)
  fpcs <- unique(scores$FPC)
  fpc.dict <- paste0("FPC", 1:length(fpcs))
  if (type=="NMF"){
    total.contrib <- sort(colSums(abs(pred[, fpcs, drop=F])), decreasing=T)
    names(fpc.dict) <- names(total.contrib)
  } else {
    scores.sd <- group_by(scores, FPC) %>% summarize(Score.SD=sd(Value)) %>% 
      arrange(desc(Score.SD))
    names(fpc.dict) <- scores.sd$FPC
  }
  colnames(pred)[colnames(pred) %in% fpcs] <- fpc.dict[colnames(pred)[colnames(pred) %in% fpcs]]
  scores$FPC <- fpc.dict[scores$FPC]
  pcs$FPC <- fpc.dict[pcs$FPC]
  list(scores=scores, pcs=pcs, pred=pred)
}

GetStats <- function(observed, predicted, print=T){
  objective <- sum(dpois(observed, 
                         predicted, log=T))
  errors <- mse.mkl(observed, predicted)
  if (print){
    cat("MSE.MKL\n")
    print(errors)
  }
  return(errors[2])
}

GetInits <- function(type, long, npc, Theta){
  nb <- ncol(Theta)
  use.subjects <- unique(long$.id)
  if (type=="NARFD"){
    init.scores <- matrix(rnorm(length(use.subjects) * npc), 
                            nrow=length(use.subjects))^2
    row.names(init.scores) <- use.subjects
    inits <- list(coefs.pcs=matrix(rnorm(nb * npc), nrow=npc)^2, 
                    scores=init.scores)
  } else if (type=="GFPCA"){
    init.scores <- matrix(rnorm(length(use.subjects) * npc), 
                            nrow=length(use.subjects))
    row.names(init.scores) <- use.subjects
    inits <- list(coefs.mean=matrix(rnorm(nb), nrow=1), 
                    coefs.pcs=matrix(rnorm(nb * npc), nrow=npc), 
                    scores=init.scores)
  }
  if (type=="GFPCA"){   # center scores for GFPCA
    for (j in 1:npc){
      mean.score <- mean(inits$scores[, j])
      inits$scores[, j] <- inits$scores[, j] - mean.score
      inits$coefs.mean <- inits$coefs.mean + mean.score * inits$coefs.pcs[j, ]
    }
  }
  for (i in 1:ncol(inits$scores)){
    s <- sd(inits$scores[, i])
    inits$scores[, i] <- inits$scores[, i] / s
    inits$coefs.pcs[i, ] <- inits$coefs.pcs[i, ] * s
  }
  inits
}

# cross validation function, returns data frame with columns Penalty, Fold, mkl
CrossValidate <- function(fold.df, f, p, long, npc, inits, Theta, 
                          type, iter, Theta.dt){
  use.subjects <- fold.df[!Which.fold==f, Curve]
  heldout.subjects <- fold.df[Which.fold==f, Curve]
  res.l <- AlternatingFit(dat=long[.id %in% use.subjects, ], 
                          npc=npc, 
                          inits=inits,
                          Theta=Theta, 
                          type=type, 
                          iter=iter,
                          lambda=p, decompose=F)
  heldout.scores <- EstimateScores(
    curr.scores=inits$scores[heldout.subjects, , drop=F], 
    curr.pc.coefs=res.l$pc.coefs, 
    curr.mean.coefs=res.l$mean.coefs, 
    dat=long[.id %in% heldout.subjects, ], 
    type=type,
    Theta.dt=Theta.dt)
  heldout.pred <- Predict(Theta.dt=Theta.dt, pc.coefs=res.l$pc.coefs, 
                          mean.coefs=res.l$mean.coefs,
                          scores=heldout.scores, 
                          grid=long[.id %in% heldout.subjects, ], 
                          type=type)
  prob.sum <- sum(heldout.pred$.pred==0 & heldout.pred$.observed>0)
  list(error=data.frame(Penalty=p, Fold=f, 
                        mkl=mse.mkl(obs=heldout.pred$.observed, 
                                    pred=heldout.pred$.pred)[2]), 
       res.l=res.l)
  
}

ProcessScores <- function(scores){
  as.data.frame(scores) %>% rownames_to_column(var=".id") %>%
    gather(FPC, Value, -.id) %>% mutate(FPC=gsub("V", "FPC", FPC))
}

ProcessPCs <- function(pcs, obs.indices=NULL){
  if (!is.null(obs.indices)) indices=obs.indices else indices=1:nrow(pcs)
  as.data.frame(pcs) %>% 
    mutate(.index=indices) %>% 
    gather(FPC, Value, -.index) %>%
    mutate(FPC=gsub("V", "FPC", FPC))
}
