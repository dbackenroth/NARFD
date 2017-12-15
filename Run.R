source("Data.R")
source("Alternating.R")

# fits a model and plots estimated FPCs
Test <- function(type="NARFD"){
  l <- Simulate(N=20, D=80, type=type, seed=1)
  long <- MakeLong(l$observations)
  if (type=="NARFD"){
    penalty.seq <- c(0.1, 1, 10)
    # I recommend starting with PenaltySequence(min=10^-3, max=10^4, num=20)
  } else {
    penalty.seq <- c(1, 10^5, 10^8)
    # I recommend starting with PenaltySequence(min=1, max=10^10, num=20)
  }
  reg.res <- RegPath(long=long, npc=2,
                     D=80,
                     type=type,
                     nbasis=20,
                     periodic=F, 
                     seed=1, 
                     folds=3,   # 3 is a bit low, it's common to use 5
                     iter=25, 
                     penalty.seq=penalty.seq)
  # plot estimated FPCs
  library(ggplot2)
  p <- ggplot(reg.res$pcs, aes(x=.index, y=Value, col=FPC)) + 
    geom_line() + 
    theme_bw() + 
    xlab("Time index")
  print(p)
}