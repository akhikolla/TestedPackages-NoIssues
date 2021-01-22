Fclust.compare <- function(VC, U, index, tnorm = c("minimum","product"))
{

  if (missing(VC))
    stop("The hard partitions VC must be given")
  if (missing(U))
    stop("The fuzzy partitions U must be given")

  if (is.null(VC))
    stop("The hard partitions VC is empty")

  if (is.null(U))
    stop("The fuzzy partitions U is empty")
 VC <- as.numeric(VC)
  U <- as.matrix(U)
  partHard=matrix(0,nrow=length(VC),ncol=length(unique(VC)))
  for (i in 1:length(VC))
  {
    partHard[i, VC[i]]=1
  }


  # if(any(dim(partHard) != dim(partFuzzy)))
  #   stop("partHard and partFuzzy must be matrix with the same dimension")

  lg=length(c("ARI.F","RI.F","JACCARD.F","ALL"))

  if (missing(index))
  {
    indexN=lg
  }
  else
  {
    if (is.null(index))
    {
      indexN=lg
    }
    else
    {
      indexN=match(toupper(index),c("ARI.F","RI.F","JACCARD.F","ALL"))
    }
  }
  if (any(is.na(indexN)))
  {
    indexN=lg
    cat("(At least one) no match is found: all the indices will be computed ",fill=TRUE)
  }
  index=indexN
  out.index=rep(0,lg-1)

  tnorm <- match.arg(tnorm, choices = eval(formals(Fclust.compare)$tnorm))

  U = as.matrix(U)


  out = partition_comp(HardClust = partHard,Fuzzy = U, t_norm = tnorm)

  if (any(index==1) || any(index==lg))
    out.index[1]=out$adjRand.F
  if (any(index==2) || any(index==lg))
    out.index[2]=out$Rand.F
  if (any(index==3) || any(index==lg))
    out.index[3]=out$Jaccard.F

  names(out.index)=c("ARI.F","RI.F","JACCARD.F")
  if (max(index)<lg)
    out.index=out.index[index]
  return(out.index)

}
