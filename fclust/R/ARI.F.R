ARI.F <-  function(VC,U, t_norm = c("minimum","product"))
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

  #partHard = as.matrix(partHard)
  U = as.matrix(U)

  # if(any(dim(partHard) != dim(partFuzzy)))
  #   stop("partHard and partFuzzy must be matrix with the same dimension")

  t_norm <- match.arg(t_norm, choices = eval(formals(ARI.F)$t_norm))
  out = partition_comp(HardClust = partHard,Fuzzy = U, t_norm = t_norm)
  return(out$adjRand.F)
}
