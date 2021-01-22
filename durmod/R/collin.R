# The following is adapted from
# http://www.eco.rug.nl/gauss/gauss99/msg00585.html
# The idea is to cholesky (or qr) decompose the
# hessian (or some other matrix)
# name the columns with (close to) zero pivot B
# and the others A, and solve the system AX = B
# This yields weights, i.e. how the variables corresponding
# to the columns of B may be written as linear combinations
# of the other variables.
collin <- function(mat,eps=NULL,weps=0.01) {
  N <- dim(mat)[1]
  nm <- colnames(mat)
  if(is.null(nm)) nm <- paste('var',1:N)
  ch <- try(chol(mat,pivot=TRUE), silent=TRUE)
  if(inherits(ch,'try-error')) stop("Can't check for multicollinearities")
  pivot <- attr(ch,'pivot')

  if(is.null(eps)) {
    rank <- attr(ch,'rank')
    okcol <- 1:rank
  } else {
    rank <- min(which(c(diag(ch),-Inf) <= eps*c(sqrt(diag(mat)),1)))-1
    okcol=1:rank
  }

  if(rank == N) return(NULL)

#  print(c(rank,N,diag(ch)[rank]))
  
  a <- ch[okcol,okcol]
  b <- ch[okcol,-okcol]
  rm(ch)
  gc()

  cf <- as.matrix(solve(a,b))
  rm(a,b)
  gc()
  badnum <- N-rank
  badnames <- nm[pivot[-okcol]]
  weightnames <- nm[pivot[okcol]]
  # Now cf is a rank x badnum matrix which describes the collinearities in the data
  # We actually want a list. Which weights should we include, i.e. which columns is it
  # a linear combination of?
  # the easy way out is to set a threshold above which we include them, but we really want
  # the most dominant ones
  # we sort them ascending in absolute value into val. Find the largest jump
  # everything after is included
  res <- list()
  for(i in 1:ncol(cf)) {
    val <- sort(abs(cf[,i]))
#    thres <- val[which.max(diff(log(val)))+1]
#    wvars <- which(abs(cf[,i]) >= thres)
    wvars <- which(abs(cf[,i]) >= weps)
    wnames <- weightnames[wvars]
    frm <- data.frame(weight=cf[wvars,i],row.names=wnames)
    res[[i]] <- frm
  }
  names(res) <- badnames
  gc()
  res
}

