#' @title Data variance estimation
#' 
#' @description The mean data variance estimation.
#' 
#' @inheritParams Sketch
#' @inheritParams DrawFreq
#' @param nblocks Number of blocks, on which the regression is performed. Default is 32.
#' @param niter Number of iterations. Default is 3. 
#' @param sigma_start An initial value of the data variance. Default is 0.1.
#' @param nparts Number of parts to split the data for the data sketch computation.
#' @param ... Additional arguments passed on to \code{\link{DrawFreq}} function. 
#' @return The estimated data variance.  
#' @examples 
#' X = matrix(rnorm(1e5), ncol=1000, nrow = 100)
#' X_FBM = bigstatsr::FBM(init = X, ncol=1000, nrow = 100)
#' sigma = EstimSigma(Data = X_FBM, ind.col = seq(1,1000, by = 2), m = 20, nblocks = 4)  
#' @seealso \code{\link{DrawFreq}}, \code{\link{Sketch}}, \code{\link{GenerateFrequencies}}
#' @note The idea of the variance estimation on the data fraction is taken from \insertRef{DBLP:journals/corr/KerivenBGP16}{chickn}. 
#' @importFrom pracma lsqnonlin 
#' @importFrom Rdpack reprompt
#' @export
EstimSigma<-function(Data,ind.col, m, nblocks=32, niter=3, sigma_start = 0.1, nparts = 1, ...){
 
  if( nblocks > m){
    nblocks = m%/%2
#    warning('number of blocks is bigger then number of frequencies')
  }
  
  blocksize = round(m/nblocks)
  ind = rbind(seq(1, m, blocksize), seq(1, m, blocksize)+ blocksize-1)
  if(m%% nblocks != 0){
    ind[2,ncol(ind)] = m
  }
  
  sigma = sigma_start;
  N = ncol(Data)
  n = nrow(Data)
  
  for (t in 1:niter){
    Wt = DrawFreq(m = m, n = n, sigma = sigma, ...)
    normWt_sorted = sort(rowSums(Wt*Wt), decreasing = TRUE, index.return = TRUE)
    Wt = Wt[normWt_sorted$ix,]
    Sk = Sketch(Data = Data, ind.col = ind.col, W = Wt, ncores = nparts)#*sqrt(m)
    abs_Sk = sqrt(Sk[1:m]^2 + Sk[(m+1):(2*m)]^2)
    
    Indmax = sapply(1:ncol(ind), function(j)  ind[1,j]-1 + which.max(abs_Sk[ind[1,j]:ind[2,j]]))
    e = abs_Sk[Indmax]
  
    R = normWt_sorted$x[Indmax]
    fun = function(sigma) e - exp(-0.5*R*sigma);
    x0= sigma
    sigma_mean = lsqnonlin(fun,x0);
    sigma = sigma_mean$x
  }
  return(sigma)
}