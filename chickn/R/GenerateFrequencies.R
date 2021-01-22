#' @title Frequency vector construction
#' 
#' @description Function performs the data variance estimation and the frequency matrix construction. 
#'
#' @details The data variance is estimated on the \code{N0} data vectors randomly selected from \code{Data} 
#' using \code{\link{EstimSigma}} function. The frequency vectors are sampled using \code{\link{DrawFreq}} function.  
#' 
#' @param Data A Filebacked Big Matrix n x N with data vectors in columns.
#' @param N0 Number of data vectors used for the variance estimation in \code{\link{EstimSigma}}.
#' @param ... Additional arguments passed on to \code{\link{EstimSigma}} and \code{\link{DrawFreq}} functions.
#' @param verbose logical that indicates whether dysplay the process steps.
#' @inheritParams DrawFreq 
#' @inheritParams EstimSigma
#' @return A list with the following attributes: 
#' \itemize{
#'            \item \code{W} is the frequency matrix with m frequency vectors in rows. 
#'            \item \code{sigma} is the estimated data variance.
#'            }
#' @examples 
#' X = matrix(rnorm(1000), ncol=100, nrow = 10)
#' X_FBM = bigstatsr::FBM(init = X, ncol=100, nrow = 10)
#' W = GenerateFrequencies(Data = X_FBM, m = 20, N0 = 100, TypeDist = "AR")$W
#' @seealso \code{\link{DrawFreq}}, \code{\link{EstimSigma}}, \code{\link{Sketch}}
#' @references \insertRef{DBLP:journals/corr/KerivenBGP16}{chickn}.    
#' @export
GenerateFrequencies<- function(Data, m, 
                               N0 = 5000,
                               TypeDist = 'AR', verbose = FALSE, ...){
  n = nrow(Data)
  N = ncol(Data)
  if(N0<=0){
    stop('N0 should be a positive number')
  }
  if(N0 > N){
    warning('Number of selected data vectors exceeds matrix dimension')
    N0 = N
  }
  IndColumn = sort(sample(N, N0))
  
  # estimate sigma
  if(verbose){message("Estimate sigma")}
  sigma = EstimSigma(Data = Data, ind.col = IndColumn, m=m, TypeDist = TypeDist, ...)
  # draw frequencies
  if(verbose){message("Draw frequencies")}
  W = DrawFreq(m = m, n = n, sigma = sigma, TypeDist = TypeDist)
  
  return(list('sigma' = sigma, 'W' = W))
}