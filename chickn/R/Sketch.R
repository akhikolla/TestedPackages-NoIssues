#' @title split_index
#' @description splits the indeces into several blocks
#' @param total_len is a total length of the index vector
#' @param block_len is a length of the block
#' @param nb_split is a number of the splits
#' @keywords internal
#' @return A list of the splits that contains lower, upper block indeces and block size  
split_index <- function(total_len, block_len,
                      nb_split = ceiling(total_len / block_len)) {
  
  #assert_one_int(total_len); assert_pos(total_len)
  
  if (nb_split > total_len) {
    nb_split <- total_len
  } else if (nb_split == 0) {  ## block_len = Inf
    nb_split <- 1
  }
  #assert_one_int(nb_split); assert_pos(nb_split)
  
  int <- total_len / nb_split
  upper <- round(1:nb_split * int)
  lower <- c(1, upper[-nb_split] + 1)
  size <- c(upper[1], diff(upper))
  
  cbind(lower, upper, size)
}
#'@title split_foreach
#' @description implements block foreach loop
#' @param FUN is a function that is executed in the loop
#' @param ind is a index vector
#' @param ... are parameters of the function FUN
#' @param .combine is a rule to combine the result
#' @param ncores is a nuber of cores
#' @param nb_split is a nuber of splits
#' @return The result of the FUN function
#' @importFrom foreach %dopar% foreach
#' @keywords internal
#' @export
split_foreach <- function(FUN, ind, ...,
                           .combine = NULL,
                           ncores =4,
                           nb_split = ncores) {
  switch (Sys.info()['sysname'],
          'Linux' = {
            cluster_Type = 'FORK'
          },
          'Windows' = {
            cluster_Type = 'PSOCK'
          },
          "Darwin" = {
            cluster_Type = 'FORK'
          },
          "SunOS" = {
            cluster_Type = 'PSOCK'
          },
          stop(paste("Package is not compatible with", Sys.info()['sysname']))
  )
  cl <- parallel::makeCluster(ncores, type = cluster_Type)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  intervals <- split_index(length(ind), nb_split = ncores)
  ic <- NULL
  res <- foreach(ic = bigstatsr::rows_along(intervals)) %dopar% {
    ind.part <- ind[seq(intervals[ic, "lower"], intervals[ic, "upper"])]
    FUN(ind = ind.part, ...)
  }
  
  `if`(is.null(.combine), res, do.call(.combine, res))
}
#'@title big_parallelize
#' @description parallel call a function on file-backed matrix
#' 
#' @param X is a file-backed data matrix
#' @param p.FUN is a function to apply
#' @param p.combine is a rule to combine the results from each core
#' @param ind is a vector of column indeces on which the computations are performed
#' @param ncores is a nuber of cores
#' @param ... are additional parameters
#' @return The result of the parallel computations
#' @keywords internal
#' @export
big_parallelize <- function(X, p.FUN,
                            p.combine = NULL,
                            ind = bigstatsr::cols_along(X),
                            ncores = 4,
                            ...) {
  
  split_foreach(
    p.FUN, ind, X, ...,
    .combine = p.combine,
    ncores = ncores)
}

#' @title Sketch
#' 
#' @description The data sketch computation.
#' @param Data A Filebacked Big Matrix n x N. Data signals are stored in the matrix columns.
#' @param W  A frequency matrix m x n. The frequency vectors are stored in the matrix rows.  
#' @param ind.col Column indeces for which the data sketch is computed. By default all matrix columns.  
#' @param ncores Number of used cores. By default 1. If \code{parallel} = FALSE, ncores defines a number of data splits
#'  on which the sketch is computed separatelly.
#' @param parallel logical parameter that indicates whether computations are performed on several cores in parallel or not. 
#' @return The data sketch vector. 
#' @details The sketch of the given data collection \eqn{x_1, \dots, x_N} is a vector of the length 2m. 
#' First m components of the data sketch vector correspond to its real part, \emph{i.e.} \eqn{\frac{1}{N} \sum_{i=1}^N \cos(W x_i)}.
#' Last m components are its imaginary part, \emph{i.e.} \eqn{\frac{1}{N} \sum_{i=1}^N \sin(W x_i)}.
#' @examples 
#' X = matrix(rnorm(1000), ncol=100, nrow = 10)
#' X_FBM = bigstatsr::FBM(init = X, ncol=100, nrow = 10)
#' W = GenerateFrequencies(Data = X_FBM, m = 20, N0 = 100, TypeDist = "AR")$W
#' SK1 = Sketch(X_FBM, W)
#' SK2 = Sketch(X_FBM, W, parallel = TRUE, ncores = 2)
#' all.equal(SK1, SK2)
#' @references \insertRef{DBLP:journals/corr/KerivenBGP16}{chickn}.    
#' @export
Sketch<-function(Data,W,ind.col = 1:ncol(Data), ncores = 1, parallel = FALSE){
  
  m = nrow(W)
  N = length(ind.col)
  switch(class(Data),
         'FBM' = {
           if(!parallel){
             weight = rep(1/N, ncol(Data))
             # SK = bigstatsr::big_apply(X = Data, a.FUN = function(X, ind, W, weight){
             #   return(exp(1i*W%*%X[,ind])%*%weight[ind])
             # }, a.combine = 'cbind',ind = ind.col,weight = weight, W = W, ncores = ncores)
             SK = bigstatsr::big_apply(X = Data, a.FUN = function(X, ind, W, weight){
               Z <- W%*%X[,ind]
               return(c(cos(Z)%*%weight[ind], sin(Z)%*%weight[ind]))
             }, a.combine = 'plus', ind = ind.col, weight = weight, W = W, ncores = ncores)
           }else{
             SK = big_parallelize(X = Data, p.FUN = function(X, ind, W, N){
               Z <- W%*%X[,ind]
               return(c(rowSums(cos(Z))/N, rowSums(sin(Z))/N))
             }, p.combine = 'cbind',ind = ind.col, W = W, N= N, ncores = ncores)
             SK = rowSums(SK[])
           }
           # if (ncores >1)
           #   Sk = matrix(rowSums(Sk), ncol = 1)
         },
         "matrix" = {
           weight = rep(1/N, N)
           SK = exp(1i*W%*%Data[])%*%weight ### QQQQ -1i or 1i????? R: it is +1i
         }, 
         "numeric" = {
           SK = exp(1i*W%*%Data[])
         }, 
         "integer" = {
           SK = exp(1i*W%*%Data[])
         }, 
         stop("Unkown class of Data"))
  
  return(SK) #/sqrt(m)
}

