
#' @title Nystrom kernel approximation
#' 
#' @description An implementation of the Nystrom kernel approximation method.
#' @param Data A Filebacked Big Matrix n x N. Data vectors are stored in the matrix columns.
#' @param c Number of columns selected for the approximation. 
#' @param l An intermediate rank l < c. 
#' @param s A target rank s < l.
#' @param gamma Kernel parameter. If it is NULL (default), the parameter is estimated using \code{\link{gamma_estimation}}.
#' @param max_neighbors Number of neigbors selected for the paramenter estimation.
#' @param DIR_output A directory for intermediate computations.
#' @param DIR_save A directory to save the result. 
#' @param ncores Number of cores. Default is 2.
#' @param ncores_svd Number of cores used for the SVD computaion. It is recommended to use 1 core (default). 
#' @param distance_type Distance function type. The available types are Wasserstein-1 ('W1') and Euclidean ('Euclide'). 
#' The default value is 'W1'.  
#' @param kernel_type Kernel function type c('Gaussian', 'Laplacian').
#' @param verbose logical that indicates whether dysplay the processing steps.  
#' @return A list with the following attributes: \itemize{
#'            \item \code{K_W1} is the Filebacked Big Matrix of the Nystrom kernel approximation. 
#'            \item \code{gamma} is the estimated kernel parameter.
#'            \item \code{RandomSample} is the data vector indices, selected for the Nystrom approximation. }
#' @details Nystrom method consists in approximating the kernel matrix \eqn{K} by \eqn{ C W^{-1} C^{\top}}, with
#'  \eqn{C \in R^{N \times c}} obtained from \eqn{K} by randomly selecting only \code{c} columns and 
#'  \eqn{W \in R^{c \times c}} obtained from \eqn{C} by selecting as well \code{c} corresponding rows. 
#'  The kernel function, based on the distance metric, is given as follows: \eqn{k(x_i,x_j) = e^{- gamma \cdot d^p(x_i,x_j)}}, 
#'  where \eqn{p} is equal to 1 for 'Laplacian' kernel and equal to 2 for 'Gaussian' kernel and 
#'  where \eqn{d(x_i,x_j)} is the distance between data vectors \eqn{x_i} and \eqn{x_j}.
#'@note This is an implemetation of the Nystrom kernel approximation method proposed in 
#'\insertRef{wang2019scalable}{chickn}.
#' @examples 
#' \donttest{
#' X = matrix(rnorm(2000), ncol=100, nrow = 20)
#' X_FBM = bigstatsr::FBM(init = X, ncol=100, nrow = 20)
#' 
#' output = Nystrom_kernel(Data = X_FBM, c = 10, l = 7, s = 5, 
#'                         max_neighbors = 3, ncores = 2)
#'                         }
#' @seealso \code{\link{W1_parallel}}, \code{\link{gamma_estimation}}, \code{\link[bigstatsr]{big_randomSVD}}, \code{\link{cumsum_parallel}}.
#' @export
Nystrom_kernel <- function(Data,c, l, s, gamma = NULL, 
                            max_neighbors = 32, DIR_output = tempfile(),DIR_save = tempfile(),
                            ncores = 2, ncores_svd = 1, distance_type = 'W1',
                           kernel_type = 'Gaussian', verbose = FALSE){
  
  RcppParallel::setThreadOptions(numThreads = ncores)
  N = ncol(Data)
  n = nrow(Data)
  # ====== Construct Nystrom sample
  set_c = sample(x = 1:N, size = c, replace = FALSE)
  set_c = sort(set_c)
  #cat(set_c, '\n')
  #save(set_c, file = file.path(DIR_output, 'set_c.RData'))
  # ====== Compute_Distance_matrix =============================================================================================
  
  switch(distance_type,
         'W1'={
           bf_cumsum = file.path(DIR_output, 'A_cumsum')
           # if(file.exists(paste0(bf_cumsum, '.bk'))){
           #   cat('Read matrix of cumulative function... ')
           #   A_cumsum = bigstatsr::big_attach(paste0(bf_cumsum, '.rds'))
           # }else{
             A_cumsum = bigstatsr::FBM(nrow = n, ncol = N, backingfile = bf_cumsum)$save()
             if(verbose){message('Compute cumulative function')}
             cumsum_parallel(Data, A_cumsum)
           #}
          
           bf_C <- file.path(DIR_output, 'C')
           # if(file.exists(paste0(bf_C, '.bk'))){
           #   file.remove(paste0(bf_C, '.bk'))
           #   file.remove(paste0(bf_C, '.rds'))
           # }
           C = bigstatsr::FBM( nrow = N, ncol = c, backingfile = bf_C)$save()   # C is n x c
           if(verbose){message('Compute Wasserstein distances')}
           W1_parallel(A_cumsum, C, set_c)
           
           # if(file.exists(paste0(bf_cumsum, '.bk'))){
           #   file.remove(paste0(bf_cumsum, '.bk'))
           #   file.remove(paste0(bf_cumsum, '.rds'))
           # }
         },
         'Euclide'={
           bf_C <- file.path(DIR_output, 'C')
           # if(file.exists(paste0(bf_C, '.bk'))){
           #   file.remove(paste0(bf_C, '.bk'))
           # }
           if(verbose){message('Compute Euclidean distances')}
           C = bigstatsr::FBM( nrow = N, ncol = c, backingfile = bf_C)$save()   # C is n x c
           E_parallel(Data, C, set_c)
         },
         stop('incorrect distance type. Choose between: W1 and E'))

  #======== Compute kernel ============================================================
  if(verbose){message("Compute kernel")}
  if(is.null(gamma)){
    gamma = gamma_estimation(C, max_neighbors, kernel_type = kernel_type)
  }

  switch(kernel_type,
         'Gaussian'={
           bigstatsr::big_apply(C, a.FUN = function(X, ind, gamma){
             X[,ind] <- exp(-gamma*X[,ind]*X[,ind])
             NULL
           }, a.combine = 'cbind', ind = bigstatsr::cols_along(C), gamma = gamma, ncores = ncores)
         },
         'Laplacian'={
           bigstatsr::big_apply(C, a.FUN = function(X, ind, gamma){
             X[,ind] <- exp(-gamma*X[,ind])
             NULL
           }, a.combine = 'cbind', ind = bigstatsr::cols_along(C), gamma = gamma, ncores = ncores)
         },
         'none'={
          message('Kernel is not computed. The distance matrix is returned.')
         }, stop("Unkown kernel type"))
  #====== compute SVD ==========================================================================
  if(verbose){message('Perform SVD')}
  bf_W = file.path(DIR_output, 'W')
   # if(file.exists(paste0(bf_W, '.bk'))){
   #   file.remove(paste0(bf_W, '.bk'))
   # }
  W = bigstatsr::big_copy(X = C, ind.row = set_c, backingfile = bf_W)$save()

  SVD_W = bigstatsr::big_randomSVD(X = W, k = l, tol = 1e-4, ncores = ncores_svd)
  
  # file.remove(paste0(bf_W, '.bk'))
  # file.remove(paste0(bf_W, '.rds'))
  
  UD_sqrt_inv = SVD_W$u %*% diag(1/sqrt(SVD_W$d))  # matrix dimension c x l

  R = bigstatsr::big_prodMat(C, UD_sqrt_inv)   # ATTENTION big_prodMat output simple matrix N x l 
  
  # file.remove(paste0(bf_C, '.bk'))
  # file.remove(paste0(bf_C, '.rds'))
  
  bf_R = file.path(DIR_output, 'R')
  # if(file.exists(paste0(bf_R, '.bk'))){
  #   file.remove(paste0(bf_R, '.bk'))
  # }
  R = bigstatsr::FBM(nrow = N, ncol = l, init = R, backingfile = bf_R)$save()
  
  SVD_R = bigstatsr::big_randomSVD(X = R, k = s, tol = 1e-4, ncores = ncores_svd)    # R_s = U_s Sigma_s  t(V_s)

  B = bigstatsr::big_prodMat(R, SVD_R$v)
  # file.remove(paste0(bf_R, '.bk'))
  # file.remove(paste0(bf_R, '.rds'))
  bf_reduced = file.path(DIR_output, 'B')
  # if(file.exists(paste0(bf_reduced, '.bk'))){
  #   file.remove(paste0(bf_reduced, '.bk'))
  # }
  B = bigstatsr::FBM(nrow = N, ncol = s, backingfile = bf_reduced, init = B)$save()

  # if(file.exists(file.path(DIR_output,'NystromKernel.bk'))){
  #   file.remove(file.path(DIR_output,'NystromKernel.bk'))
  # }
  K_W1 = bigstatsr::big_transpose(X = B, backingfile = file.path(DIR_save,'NystromKernel'))$save()
  
  # file.remove(paste0(bf_reduced, '.bk'))
  # file.remove(paste0(bf_reduced, '.rds'))
  
  bigstatsr::big_apply(K_W1, a.FUN = function(X, ind){
    X[,ind]<- X[,ind]/matrix(sqrt(colSums(X[,ind]^2)), nrow = nrow(X), ncol = ncol(X[,ind]), byrow = TRUE)
    NULL
  }, a.combine = 'c', ind = bigstatsr::cols_along(K_W1), ncores = ncores)
  
  retList <- list("K_W1" = K_W1, 'gamma' = gamma, 'RandomSample' = set_c)
  return(retList)
}
