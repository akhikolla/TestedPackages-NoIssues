#' @title Compressive Orthogonal Matching Pursuit with Replacement
#' @description An implementation of the Compressive Orthogonal Matching Pursuit with Replacement algorithm
#' 
#' @details COMPR is an iterative greedy method, which alternates between expanding 
#' the cluster centroid set \eqn{C} with a new element \eqn{c_i}, whose sketch is the most correlated to the residue and 
#' the global minimization with respect to cluster centroids \eqn{c_1, \dots, c_K} and their weights \eqn{w_1, \dots, w_K}.   
#' It clusters the data collection into K groups 
#' by minimizing the difference between the compressed data version (data sketch) and 
#' a linear combination of cluster centroid sketches, \emph{i.e.} \eqn{\|Sk(Data) - \sum_{i=1}^K w_i \cdot Sk(c_i)\|}.

#' @param Data A Filebacked Big Matrix n x N, data vectors are stored in the matrix columns.
#' @param ind.col Column indeces, which indicate which data vectors are considered for clustering. By default the entire \code{Data} matrix.   
#' @param K Number of clusters.
#' @param Frequencies A frequency matrix m x n with frequency vectors in rows. 
#' @param lower_b A vector of the lower boundary of data. 
#' @param upper_b A vector of the upper boundary.
#' @param SK_Data Data sketch vector of the length 2m. It can be computed using \code{\link{Sketch}}. 
#' @param maxIter Maximum number of iterations in the global optimization with respect to cluster centroid vectors and their weights. Default is 300.
#' @param HardThreshold logical that indicates whether to perform the replacement. Default is TRUE.
#' @param options List of optimization parameters: 
#'              \itemize{
#'               \item \code{tol_centroid} is a tolerance value for the centroid optimization. Default is 1e-8.
#'               \item \code{nIterCentroid} is a maximum number of iterations in the centroid optimization (default is 1500). 
#'                \item \code{min_weight} is a lower bound for centroid weights (default is 0).
#'                \item \code{max_weight} is an upper bound for centroids weights (default is Inf)
#'                \item \code{nIterLast} is a number of iteration in the global optimization at the last algorithm iteration. Default is 1000.
#'                \item \code{tol_global} is a tolerance value for the global optimization. Default is 1e-12.
#'              }
#' @return A matrix n x K with cluster centroid vectors in columns.
#' @examples 
#' X = matrix(rnorm(1e5), ncol=1000, nrow = 100)
#' lb = apply(X, 1, min)
#' ub = apply(X, 1, max)
#' X_FBM = bigstatsr::FBM(init = X, ncol=1000, nrow = 100)
#' out = GenerateFrequencies(Data = X_FBM, m = 20, N0 = ncol(X_FBM))
#' SK = Sketch(Data = X_FBM, W = out$W)
#' C  <- COMPR(Data = X_FBM, K = 2, Frequencies = out$W, lower_b = lb, upper_b = ub, SK_Data = SK)
#' @note This method is also referred to as Compressive K-means and it has been published in 
#' \insertRef{DBLP:journals/corr/KerivenTTG16}{chickn}.
#' @importFrom nloptr nloptr
#' @importFrom stats optim
#' @importFrom pracma lsqnonneg
#' @export
COMPR <-function(Data, 
                 ind.col = 1:ncol(Data), 
                 K, 
                 Frequencies, 
                 lower_b, upper_b, 
                 SK_Data, 
                 maxIter = 300, 
                 HardThreshold = TRUE,
                 options = list('tol_centroid' = 1e-8,
                                'nIterCentroid' = 1500, 
                                'min_weight' = 0, 'max_weight' = Inf,
                                'nIterLast' = 1000, 'tol_global' = 1e-12)){
    
    if(nrow(Data) != ncol(Frequencies)){
      stop('Error: incompatible data and frequency matrix dimensions')
    }
    if(length(SK_Data) != 2*nrow(Frequencies)){
      stop(paste0('Error: incompatible sketch length ', length(SK_Data), 
                  ' and frequency matrix dimensions', nrow(Frequencies)))
    }
  #=========== Initialization ================================================================================================
    m = nrow(Frequencies)
    d = nrow(Data)
    N = length(ind.col)
    r = SK_Data

    if(HardThreshold){
      iters = c(1:(2*K))
      C = matrix(NA, ncol = K+1, nrow = d) # K+1 column is reserved to hard thresholding
      SK_C = matrix(NA, nrow = 2*m, ncol = K+1)
      prodWC = matrix(NA, nrow = m, ncol = K+1)
      norm_SK_C = matrix(NA, ncol = K+1, nrow = 1)
      ObjF = rep(NA, 2*K)
      weight = rep(NA, K+1)
    }else{
      iters = 1:K
      C = matrix(NA, ncol = K, nrow = d) # K column 
      SK_C = matrix(NA, ncol = K, nrow = 2*m)
      prodWC = matrix(NA, nrow = m, ncol = K)
      ObjF = rep(NA, K)
      weight = rep(NA, K)
    }
    
    for( t in iters){
      dim_C_t = min(t, K+1)
      #=== New centroid ======================================== 
      res_C <- optim(par = Data[,sample(ind.col, 1)], fn = ObjFun_OMP_cpp, gr = Gradient_OMP_cpp, 
                     method = 'L-BFGS-B', lower = lower_b, upper = upper_b, W = Frequencies, 
                     residue = r, control = list("maxit" = options$nIterCentroid, 'trace' = 0))
      C[,dim_C_t] = res_C$par
     
      #=== compute sketch of C =================================
      prodWC[,dim_C_t] = Frequencies%*%C[,dim_C_t] # Need to be recomputed at each iteration!!! at the end
      SK_C[,dim_C_t] = c(cos(prodWC[,dim_C_t]), sin(prodWC[,dim_C_t]))

      # ====== HardThreshold ===================================
      if(HardThreshold){
        norm_SK_C[dim_C_t] = sqrt(sum(SK_C[,dim_C_t]*SK_C[,dim_C_t]))
      }
      if(t > K & HardThreshold){
        res_weight_ht = lsqnonneg(as.matrix(SK_C), as.vector(SK_Data))
        beta = res_weight_ht$x*norm_SK_C
        beta_sort = sort.int(beta, index.return = TRUE, decreasing = TRUE)
        if(max(beta_sort$ix[1:K]) > K){
          ind = beta_sort$ix[1:K]
          C[,1:K] = C[, ind]
          SK_C[,1:K] = SK_C[,ind]
        }
      }
      #=== compute weights ===
      res_weight = lsqnonneg(as.matrix(SK_C[,1:min(t, K)]), as.vector(SK_Data))
      weight[1:min(t, K)] = res_weight$x
      #======= Global optimization ========
      x0 = c(as.vector(C[,1:min(t, K)]), weight[1:min(t, K)])
      lb_global = c(rep(lower_b, min(t, K)), rep(options$min_weight, min(t, K)))
      ub_global = c(rep(upper_b, min(t, K)), rep(options$max_weight, min(t, K)))
      if(any(x0 > ub_global)){
        index_error = which((x0 -ub_global)>0)
        epsilon = abs((x0 - ub_global)[index_error])
        if(any(epsilon < 1e-6)){
          x0[index_error] = x0[index_error]-2*epsilon
        }else{
          stop('at least one element in x0 > ub ')
        }
      }
      if(any(x0 < lb_global)){
        index_error = which((lb_global - x0) > 0)
        epsilon = abs((lb_global- x0)[index_error])
        if(any(epsilon < 1e-6)){
          x0[index_error] = x0[index_error]+2*epsilon
        }else{
          stop('at least one element in x0 < ub ')
        }
      }
      res_global = nloptr(x0 = x0, eval_f = Gradient_cpp, W = Frequencies, SK = SK_Data, K= min(t, K), 
                          opts = list('algorithm' = 'NLOPT_LD_LBFGS', 'maxeval' = maxIter, 'print_level' =0),
                          lb = lb_global, ub = ub_global)
      x = res_global$solution
      C[, 1:min(t, K)] = matrix(x[1:(min(t,K)*ncol(Frequencies))], nrow = ncol(Frequencies))
      weight[1:min(t, K)] = x[(min(t,K)*ncol(Frequencies)+1):length(x)]
      
      #===recompute sketch centroids=====
      prodWC[,1:min(t, K)] = Frequencies%*%C[,1:min(t, K)] # Need to be recomputed at each iteration at the end, see (*)!!! 
      SK_C[1:m, 1:min(t, K)] = cos(prodWC[,1:min(t, K)])
      SK_C[(m+1):(2*m), 1:min(t, K)] = sin(prodWC[,1:min(t, K)])
      

      r = SK_Data - as.matrix(SK_C[,1:min(t, K)])%*%weight[1:min(t, K)]
      #cat(t(r)%*%r, '\n')
    }
    
    x0 = c(as.vector(C[,1:min(t, K)]), weight[1:min(t, K)])
    lb_global = c(rep(lower_b, min(t, K)), rep(options$min_weight, min(t, K)))
    ub_global = c(rep(upper_b, min(t, K)), rep(options$max_weight, min(t, K)))
    if(any(x0 > ub_global)){
      index_error = which((x0 -ub_global)>0)
      epsilon = abs((x0 -ub_global)[index_error])
      if(any(epsilon < 1e-6)){
        x0[index_error] = x0[index_error]-2*epsilon
      }else{
        stop('at least one element in x0 > ub ')
      }
    }
    if(any(x0 < lb_global)){
      index_error = which((lb_global - x0) > 0)
      epsilon = abs((lb_global- x0)[index_error])
      if(any(epsilon < 1e-6)){
        x0[index_error] = x0[index_error]+2*epsilon
      }else{
        stop('at least one element in x0 < ub ')
      }
    }
   # cat('Global optimization last iteration...')
    res_global = nloptr(x0 = x0, eval_f = Gradient_cpp, W = Frequencies, SK = SK_Data, K= min(t, K), 
                        opts = list('algorithm' = 'NLOPT_LD_LBFGS', 'maxeval' = options$nIterLast, 'print_level' =0),
                        lb = lb_global, ub = ub_global)

    x = res_global$solution
    C = matrix(x[1:(min(t,K)*d)], nrow = d)

    return(list('C' = C, 'weight' = weight))
}

