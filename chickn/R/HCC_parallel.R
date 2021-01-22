#' @title hcc_parallel
#' @description Compressed Hierarchical Clustering. 
#' 
#' @param Data A Filebacked Big Matrix n x N. Data signals are stored in the matrix columns. 
#' @param W A frequency matrix m x n with frequency vectros in rows. 
#' @param K Number of clusters at each call of the clustering algorithm.
#' @param maxLevel Maximum number of hierarchical levels.
#' @param ncores Number of cores. By default 4.
#' @param DIR_output An output directory. 
#' @param hybrid logical parameter. If TRUE \code{K} decreases progressively over 
#' hierarchical levels as \eqn{\lceil \frac{K}{level} \rceil}. Default is FALSE.   
#' @param ... Additional arguments passed on to \code{\link{COMPR}}.
#' @param verbose logical that indicates whether dysplay the processing steps.
#' @return The cluster assignment as a list of clusters with corresponding data vector indeces.
#' @details This function provides a divisive hierarchical implementation of \code{\link{COMPR}}. 
#' Parallel computations are performed using 'FORK' clusters (Linux-like platform) or 'PSOCK' clusters (Windows platform) using the \code{parallel} package. 
#' This function generates in the \code{DIR_output} directory the following files:
#' \itemize{
#' \item 'Cluster_assign_out.bk' is a Filebacked Big Matrix N x \code{maxLevel}+1, which stores the cluster assignment at each hierarchical level.
#' \item 'Centroids_out.bk' is a Filebacked Big Matrix with the resulting cluster centroids in columns. 
#'  } 
#'@examples
#'\donttest{
#' data("UPS2")
#' N = ncol(UPS2)
#' n= nrow(UPS2)
#' X_FBM = bigstatsr::FBM(init = UPS2, ncol=N, nrow = n)$save()
#' K_W1 = Nystrom_kernel(Data = X_FBM, c = 14, l = 7, s = 5, 
#'                       max_neighbors = 3, ncores = 1, kernel = 'Gaussian')$K_W1
#' W = GenerateFrequencies(Data = K_W1, m = 20, N0 = ncol(X_FBM))$W
#' C = hcc_parallel(Data = K_W1, W = W, K = 2, maxLevel = 4, 
#'                  DIR_output = tempfile(), ncores = 2)
#'}
#' @seealso \code{\link{COMPR}}
#' @references  \insertRef{DBLP:journals/corr/KerivenTTG16}{chickn}
#' @importFrom utils unstack
#' @importFrom foreach %dopar%      
#' @export
hcc_parallel<- function(Data, W, K, maxLevel, ncores = 2, DIR_output = tempfile(),
                        hybrid = FALSE, verbose = FALSE,...){
  bf_cl_assign =  file.path(DIR_output, 'Cluster_assign_out')
  if(file.exists(paste0(bf_cl_assign, '.bk'))){
    file.remove(paste0(bf_cl_assign, '.bk'))
  }
  Output =bigstatsr::FBM(nrow =ncol(Data), ncol = maxLevel+1, backingfile = bf_cl_assign)$save()
  Output[,1] = 1
  K_global = K
  bf_C_out = file.path(DIR_output, 'Centroids_out')
  if(file.exists(paste0(bf_C_out, '.bk'))){
    file.remove(paste0(bf_C_out, '.bk'))
  }
  if(hybrid){
    size_C_out = prod(sapply(1:(maxLevel+1), function(i) max(floor(K_global/i),2)))
  }else{
    size_C_out = K^maxLevel
  }
  C_out = bigstatsr::FBM(nrow = nrow(Data), ncol = size_C_out, backingfile = bf_C_out)$save()
  
  #===== The first level sketch computation is parallel ========
  if(verbose){message('level 1')}
  SK <- Sketch(Data = Data, ind.col = seq(ncol(Data)), W = W, ncores = ncores, parallel = TRUE)
  #SK = c(Re(SK),Im(SK))
  
  #================= Perform clustering ====================================================
  rowmaxs<- bigstatsr::big_apply(Data, a.FUN = function(X, ind){
    apply(X[, ind],1, max)
  }, a.combine = 'cbind', ind = seq(ncol(Data)))
  rowmins<-bigstatsr::big_apply(Data, a.FUN =function(X, ind){
    apply(X[, ind],1, min)}, a.combine = 'cbind', ind = seq(ncol(Data)))
  rowmaxs = apply(rowmaxs,1, max)
  rowmins = apply(rowmins,1, min)
  CKM_out <- COMPR(Data = Data, ind.col = 1:ncol(Data), K = K, Frequencies = W, 
                   lower_b = rowmins, upper_b = rowmaxs,
                   SK_Data = SK, HardThreshold = TRUE, ...)
  non_null_weights <-which(CKM_out$weight > 0)
  Centroids <-as.matrix(CKM_out$C[,non_null_weights], ncol = length(non_null_weights))
  #================= Compute distance between centroids and data ===========================
  if(length(non_null_weights)>1){
    Centroids_norm = apply(Centroids, 2, function(x) x/sqrt(sum(x*x)))
  }else{
    Centroids_norm = as.matrix(Centroids/sqrt(sum(Centroids*Centroids)), ncol = 1)
  }
  K_c_norm = bigstatsr::big_apply(Data, a.FUN = function(X, ind, Y){
    return(abs(Y%*%X[,ind]))
  } , a.combine = 'cbind', ind = 1:ncol(Data), Y = t(Centroids_norm))#, ncores = nbrc)
  #================ Compute cluster partition ==============================================
  Neighbours = apply(K_c_norm, 2, function(x) which.max(x))
  Output[,2] = Neighbours 
  
  if(maxLevel !=1){
    index = 1:ncol(Data)
    cluster_size = table(Neighbours)
    small_cluster = as.numeric(names(cluster_size)[which(cluster_size <= K)])
    for(s_cl in small_cluster){
      small_cl_index = index[which(Neighbours == s_cl)]
      if(!hybrid){
        Output[small_cl_index, maxLevel + 1] = K^(maxLevel - 1)* (Output[small_cl_index,2]-1) +1
      }else{
        coeff = prod(sapply(3:(maxLevel+1), function(i) max(floor(K_global/i),2)))
        Output[small_cl_index, maxLevel + 1] = coeff* (Output[small_cl_index,2]-1) +1
      }
      c_nbr_in_C_out = unique(Output[small_cl_index, maxLevel+1])
      C_out[, c_nbr_in_C_out] = Centroids[,s_cl]
    }
  }else{
    c_nbr_in_C_out = sort(unique(Output[, level+1]))
    c_nbr = sort(unique(Neighbours))
    C_out[, c_nbr_in_C_out] = Centroids[,c_nbr]
  }
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
  core <- parallel::makeCluster(ncores, type = cluster_Type)
  doParallel::registerDoParallel(core)
  doRNG::registerDoRNG()
  on.exit(parallel::stopCluster(core))
 
  for(level in 2:maxLevel){
      if(hybrid){
        K = max(floor(K_global/level), 2)
      }
      clusters_level = sort(unique(Output[,level]))
      if(verbose){message('level ',level)}
      cl <- NULL
      tmp <- foreach::foreach(cl = clusters_level) %dopar% {
        # for(cl in clusters_level){
        #cat('level ',level,' cluster ', cl, '\n')
        index = which(Output[,level]==cl)
        if(length(index)> K){
          SK <- Sketch(Data = Data, ind.col = index, W = W)
          #================= Perform clustering ====================================================
          rowmaxs<- bigstatsr::big_apply(Data, a.FUN = function(X, ind){
            apply(X[, ind],1, max)
          }, a.combine = 'cbind', ind = index)
          rowmins<- bigstatsr::big_apply(Data, a.FUN =function(X, ind){
            apply(X[, ind],1, min)}, a.combine = 'cbind', ind = index)
          rowmaxs = apply(rowmaxs,1, max)
          rowmins = apply(rowmins,1, min)
          #cat('Clustering\n')
          CKM_out <- COMPR(Data = Data, ind.col = index, K = K, Frequencies = W, 
                           lower_b = rowmins, upper_b = rowmaxs,
                           SK_Data = SK,
                           HardThreshold = TRUE, ...)
          non_null_weights <-which(CKM_out$weight > 0)
          Centroids <-as.matrix(CKM_out$C[,non_null_weights], ncol = length(non_null_weights))
          #================= Compute distance between centroids and data ===========================
          if(length(non_null_weights)>1){
            Centroids_norm = apply(Centroids, 2, function(x) x/sqrt(sum(x*x)))
          }else{
            Centroids_norm = as.matrix(Centroids/sqrt(sum(Centroids*Centroids)), ncol = 1)
          }
          K_c_norm = bigstatsr::big_apply(Data, a.FUN = function(X, ind, Y){
            return(abs(Y%*%X[,ind]))
          } , a.combine = 'cbind', ind = index, Y = t(Centroids_norm))#, ncores = nbrc)
          #================ Compute cluster partition ==============================================
          Neighbours = apply(K_c_norm, 2, function(x) which.max(x))
          Output[index,level+1] = Neighbours + K*(Output[index, level]-1)
          
          if(level < maxLevel){
            cluster_size = table(Neighbours)
            small_cluster = as.numeric(names(cluster_size)[which(cluster_size <= K)])
            #cat('Small clusters\n')
            for(s_cl in small_cluster){
              small_cl_index = index[which(Neighbours == s_cl)]
              if(hybrid){
                coeff = prod(sapply(min((level+1), maxLevel+1):(maxLevel+1), function(i) max(floor(K_global/i),2)))
                Output[small_cl_index, maxLevel + 1] = coeff* (Output[small_cl_index,level+1]-1) +1
              }else{
                Output[small_cl_index, maxLevel + 1] = K^(maxLevel - level)* (Output[small_cl_index,level+1]-1) +1
              }

              c_nbr_in_C_out = unique(Output[small_cl_index, maxLevel+1])
              C_out[, c_nbr_in_C_out] = Centroids[,s_cl]
            }
          }else{
            c_nbr_in_C_out = sort(unique(Output[index, level+1]))
            c_nbr = sort(unique(Neighbours))
            C_out[, c_nbr_in_C_out] = Centroids[,c_nbr]
          }
        }else{
          Output[index,level+1] = 1 + K*(Output[index, level]-1)
        }
        NULL
      }
  }
  clusters = data.frame('row' = 1:ncol(Data), 'col' = Output[,maxLevel+1])
  Cluster_partition = unstack(clusters)
  names(Cluster_partition) = 1:length(Cluster_partition)
  return(Cluster_partition)
}