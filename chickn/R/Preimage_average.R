#' @title Preimage
#' @description Consensus chromatogram computation.
#' @param Data A Filebacked Big Matrix n x N. Data signals are stored in the matrix columns. 
#' @param K_W1 A Filebacked Big Matrix of the Nystrom kernel matrix \eqn{s \times N},
#'  where N is the number of signal in the training collection and s is the Nystrom sample size. 
#' @param C_out A Filebacked Big Matrix of cluster centroids.
#' @param Cl_assign A Filebacked Big matrix of the cluster assignment.
#' @param max_Nsize Maxumum number of neighbors.
#' @param ncores Number of cores.
#' @param DIR_out A directory to save the result, by default it is the working directory.
#' @param BIG logical parameter that controls whether the resulting consensus chromatograms are stored as a Filebacked Big Matrix ('Centroid_preimage.bk').  
#' Default is FALSE.    
#' @return A matrix or a Filebacked Big Matrix if \code{BIG} = TRUE.
#' @export
Preimage <- function(Data,K_W1, C_out, Cl_assign, max_Nsize = 32, ncores = 4, DIR_out = getwd(), 
                     BIG = FALSE){
  if(BIG){
    bf_preimage =  file.path(DIR_out, 'Centroid_preimage')
    if(file.exists(paste0(bf_preimage, '.bk'))){
      file.remove(paste0(bf_preimage, '.bk'))
    }
  }

  #========
  clusters = Cl_assign[,ncol(Cl_assign)]
  centroid_IDs = sort(unique(clusters))
  Centroids = C_out[,centroid_IDs]
  Centroids_norm = apply(Centroids, 2, function(x) x/sqrt(sum(x*x)))
  
  if(BIG){
    C_preimage = bigstatsr::FBM(nrow =nrow(Data), ncol = ncol(Centroids_norm), backingfile = bf_preimage)$save()
  }else{
    C_preimage <- matrix(0, nrow = nrow(Data), ncol = ncol(Centroids_norm))
  }
  
  if(BIG){
    for(i in 1:ncol(Centroids_norm)){
      index = which(clusters == centroid_IDs[i])
     # cat('Cluster', i, 'from', ncol(Centroids_norm), 'length', length(index), '\n')
      
      kernel_c<- bigstatsr::big_apply(K_W1, a.FUN = function(X, ind, Y){
        return(abs(Y%*%X[,ind]))
      } , a.combine = 'c', ind = index, Y = t(Centroids_norm[,i]))
      if(length(index)> 1){
        if(length(index)<= max_Nsize){
          neigbours_IDs = index
        }else{
          neigbours_IDs = index[sort.int(kernel_c, index.return = TRUE, decreasing = TRUE)$ix[1:max_Nsize]]
        }
        
        cluster_sum <- bigstatsr::big_apply(Data, a.FUN = function(X, ind){
          return( rowSums(X[,ind]))
          }, a.combine = 'plus', ind = neigbours_IDs)
        C_preimage[, i]<- cluster_sum/ min(length(index), max_Nsize)# min if cluster size smaller than max_Nsize
      }else{
        C_preimage[,i]<- Data[,index]
      }
    }
  }else{
    for(i in 1:ncol(Centroids_norm)){
      #cat('Centroid', i, '\n')
      index = which(clusters == centroid_IDs[i])
      kernel_c<- bigstatsr::big_apply(K_W1, a.FUN = function(X, ind, Y){
        abs(Y%*%X[,ind])
      } , a.combine = 'c', ind = index, Y = t(Centroids_norm[,i]))
      if(length(index)> 1){
        neigbours_IDs = sort(index[sort.int(kernel_c, index.return = TRUE, decreasing = TRUE)$ix[1:min(length(index), max_Nsize)]])
        C_preimage[, i]<- rowSums(Data[,neigbours_IDs])/min(length(index), max_Nsize)
      }else{
        C_preimage[,i]<- Data[,index]
      }
    }
  }
  return(C_preimage)
}
