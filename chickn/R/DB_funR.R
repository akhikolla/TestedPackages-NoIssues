#'@title DB_funR
#'
#'@description Davies Bouldin index computation.
#'
#' @param Data_cumsum A Filebacked Big Matrix with signal cumultive sums in columns
#' @param Cluster_partition A list of the resulting cluster assignment. Each list element correspond to a cluster with a list of data vector indices.
#' @param C A matrix of cluster centroids.
#' @param DIR_output An output directory for temporal files.
#' @param ncores Number of cores. Default is 4.
#' @return A positive scalar.  
#' @export
#' @keywords internal
DB_funR <- function(Data_cumsum, Cluster_partition, C, DIR_output, ncores = 4){
  
  K = length(Cluster_partition)
  if(file.exists(file.path(DIR_output, 'dist_intra.bk'))){
    file.remove(file.path(DIR_output, 'dist_intra.bk'))
  }
  
  dist_intr <- bigstatsr::FBM(ncol = K, nrow = 1, backingfile = file.path(DIR_output, 'dist_intra.bk'))$save()
  
  if(file.exists(file.path(DIR_output, 'C_cumsum.bk'))){
    file.remove(file.path(DIR_output, 'C_cumsum.bk'))
  }
  C_cumsum <- bigstatsr::FBM(ncol = ncol(C), nrow = nrow(C), backingfile = file.path(DIR_output, 'C_cumsum'))$save()
  cumsum_parallel(C, C_cumsum)
  
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
  blocks <- split_index(K, 500)
  system.time(for(i in 1:2){
    
    clInd <-  seq(blocks[i, "lower"], blocks[i, "upper"])
    
    system.time(tmp<- foreach::foreach(i = clInd)%dopar% {
      return(IntraDist(Data_cumsum, bigstatsr::rows_along(Data_cumsum), Cluster_partition[[i]], C_cumsum[,i]))
    })
    
  })
  parallel::stopCluster(core)
  
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
  system.time(dist_intra<- foreach(i = 1:500, .combine = 'c')%dopar% {
    return(IntraDist(Data_cumsum, bigstatsr::rows_along(Data_cumsum), Cluster_partition[[i]], C_cumsum[,i]))
  })
  
  system.time(dist_intra<- foreach::foreach(i = 501:1000, .combine = 'c')%dopar% {
    return(IntraDist(Data_cumsum, bigstatsr::rows_along(Data_cumsum), Cluster_partition[[i]], C_cumsum[,i]))
  })
  parallel::stopCluster(core)
  
  if(file.exists(file.path(DIR_output, 'C_cumsum.bk'))){
    file.remove(file.path(DIR_output, 'C_cumsum.bk'))
  }
}