centroid <- function(lot,  
                dist_method = c("norm1", "norm2", "norm2_square"),
                step_pattern = c("symmetric2", "symmetric1"),
                normalize = TRUE, ws = NULL, ncores = NULL,
                useRcppParallel = TRUE){
   
   dm <- dtw_dismat(lot, dist_method = dist_method,
                    step_pattern = step_pattern,
                    normalize = normalize, ws = ws, 
                    return_matrix = TRUE, ncores = ncores, 
                    useRcppParallel = useRcppParallel)
   centroid_index <- as.numeric(which.min(colSums(dm$dismat))[1])
   return(list(centroid_index = centroid_index,
               dismat_result = dm))
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


DBA <- function(lot, m0 = NULL, iterMax = 10, eps = NULL, 
                dist_method = c("norm1", "norm2", "norm2_square"),
                step_pattern = c("symmetric2", "symmetric1"),
                ws = NULL,
                iter_dist_method = c("dtw_norm1", "dtw_norm2", 
                                     "norm1","norm2", "max", "min"), 
                plotit = FALSE){
   
   .Deprecated("dba")
   return(
      dba(lot = lot, m0 = m0, iterMax = iterMax, eps = eps, dist_method = dist_method,
       step_pattern = step_pattern, ws = ws, iter_dist_method = iter_dist_method, plotit = plotit)
   )
}

dba <- function(lot, m0 = NULL, iterMax = 10, eps = NULL, 
                dist_method = c("norm1", "norm2", "norm2_square"),
                step_pattern = c("symmetric2", "symmetric1"),
                ws = NULL,
                iter_dist_method = c("dtw_norm1", "dtw_norm2", 
                                     "norm1","norm2", "max", "min"), 
                plotit = FALSE){
   # implementation of the paper:
   # "A global averaging method for dynamic time warping, 
   # with applications to clustering"
   
   
   
   dist_method <- match.arg(dist_method)
   step_pattern <- match.arg(step_pattern)
   iter_dist_method <- match.arg(iter_dist_method)
   iter_dist_foo <- match_idm(iter_dist_method)
   ws <- ws_Inf2NULL(ws)
   
   input <- list(iterMax = iterMax, eps = eps, step_pattern = step_pattern,
                 dist_method = dist_method, ws = ws, 
                 iter_dist_method = iter_dist_method)
   
   if(is.null(m0)){
      tmp <- centroid(lot = lot, dist_method = dist_method, 
                      step_pattern = step_pattern, normalize = TRUE, ws = ws)
      m0 <- lot[[tmp$centroid_index]]
   }
   
   # transform vectors to matrices
   if(is.vector(m0)){
      m0 <- as.matrix(m0, ncol=1)
      lot <- lapply(lot, function(xx){as.matrix(xx, ncol=1)})   
   } else{
      plotit <- FALSE
   }
   
   results <- list(input = input,
                   m1 = NA,
                   iterations = list(m0),
                   iterDist_m2lot = list(),
                   iterDist_m2lot_norm = list(),
                   iterDist_m2m = c()
   )
   
   if(plotit) plot(m0, type="l")
   m1 <- m0
   
   for(iter in 1:iterMax){
      resList <- lapply(lot, function(xx){
         IncDTW::dtw(Q = m1, C = xx, ws = ws, 
                     dist_method = dist_method, return_wp = TRUE, step_pattern = step_pattern)
      })
      m2 <- m1
      for(i in 1:nrow(m0)){
         tmp <- lapply(seq_along(lot), function(j){
            ixQ <- resList[[j]]$ii == i
            ixC <- resList[[j]]$jj[ixQ]
            lot[[j]][ixC, , drop = FALSE]
         })
         tmp <- do.call(rbind, tmp)
         tmp <- rbind(tmp, m1[i, , drop = FALSE])
         m2[i, ] <- colMeans(tmp)
      }
      
      iterDist <- iter_dist_foo(m1, m2)
      if(plotit) lines(m2, col="red")
      if(!is.null(eps)){
         if(iterDist < eps) break
      }
      results$iterations[[iter + 1]] <- m2
      results$iterDist_m2lot[[iter]] <- sapply(resList, "[[", "distance")
      results$iterDist_m2lot_norm[[iter]] <- sapply(resList, "[[", "normalized_distance")
      results$iterDist_m2m[iter] <- iterDist
      m1 <- m2
   }#end iter for
   
   results$m1 <- results$iterations[[length(results$iterations)]]
   class(results) <- append(class(results), "dba", 0)
   return(results)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


match_idm <- function(iter_dist_method = c("dtw_norm1", "dtw_norm2", "norm1",
                                           "norm2", "max", "min")){
   # match  iter_dist_method
   iter_dist_method <- match.arg(iter_dist_method)
   
   if(iter_dist_method == "dtw_norm1"){
      return(idm_dtw1)
      
   }else if(iter_dist_method == "dtw_norm2"){
      return(idm_dtw2)
      
   }else if(iter_dist_method == "norm1"){
      return(idm_norm1)
      
   }else if(iter_dist_method == "norm2"){
      return(idm_norm2)
      
   }else if(iter_dist_method == "max"){
      return(idm_max)
      
   }else if(iter_dist_method == "min"){
      return(idm_min)
      
   }else {
      stop("This should not happen. wrong iter_dist_method")
   }
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


idm_dtw1 <- function(m1, m2){
   dtw2vec(m1, m2, dist_method = "norm1", step_pattern = "symmetric2")$normalized_distance
}

idm_dtw2 <- function(m1, m2){
   dtw2vec(m1, m2, dist_method = "norm2", step_pattern = "symmetric2")$normalized_distance
}

idm_norm1 <- function(m1, m2){
   sum(abs(m1-m2))/(ncol(m1) * 2 * nrow(m1))
}

idm_norm2 <- function(m1, m2){
   sqrt(sum((m1-m2)^2))/(ncol(m1) * 2 * nrow(m1))
}

idm_max <- function(m1, m2){
   max(abs(m1-m2))
}

idm_min <- function(m1, m2){
   min(abs(m1-m2))
}

