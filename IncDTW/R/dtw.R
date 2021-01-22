
cm <- function(Q, C, dist_method = c("norm1", "norm2", "norm2_square"), ws = NULL, ...){
   
   if(!is.matrix(Q)){ Q <- as.matrix(Q) }
   if(!is.matrix(C)){ C <- as.matrix(C) }
   
   if(class(dist_method) == "function"){
      
      nq <- nrow(Q); nc <- nrow(C)
      mat <- matrix(0, nrow = nq, ncol = nc)
      for(i in 1:nq){
         for(j in 1:nc){
            mat[i,j] <- dist_method(Q[i, , drop = FALSE], C[j, , drop = FALSE], ...)
         }
      }
      return(mat)   
      
   }else{
   
      dist_method <- match.arg(dist_method)
      ws_cpp <- ifelse(is.null(ws), yes = -1, no = ws)
      return(cpp_cm(Q, C, dist_method = dist_method, ws = ws_cpp, nPrevObs = 0))   
   }
   
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


dtw <- function(Q, C, dist_method = c("norm1", "norm2", "norm2_square"), 
                step_pattern = c("symmetric2", "symmetric1"), ws = NULL, 
                return_cm = FALSE,
                return_diffM = FALSE,
                return_wp = FALSE,
                return_diffp = FALSE,
                return_QC = FALSE){
   
   dist_method <- match.arg(dist_method)
   step_pattern <- match.arg(step_pattern)
   if(return_diffp) return_wp <- TRUE
   
   if(is.character(C)){
      if(C == "diffM"){
         cm <- abs(Q)
         n <- nrow(cm)
         m <- ncol(cm)
         diffM <- Q
         rm(list=c("Q"))
         
      } else if(C == "cm"){
         cm <- Q
         n <- nrow(cm)
         m <- ncol(cm)
         return_diffp <- FALSE
         diffM <- NA
      } else{
         stop("C needs to be a vector, matrix or one of the two strings: 'diffM' for difference Matrix or 'cm' for cost Matrix")
      }
   } else {
      if(is.vector(Q)){
         if(dist_method != "norm1"){
            dist_method <- "norm1"
            warning("dist_method is set to 'norm1' for the univariate case")
         }
         n <- length(Q)
         m <- length(C)   
         
      }else{
         if(ncol(Q) != ncol(C)){
            stop("C and Q are matrices and must have the same number of columns.")
         }
         return_diffp <- FALSE
         diffM <- NA
         n <- nrow(Q)
         m <- nrow(C)   
      }
   }
  
   #--- initial checking
   ws <- initial_ws_check2(nQ = n, nC = m, ws = ws)
   
   #--- preparation
   if(!is.character(C)){
      ws_cpp <- ifelse(is.null(ws), yes = -1, no = ws)
      if(is.vector(Q)){
         # univariate case
         diffM <- cpp_diffm(Q, C, ws = ws_cpp, nPrevObs = 0)
         cm <- abs(diffM)
      }else{
         # multivariate case
         cm <- cpp_cm(Q, C, dist_method = dist_method, ws = ws_cpp, nPrevObs = 0)
      }
   }
   
   #--- calculation in C++
   if(is.null(ws)){
      ret <- GCM_cpp(cm, step_pattern = step_pattern)
   } else {
      ret <- GCM_Sakoe_cpp(cm, ws, step_pattern = step_pattern)
   }
   
   #--- get warping path and diff-path
   if(return_wp){
      if(return_diffp){
         if(is.integer(diffM[2,2])){
            tmp <- BACKTRACK2II_cpp(ret$dm, diffM)   
         }else{
            tmp <- BACKTRACK2IN_cpp(ret$dm, diffM)
         }
         ret <- c(ret, list(diffp = tmp$diffp))
      }else {
         tmp <- BACKTRACK_cpp(ret$dm)
      }
      ret <- c(ret, list(ii = rev(tmp$ii),
                         jj = rev(tmp$jj),
                         wp = rev(tmp$wp)))
   }
   
   #--- add dtw_distance
   ret <- c(ret, list(distance = ret$gcm[n, m]))
   
   #Normalization
   if(step_pattern == "symmetric2"){
      ret <- c(ret, normalized_distance = ret$distance/(n+m) )
   }else{
      ret <- c(ret, normalized_distance = NA )
   }

   if(return_cm) ret <- c(ret, list(cm = cm))
   if(return_diffM) ret <- c(ret, list(diffM = diffM))
   if(return_QC) ret <- c(ret, list(Q = Q, C = C))
   class(ret) <- append(class(ret), "idtw", 0)
   return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


dtw_partial <- function(x, partial_Q = TRUE, partial_C = TRUE, reverse = FALSE){
   
   
   if(!is.null(x$gcm)){
      nQ <- nrow(x$gcm)
      nC <- ncol(x$gcm)
      lc <- rep(Inf, nQ) # last column
      lr <- rep(Inf, nC) # last row
      
      if(partial_Q) lc <- x$gcm[, nC]
      if(partial_C) lr <- x$gcm[nQ, ]
      
   }else if(!is.null(x$gcm_lc_new)){
      nQ <- length(x$gcm_lc_new)
      nC <- length(x$gcm_lr)
      lc <- rep(Inf, nQ) # last column
      lr <- rep(Inf, nC) # last row
      
      if(partial_Q) lc <- x$gcm_lc_new
      if(partial_C) lr <- x$gcm_lr
      
   }else{
      stop("Error, x needs to have either the attribute 'gcm' or 'gcm_lr_new' as a result
           of dtw() or idtw2vec()!")
   }
   lr[is.na(lr)] <- Inf
   lc[is.na(lc)] <- Inf
   
   vals <- list( partialC = lr/(  nQ + 1:nC),# partial matching of C
                 partialQ = lc/(1:nQ +   nC))# partial matching of Q
   mins <- sapply(vals, min)
   ix <- which.min(mins)[1]
   ixb <- which.min(vals[[ix]])[1]
   
   
   if(reverse){# open begin matching
      if(ix == 1){# partial matching of C
         rangeQ <- c(1, nQ)
         rangeC <- c(nC-ixb+1, nC)
         
      } else if(ix == 2){# partial matching of Q
         rangeQ <- c(nQ-ixb+1, nQ)
         rangeC <- c(1, nC)
      }
   }else{# open end matching
      if(ix == 1){# partial matching of C
         rangeQ <- c(1, nQ)
         rangeC <- c(1, ixb)
         
      } else if(ix == 2){# partial matching of Q
         rangeQ <- c(1, ixb)
         rangeC <- c(1, nC)
      }
   }
   
   return(list(rangeQ = rangeQ, rangeC = rangeC, 
               normalized_distance = vals[[ix]][ixb]))
   
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


dec_dm <- function(dm, Ndec, diffM = NULL){

   if(!is.integer(dm[2,2])){
      warning("The direction matrix is no integer matrix -> takes longer to process")
   }
   
  Nnew <- ncol(dm) - Ndec
  if(Nnew == 1){
     warning("Nnew = 1, calculation is maybe meaningless")
     ret <- list(ii = nrow(dm):1,
                 jj = rep(1, nrow(dm)),
                 wp = rep(3, nrow(dm)-1) )
  } else {
     if(is.null(diffM)){
        tmp <- BACKTRACK_cpp(dm[, 1:Nnew])
        ret <- list(ii = rev(tmp$ii),
                    jj = rev(tmp$jj),
                    wp = rev(tmp$wp))   
     }else{
        if(is.integer(diffM[2,2])){
           tmp <- BACKTRACK2II_cpp(dm[, 1:Nnew], diffM)
        }else{
           tmp <- BACKTRACK2IN_cpp(dm[, 1:Nnew], diffM)
        }
        ret <- list(ii = rev(tmp$ii),
                    jj = rev(tmp$jj),
                    wp = rev(tmp$wp),
                    diffp = rev(tmp$diffp))   
     }
     
  }
  return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


dtw2vec <- function(Q, C, dist_method = c("norm1", "norm2", "norm2_square"),
                    step_pattern = c("symmetric2", "symmetric1"),
                    ws = NULL, threshold = NULL){
   
   dist_method <- match.arg(dist_method)
   step_pattern <- match.arg(step_pattern)
   initial_dim_check(Q = Q, C = C)

   
   
   if(is.character(C)){
      if(C != "cm"){
         stop("If Q is the costmatrix, C needs to be equal 'cm'.")
      }
      return(dtw2vec_cm(cm = Q, step_pattern = step_pattern, ws = ws, threshold = threshold))      
   }
   
   
   if(is.vector(Q)){
      if(dist_method != "norm1"){
         warning("dist_method is set to 'norm1' for the univariate case")
      }
      return(dtw2vec_univ(Q = Q, C = C, step_pattern = step_pattern, ws = ws, threshold = threshold))
   }else{
      if(ncol(Q) == 1 & ncol(C) == 1){
         if(dist_method != "norm1"){
            warning("dist_method is set to 'norm1' for the univariate case")
         }
         return(dtw2vec_univ(Q = Q, C = C, step_pattern = step_pattern, 
                             ws = ws, threshold = threshold))
      }else{
         return(dtw2vec_multiv(Q = Q, C = C, dist_method = dist_method, step_pattern = step_pattern, 
                               ws = ws, threshold = threshold))   
      }
   }
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


dtw2vec_cm <- function(cm, step_pattern = c("symmetric2", "symmetric1"), ws = NULL, threshold = NULL){
   
   # vector based implementation with pre-calculated cost matrix
   if(is.null(ws) & is.null(threshold)){
      ret <- cpp_dtw2vec_cm(cm = cm, step_pattern = step_pattern)
      
   }else {
      # ea and sakoe chiba
      ws <- initial_ws_check2(nrow(cm), ncol(cm), ws)
      if(is.null(ws)){ws <- max(c(nrow(cm), ncol(cm)))}
      if(is.null(threshold)) threshold <- Inf
      ret <- cpp_dtw2vec_cm_ws_ea(cm = cm, step_pattern = step_pattern, 
                                  ws = ws, threshold = threshold)
   }
   
   #Normalization
   if(step_pattern == "symmetric2"){
      ret <- list(distance = ret, normalized_distance = ret/(nrow(cm) + ncol(cm)) )
   }else{
      ret <- list(distance = ret, normalized_distance = NA )
   }
   return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


dtw2vec_univ <- function(Q, C, step_pattern = c("symmetric2", "symmetric1"), ws = NULL, threshold = NULL){
   # wrapper function for the C++ implementation
   
   step_pattern <- match.arg(step_pattern)
   initial_dim_check(Q = Q, C = C)
   # params <-  get_params(step_pattern = step_pattern,
   #                       ws = ws, threshold = threshold)
   ws <- initial_ws_check2(nQ = length(Q), nC = length(C), ws = ws)
   
   
   if(is.null(ws) & is.null(threshold)){
      # fastest implementation for unrestricted warp
      # if(step_pattern == "symmetric1"){
      #    ret <- cpp_dtw2vec_v32(Q, C)   
      # }else{
      #    ret <- cpp_dtw2vec(Q, C, step_pattern = step_pattern)
      # }
      ret <- cpp_dtw2vec(Q, C, step_pattern = step_pattern)
      
   }else if (!is.null(ws) & is.null(threshold)){
      # sakoe chiba warping window
      ret <- cpp_dtw2vec_ws(x=Q,y=C, step_pattern = step_pattern, ws = ws)
      
   }else if (is.null(ws) & !is.null(threshold)){
      # early abandoning if the costs exceeds threshold
      ret <- cpp_dtw2vec_ea(x=Q, y=C, step_pattern = step_pattern, threshold = threshold)
      
   }else{
      # ea and sakoe chiba
      ret <- cpp_dtw2vec_ws_ea(x=Q, y=C, step_pattern = step_pattern, ws = ws, threshold = threshold)   
      
   }
   
   #Normalization
   if(step_pattern == "symmetric2"){
      ret <- list(distance = ret, normalized_distance = ret/(length(Q) + length(C)) )
   }else{
      ret <- list(distance = ret, normalized_distance = NA )
   }
   return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


dtw2vec_multiv <- function(Q, C, dist_method = c("norm1", "norm2", "norm2_square"),
                           step_pattern = c("symmetric2", "symmetric1"), ws = NULL, threshold = NULL){
   # wrapper function for the C++ implementation
   
   dist_method <- match.arg(dist_method)
   step_pattern <- match.arg(step_pattern)
   initial_dim_check(Q = Q, C = C)
   ws <- initial_ws_check2(nQ = nrow(Q), nC = nrow(C), ws = ws)
   
   if(is.null(ws) & is.null(threshold)){
      # fastest implementation for unrestricted warp
      ret <- cpp_dtw2vec_mv(Q, C, step_pattern = step_pattern, dist_method = dist_method)
      
   }else {
      # ea and sakoe chiba
      if(is.null(ws)) ws <- max(c(nrow(Q), nrow(C)))
      if(is.null(threshold)) threshold <- Inf
      ret <- cpp_dtw2vec_mv_ws_ea(x=Q, y=C, step_pattern = step_pattern, 
                                  dist_method = dist_method,ws = ws, threshold = threshold)
   }
   
   #Normalization
   if(step_pattern == "symmetric2"){
      ret <- list(distance = ret, normalized_distance = ret/(nrow(Q) + nrow(C)) )
   }else{
      ret <- list(distance = ret, normalized_distance = NA )
   }
   
   return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


dtw_dismat <- function(lot, dist_method = c("norm1", "norm2", "norm2_square"),
                                step_pattern = c("symmetric2", "symmetric1"),
                                normalize = TRUE, ws = NULL, threshold = NULL,
                                return_matrix = FALSE, ncores = NULL, useRcppParallel = TRUE){
   
   dist_method  <- match.arg(dist_method)
   step_pattern <- match.arg(step_pattern)
   
   input <- list(dist_method = dist_method, ws = ws,  normalize = normalize,
                 threshold = threshold, step_pattern = step_pattern,
                 ncores = ncores)
   
   if(normalize == TRUE & step_pattern == "symmetric1"){
      stop("It is recommended to use step_pattern 'symmetric2' to compare normalized 
           DTW distances for time series of differnt lengths. Either set 'normalized' to 'FALSE' 
           or 'step_pattern' to 'symmetric2'.")
   }
   
   
   # initial dim check
   lapply(lot[-1], function(x){initial_dim_check(Q = lot[[1]], C = x)})
   ws <- ws_Inf2NULL(ws)
   
   if(is.vector(lot[[1]])){
      is_vector <- TRUE 
      lapply(lot, function(x){
         lapply(lot, function(y){initial_ws_check(nQ = length(x), nC = length(y), ws = ws)})   
      })
      
   } else{
      is_vector <- FALSE   
      lapply(lot, function(x){
         lapply(lot, function(y){initial_ws_check(nQ = nrow(x), nC = nrow(y), ws = ws)})   
      })
      
   }
   
   
   if(is.vector(lot[[1]])){
      is_vector <- TRUE 
   } else{
      is_vector <- FALSE   
   }
   

   if(is.null(ncores)){
      ncores <- detectCores() - 1 
   } else if (ncores <= 0){
      stop("Error: ncores needs to be >= 1")
   }
   if(ncores > 1){
      RcppParallel::setThreadOptions(numThreads = ncores)   
   }
   
   
   N <- length(lot)
   NN <- (N*(N-1))/2
   jj <- sapply(1:(N-1), function(j){rep(j, N-j)})
   ii <- sapply(1:(N-1), function(i){(i+1):N})
   # regard labels of lot
   tmp <- rep(0,N)
   if(!is.null(names(lot))){
      names(tmp) <- names(lot)
   }
   dis <- dist(tmp)

   
   if(ncores > 1){
      if(useRcppParallel){
         # PARALLEL EXECUTION
         ii <- unlist(ii) - 1 # minus 1 for C compatibility
         jj <- unlist(jj) - 1
         if(is.null(threshold)) threshold <- Inf
         
         if(is_vector){
            # univariate, pass vectors
            if(is.null(ws)) ws <- ws <- max(sapply(lot, length))
            dis[1:NN] <- parallel_dm_dtw(lot = lot, ii, jj, normalize = normalize, 
                                         step_pattern = step_pattern, 
                                         ws = ws, threshold = threshold)
         }else{
            # multivariate, pass matrices 
            if(is.null(ws)) ws <- ws <- max(sapply(lot, nrow))
            dis[1:NN] <- parallel_dm_dtw_mv(lot = lot, ii, jj, normalize = normalize, 
                                            step_pattern = step_pattern, dist_method = dist_method, 
                                            ws = ws, threshold = threshold)
         }
      
         
         
      }else{
         ii <- unlist(ii) 
         jj <- unlist(jj) 
         cl <- makeCluster(ncores)
         clusterExport(cl, envir = environment(), list("lot", "ii", "jj", "dist_method", "ws",
                                                       "threshold", "step_pattern"))
         if(is_vector){
            # univariate, pass vectors
            dis[1:NN] <- parSapply(cl,1:NN, function(k){
               tmp <- IncDTW::dtw2vec_univ(Q = lot[[ ii[k] ]],
                                    C = lot[[ jj[k] ]],
                                    ws = ws, threshold = threshold, step_pattern = step_pattern)
               ifelse(normalize, yes = tmp$normalized_distance, no = tmp$distance)
            })
         }else{
           # multivariate, pass matrices   
            dis[1:NN] <- parSapply(cl,1:NN, function(k){
               tmp <- IncDTW::dtw2vec_multiv(Q = lot[[ ii[k] ]],
                                   C = lot[[ jj[k] ]],
                                   dist_method = dist_method, ws = ws, threshold = threshold,
                                   step_pattern = step_pattern)
               ifelse(normalize, yes = tmp$normalized_distance, no = tmp$distance)
            })
         }
         stopCluster(cl)
      }
      
      
   }else{
      ii <- unlist(ii) 
      jj <- unlist(jj) 
      
      # !!!NO!!! PARALLEL EXECUTION
      if(is_vector){
         # univariate, pass vectors
         dis[1:NN] <- sapply(1:NN, function(k){
            tmp <- IncDTW::dtw2vec_univ(Q = lot[[ ii[k] ]],
                                 C = lot[[ jj[k] ]],
                                 ws = ws, threshold = threshold, step_pattern = step_pattern)
            ifelse(normalize, yes = tmp$normalized_distance, no = tmp$distance)
         })
      }else{
         # multivariate, pass matrices   
         dis[1:NN] <- sapply(1:NN, function(k){
            tmp <- IncDTW::dtw2vec_multiv(Q = lot[[ ii[k] ]],
                                   C = lot[[ jj[k] ]],
                                   dist_method = dist_method, ws = ws, threshold = threshold,
                                   step_pattern = step_pattern)
            ifelse(normalize, yes = tmp$normalized_distance, no = tmp$distance)
         })
      }
   }
   
   
   if(return_matrix){
      dis <- as.matrix(dis)
      return(list(input = input, dismat = dis))
   }else{
      return(list(input = input, dismat = dis))
   }
   
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


dtw_disvec <- function(Q, lot, dist_method = c("norm1", "norm2", "norm2_square"),
                                step_pattern = c("symmetric2", "symmetric1"),
                                normalize = TRUE, ws = NULL, threshold = NULL,
                                ncores = NULL){
   
   dist_method <- match.arg(dist_method)
   step_pattern <- match.arg(step_pattern)
   
   input <- list(dist_method = dist_method, ws = ws,  normalize = normalize,
                 threshold = threshold, step_pattern = step_pattern,
                 ncores = ncores)
   
   if(normalize == TRUE & step_pattern == "symmetric1"){
      stop("It is recommended to use step_pattern 'symmetric2' to compare normalized 
           DTW distances for time series of differnt lengths. Either set 'normalized' to 'FALSE' 
           or 'step_pattern' to 'symmetric2'.")
   }
   
   # initial dim check
   lapply(lot, function(x){initial_dim_check(Q = Q, C = x)})
   ws <- ws_Inf2NULL(ws)
   
   if(is.vector(lot[[1]])){
      is_vector <- TRUE 
      lapply(lot, function(x){initial_ws_check(nQ = length(Q), nC = length(x), ws = ws)})
      
   } else{
      is_vector <- FALSE   
      lapply(lot, function(x){initial_ws_check(nQ = nrow(Q), nC = nrow(x), ws = ws)})
   }
   
   
   if(is.null(ncores)){
      ncores <- detectCores() - 1 
   } else if (ncores <= 0){
      stop("Error: ncores needs to be >= 1")
   }
   if(ncores > 1){
      RcppParallel::setThreadOptions(numThreads = ncores)   
   }
   
   if(ncores > 1){
      if(is.null(threshold)) threshold <- Inf
      # PARALLEL EXECUTION

      if(is_vector){
         # univariate, pass vectors
         if(is.null(ws)) ws <- max(max(sapply(lot, length)), length(Q))
         ret <- parallel_dv_dtw(Q = Q, lot = lot, normalize = normalize, 
                                step_pattern = step_pattern, 
                                ws = ws, threshold = threshold)
         
      }else{
         # multivariate, pass matrices 
         # if(  !( is.null(ws) & is.null(threshold) )  ){
            if(is.null(ws)) ws <- max(max(sapply(lot, nrow)), nrow(Q))
            
         # }  
         ret <- parallel_dv_dtw_mv(Q = Q, lot = lot, normalize = normalize, 
                                   step_pattern = step_pattern, dist_method = dist_method, 
                                   ws = ws, threshold = threshold)
      }
      
      
   }else{
      # !!!NO!!! PARALLEL EXECUTION
      if(is_vector){
         # univariate, pass vectors
         ret <- sapply(1:length(lot), function(k){
            tmp <- IncDTW::dtw2vec_univ(Q = Q,
                                        C = lot[[ k ]],
                                        ws = ws, threshold = threshold, step_pattern = step_pattern)
            ifelse(normalize, yes = tmp$normalized_distance, no = tmp$distance)
         })
      }else{
         # multivariate, pass matrices   
         ret <- sapply(1:length(lot), function(k){
            tmp <- IncDTW::dtw2vec_multiv(Q = Q,
                                          C = lot[[ k ]],
                                          dist_method = dist_method, ws = ws, threshold = threshold,
                                          step_pattern = step_pattern)
            ifelse(normalize, yes = tmp$normalized_distance, no = tmp$distance)
         })
      }
   }
   
   
   #regard labels of lot
   if(!is.null(names(lot))){
      names(ret) <- names(lot)
   }
   
   return(list(input = input, disvec = ret))   
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

norm <- function(x, type = c("z", "01"), 
                 threshold = 1e-5,
                 xmean = NULL, xsd = NULL, xmin = NULL, xmax = NULL){
   
   .Deprecated("scale")
   return(
      scale(x = x, type = type, threshold = threshold, xmean = xmean, xsd = xsd, xmin = xmin, xmax = xmax)
   )
}

scale <- function(x, type = c("z", "01"), threshold = 1e-5,
                 xmean = NULL, xsd = NULL, xmin = NULL, xmax = NULL){
   
   type <- match.arg(type)
   if(type == "z"){
      if(is.vector(x)){
         return(cpp_znorm(x, sd_threshold = threshold, mu_in = xmean, sd_in = xsd))   
      }else{
         return(apply(x, 2, function(y){cpp_znorm(x = y, sd_threshold = threshold, mu_in = xmean, sd_in = xsd)}))   
      }
      
   }else if(type == "01"){
      if(is.vector(x)){
         return(cpp_norm01(x, sd_threshold = threshold, min_in = xmin, max_in = xmax))
      }else{
         return(apply(x, 2, function(y){cpp_norm01(x = y, sd_threshold = threshold, min_in = xmin, max_in = xmax)}))
      }
      
   }else{
      stop("This should not happen. Invalid type parameter")
   }
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


initial_dim_check <- function(Q=Q, C=C){
   is_error <- TRUE
   if(!is.character(C)){
      if(is.vector(Q) & is.vector(C)){
         is_error <- FALSE
         
      }else if(is.matrix(Q) & is.matrix(C)){
         if(ncol(Q) == ncol(C)) is_error <- FALSE
         
      }
   }else{
      if(is.matrix(Q) | is.vector(Q)) is_error <- FALSE
   }
   
   if(is_error){
      stop("Wrong input for Q and C, see the help pages")
   }
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


initial_ws_check <- function(nQ, nC, ws){
   is_error <- TRUE
   if(is.null(ws)){
      is_error <- FALSE

   }else{
      if(abs(nQ - nC) <= ws){
         is_error <- FALSE
      }
   }


   if(is_error){
      stop("window size is too small, no warping path can be found")
   }
}


initial_ws_check2 <- function(nQ, nC, ws){
   
   if(!is.null(ws)){
      if(ws == Inf){
         ws <- NULL
      }else{
         if(abs(nQ - nC) > ws){
            stop("window size is too small, no warping path can be found")
         }   
      }
   }
   return(ws)
}


ws_Inf2NULL <- function(ws){
   
   if(!is.null(ws)){
      if(ws == Inf){
         ws <- NULL
      }
   }
   return(ws)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


