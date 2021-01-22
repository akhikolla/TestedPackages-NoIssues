idtw <- function(Q, C, newObs, gcm, dm, dist_method = c("norm1", "norm2", "norm2_square"),
                 step_pattern = c("symmetric2", "symmetric1"),
                 diffM = NULL, ws = NULL, 
                 return_cm = FALSE,
                 return_diffM = FALSE,
                 return_wp = FALSE,
                 return_diffp = FALSE,
                 return_QC = FALSE){

   dist_method <- match.arg(dist_method)
   step_pattern <- match.arg(step_pattern)
   if(return_diffp) return_wp <- TRUE
   
   
   if(is.character(C)){
      return_QC <- FALSE
      if(C == "diffM_add"){
         cm_add <- abs(Q)
         n <- nrow(cm_add)#should be equal nrow(gcm)
         m <- ncol(gcm)
         o <- ncol(cm_add)
         
         if(return_diffM) diffM <- Q
         rm(list=c("Q"))
         
      } else if(C == "cm_add"){
         cm_add <- Q
         n <- nrow(cm_add)#should be equal nrow(gcm)
         m <- ncol(gcm)
         o <- ncol(cm_add)
         
         rm(list=c("Q"))
         return_diffp <- FALSE
         diffM <- NA
      } else{
         stop("C needs to be a vector, matrix or one of the two strings: 
              'diffM_add' for difference Matrix or
              'cm_add' for cost Matrix")
      }
      } else {
         if(is.vector(Q)){
            if(dist_method != "norm1"){
               dist_method <- "norm1"
               warning("dist_method is set to 'norm1' for the univariate case")
            }
            n <- length(Q)
            m <- ncol(gcm)
            o <- length(newObs)
         }else{
            if(ncol(Q) != ncol(C) | ncol(C) != ncol(newObs)){
               stop("newObs, C and Q are matrices and must have the same number of columns.")
            }
            return_diffp <- FALSE
            diffM <- NA
            n <- nrow(Q)
            m <- ncol(gcm)
            o <- nrow(newObs)
         }
      }
   
   
   #--- initial checking
   m2 <- o + m
   ws <- initial_ws_check2(nQ = n, nC = m2, ws = ws)
   
   
   #--- preparation
   gcm <- cbind(gcm, matrix(NA, nrow = n, ncol = o))
   dm  <- cbind(dm, matrix(NA, nrow = n, ncol = o))
   if(!is.character(C)){
      ws_cpp <- ifelse(is.null(ws), yes = -1, no = ws)
      if(is.vector(Q)){
         diffM_add <- cpp_diffm(Q, newObs, ws = ws_cpp, nPrevObs = m)
         cm_add <- abs(diffM_add)
      }else{
         cm_add <- cpp_cm(Q, newObs, dist_method = dist_method, ws = ws_cpp, nPrevObs = m)
      }
   }
   
   #--- calculation in C++
   if(is.null(ws)){
      ret <- IGCM_cpp(gcmN = gcm, dmN = dm, cmN = cm_add, step_pattern = step_pattern)
   } else {
      ret <- IGCM_Sakoe_cpp(gcmN = gcm, dmN = dm, cmN = cm_add, ws = ws, step_pattern = step_pattern)
   }
   
   #--- get warping path and diff-path
   if(return_diffp | return_diffM) diffM <- cbind(diffM, diffM_add)
   if(return_wp){
      if(return_diffp){
         if( is.integer(diffM[2,2]) & is.integer(diffM_add[2,1]) ){
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
   ret <- c(list(distance = ret$gcm[n, m2]), ret)
   
   #Normalization
   if(step_pattern == "symmetric2"){
      ret <- c(ret, normalized_distance = 
                  ret$distance/(nrow(gcm) + ncol(gcm) + o) )
   }else{
      ret <- c(ret, normalized_distance = NA )
   }
   
   if(return_cm) ret <- c(ret, list(cm = cm_add))
   if(return_diffM) ret <- c(ret, list(diffM = diffM))
   if(return_QC) ret <- c(ret, list(Q = Q, C = c(C, newObs)))
   class(ret) <- append(class(ret), "idtw", 0)
   return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


initialize_plane <- function(Q, C, dist_method = c("norm1", "norm2", "norm2_square"),
                     step_pattern = c("symmetric2", "symmetric1"), ws = NULL){
   
   dist_method <- match.arg(dist_method)
   step_pattern <- match.arg(step_pattern)
   
   
   if(is.vector(Q)){
      if(dist_method != "norm1"){
         dist_method <- "norm1"
         warning("dist_method is set to 'norm1' for the univariate case")
      }
      nQ <- length(Q)
      nC <- length(C)   
      
   }else if(is.matrix(Q)){
      if(ncol(Q) != ncol(C)){
         stop("C and Q must be vectors or matrices having the same number of columns.")
      }
      nQ <- nrow(Q)
      nC <- nrow(C)   
   }else{
      stop("C and Q must be vectors or matrices having the same number of columns.")
   }

   ret <- idtw2vec(Q = Q,
                   newObs = C, 
                   dist_method = dist_method, 
                   step_pattern = step_pattern, 
                   gcm_lc = NULL, 
                   gcm_lr = NULL, 
                   nC = NULL,
                   ws = ws)
   
   control <- list(dist_method = dist_method, step_pattern = step_pattern,
                   nQ = nQ, nC = nC, ws = ws, reverse = FALSE)
   
   ret$control <- control
   ret$Q <- Q
   ret$C <- C
   # ret$normalized_distance <- ret$distance/(ret$control$nC + ret$control$nQ)
   class(ret) <- append(class(ret), "planedtw", 0)
   return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


reverse <- function(x, ...) {
   UseMethod("reverse")
}

reverse.planedtw <- function(x, ...){
   
   
   if(is.vector(x$Q)){
      ret <- idtw2vec(Q = x$Q[x$control$nQ:1],
                      newObs = x$C[x$control$nC:1], 
                      dist_method = x$control$dist_method, 
                      step_pattern = x$control$step_pattern, 
                      gcm_lc = NULL, 
                      gcm_lr = NULL, 
                      nC = NULL,
                      ws = x$control$ws)
      ret$Q <- x$Q[x$control$nQ:1]
      ret$C <- x$C[x$control$nC:1]
      
   }else if(is.matrix(x$Q)){
      ret <- idtw2vec(Q = x$Q[x$control$nQ:1, ,drop = FALSE],
                      newObs = x$C[x$control$nC:1, ,drop = FALSE], 
                      dist_method = x$control$dist_method, 
                      step_pattern = x$control$step_pattern, 
                      gcm_lc = NULL, 
                      gcm_lr = NULL, 
                      nC = NULL,
                      ws = x$control$ws)
      ret$Q <- x$Q[x$control$nQ:1, ,drop = FALSE]
      ret$C <- x$C[x$control$nC:1, ,drop = FALSE]
      
   }else{
      stop("C and Q must be vectors or matrices having the same number of columns.")
   }
   
   
   ret$control <- x$control
   ret$control$reverse <- as.logical(1 - x$control$reverse)
   # ret$normalized_distance <- ret$distance/(ret$control$nC + ret$control$nQ)
   
   class(ret) <- append(class(ret), "planedtw", 0)
   return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


refresh <- function(x, ...) {
   UseMethod("refresh")
}

refresh.planedtw <- function(x, ...){
   
   ret <- idtw2vec(Q = x$Q,
                   newObs = x$C,
                   dist_method = x$control$dist_method, 
                   step_pattern = x$control$step_pattern, 
                   gcm_lc = NULL, 
                   gcm_lr = NULL, 
                   nC = NULL,
                   ws = x$control$ws)
   ret$Q <- x$Q
   ret$C <- x$C

   ret$control <- x$control
   # ret$normalized_distance <- ret$distance/(ret$control$nC + ret$control$nQ)
   class(ret) <- append(class(ret), "planedtw", 0)
   return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


increment <- function(x, newObs, direction = c("C", "Q"), ...){
   UseMethod("increment")
}


increment.planedtw <- function(x, newObs, direction = c("C", "Q"), ...){
   direction <- match.arg(direction)

   #initial check
   if(is.matrix(newObs)){
      if(is.matrix(x$Q)){
         if(ncol(x$Q) != ncol(newObs)){
            stop("newObs must have same number of columns as Q and C.")
         }
      }else{
         stop("Q and C are vectors but newObs is a matrix.")
      }
   }


   if(direction == "C"){
      #--- increment C
      ret <- idtw2vec(Q = x$Q,
                      newObs = newObs,
                      dist_method = x$control$dist_method,
                      step_pattern = x$control$step_pattern,
                      gcm_lc = x$gcm_lc_new,
                      gcm_lr = x$gcm_lr_new,
                      nC = x$control$nC,
                      ws = x$control$ws)

      if(is.matrix(x$Q)){
         nC_new <- x$control$nC + nrow(newObs)
         ret$C <- rbind(x$C, newObs)
      }else{
         nC_new <- x$control$nC + length(newObs)
         ret$C <- c(x$C, newObs)
      }
      ret$Q <- x$Q
      ret$control <- x$control
      ret$control$nQ <- x$control$nQ
      ret$control$nC <- nC_new

   }else if(direction == "Q"){
      #--- increment Q
      ret <- idtw2vec(Q = x$C,
                      newObs = newObs,
                      dist_method = x$control$dist_method,
                      step_pattern = x$control$step_pattern,
                      gcm_lc = x$gcm_lr_new,
                      gcm_lr = x$gcm_lc_new,
                      nC = x$control$nQ,
                      ws = x$control$ws)

      if(is.matrix(x$Q)){
         nQ_new <- x$control$nQ + nrow(newObs)
         ret$Q <- rbind(x$Q, newObs)
      }else{
         nQ_new <- x$control$nQ + length(newObs)
         ret$Q <- c(x$Q, newObs)
      }
      ret$C <- x$C
      ret$control <- x$control
      ret$control$nC <- x$control$nC
      ret$control$nQ <- nQ_new

      nms <- names(ret)
      nms <- switch_names(nms, "gcm_lr_new", "gcm_lc_new")
      names(ret) <- nms
   }

   # ret$normalized_distance <- ret$distance/(ret$control$nC + ret$control$nQ)
   class(ret) <- append(class(ret), "planedtw", 0)
   return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


decrement <- function(x, direction = c("C", "Q", "both"), refresh_dtw = FALSE, nC = NULL, nQ = NULL, ...) {
   UseMethod("decrement")
}

decrement.planedtw <- function(x, direction = c("C", "Q", "both"), refresh_dtw = FALSE, nC = NULL, nQ = NULL, ...){
   
   # primitive decrement 
   if(!is.null(nC) | !is.null(nQ)){
     
       if(!is.null(nC)){
         if(is.matrix(x$C)){
            x$C <- x$C[1:nC, ,drop = FALSE]
         }else{
            x$C <- x$C[1:nC]
         }
         x$control$nC <- nC
      }
      
      if(!is.null(nQ)){
         if(is.matrix(x$Q)){
            x$C <- x$Q[1:nQ, ,drop = FALSE]
         }else{
            x$C <- x$Q[1:nQ]
         }
         x$control$nQ <- nQ
      }
      
      
      x$distance <- NA
      x$normalized_distance <- NA   
      
      # try to keep the normalized 
      if(!is.null(nC) & is.null(nQ)){
         if(!is.null(x$gcm_lr_new)){
            x$distance <- x$gcm_lr_new[nC]
            x$normalized_distance <- x$distance/(x$control$nC + x$control$nQ)
         }
      } else if (is.null(nC) & !is.null(nQ)){
         if(!is.null(x$gcm_lc_new)){
            x$distance <- x$gcm_lc_new[nQ]
            x$normalized_distance <- x$distance/(x$control$nC + x$control$nQ)
         }
      }
      
      # delete invalid parts
      x$gcm_lr_new <- NULL
      x$gcm_lc_new <- NULL
      
       
   } else{
      
      # decrement with help of dtw_partial
      direction <- match.arg(direction)
      
      if(direction %in% c("both", "C")){
         partial_C <- TRUE
      }else{
         partial_C <- FALSE   
      }
      if(direction %in% c("both", "Q")){
         partial_Q <- TRUE
      }else{
         partial_Q <- FALSE   
      }
      
      if(is.null(x$gcm_lc_new) | is.null(x$gcm_lr_new)){
         return(x)
      }
      
      y <- dtw_partial(x, partial_Q = partial_Q, 
                       partial_C = partial_C, reverse = x$control$reverse)
      
      nQ_new <- y$rangeQ[2] - y$rangeQ[1] + 1
      nC_new <- y$rangeC[2] - y$rangeC[1] + 1
      
      if(nQ_new == x$control$nQ & nC_new == x$control$nC){
         
         # no changes: full alignment is best
         return(x)
         
      }else if(nQ_new != x$control$nQ & nC_new == x$control$nC){
      
         if(is.matrix(x$Q)){
            x$Q <- x$Q[y$rangeQ[1]:y$rangeQ[2], , drop = FALSE]
         }else{
            x$Q <- x$Q[y$rangeQ[1]:y$rangeQ[2]]
         }  
         x$distance <- x$gcm_lc_new[nQ_new]
            
      }else if(nQ_new == x$control$nQ & nC_new != x$control$nC){
         
         if(is.matrix(x$C)){
            x$C <- x$C[y$rangeC[1]:y$rangeC[2], , drop = FALSE]
         }else{
            x$C <- x$C[y$rangeC[1]:y$rangeC[2]]
         }
         x$distance <- x$gcm_lr_new[nC_new]
      }
      
      x$control$nQ <- nQ_new
      x$control$nC <- nC_new
      x$gcm_lr_new <- NULL
      x$gcm_lc_new <- NULL
      x$normalized_distance <- y$normalized_distance
   }
   
   
   if(refresh_dtw){
      return(refresh(x))
   }else{
      return(x)   
   }
   
   # if(refresh_dtw){
   #    return(refresh(x))
   #    ret <- idtw2vec(Q = x$Q,
   #                    newObs = x$C,
   #                    dist_method = x$control$dist_method,
   #                    step_pattern = x$control$step_pattern,
   #                    gcm_lc = NULL,
   #                    gcm_lr = NULL,
   #                    nC = NULL,
   #                    ws = x$control$ws)
   # 
   #    ret$Q <- x$Q
   #    ret$C <- x$C
   #    ret$control <- x$control
   #    ret$control$nQ <- nQ_new
   #    ret$control$nC <- nC_new
   #    # ret$normalized_distance <- ret$distance/(ret$control$nC + ret$control$nQ)
   #    class(ret) <- append(class(ret), "planedtw", 0)
   #    return(ret)
   #    
   # }else{
   #    x$control$nQ <- nQ_new
   #    x$control$nC <- nC_new
   #    x$gcm_lr_new <- NULL
   #    x$gcm_lc_new <- NULL
   #    x$distance <- NA
   #    x$normalized_distance <- y$normalized_distance
   #    return(x)
   # }
   
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


idtw2vec <- function(Q, newObs, dist_method = c("norm1", "norm2", "norm2_square"),
                     step_pattern = c("symmetric2", "symmetric1"),
                     gcm_lc = NULL, gcm_lr = NULL, nC = NULL, ws = NULL){
   
   dist_method <- match.arg(dist_method)
   step_pattern <- match.arg(step_pattern)
   initial_dim_check(Q = Q, C = newObs)
   
   if(is.character(newObs)){
      
      #--- precalculated cost matrix
      if(newObs != "cm"){
         stop("If Q is the costmatrix, C needs to be equal 'cm'.")
      }
      return(idtw2vec_cm(cm = Q, step_pattern = step_pattern,
                         gcm_lc = gcm_lc, gcm_lr = gcm_lr, nC = nC, ws = ws))
      
   }else{
      
      #--- regular case - no precalc cost matrix
      if(is.vector(Q)){
         if(dist_method != "norm1"){
            warning("dist_method is set to 'norm1' for the univariate case")
         }
         return(idtw2vec_univ(Q = Q, newObs = newObs, gcm_lc = gcm_lc,
                              gcm_lr = gcm_lr, nC = nC, ws = ws, 
                              step_pattern = step_pattern))
      }else{
         if(ncol(Q) == 1 & ncol(newObs) == 1){
            if(dist_method != "norm1"){
               warning("dist_method is set to 'norm1' for the univariate case")
            }
            return(idtw2vec_univ(Q = Q, newObs = newObs, gcm_lc = gcm_lc, 
                                 gcm_lr = gcm_lr, nC = nC, ws = ws, 
                                 step_pattern = step_pattern))
         }else{
            return(idtw2vec_multiv(Q = Q, newObs = newObs, dist_method = dist_method,
                                   gcm_lc = gcm_lc, gcm_lr = gcm_lr, 
                                   nC = nC, ws = ws, step_pattern = step_pattern))
         }
      }
   }
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


idtw2vec_cm <- function(cm, step_pattern = c("symmetric2", "symmetric1"),
                        gcm_lc = NULL, gcm_lr = NULL, nC = NULL, ws = NULL){
   
   step_pattern <- match.arg(step_pattern)
   ws <- ws_Inf2NULL(ws)
   
   if(is.null(gcm_lc)){
      # initial case
      if(ncol(cm) <= 1){
         stop("ERROR: in the initial calculation the cost matrix 
              cm needs to have at least 2 columns")
         return(NA)
      }
      
      init_gcm_lc <- cumsum(cm[,1])
      if(is.null(ws)){
         ret <- cpp_dtw2vec_cm_inc(gcm_lc = init_gcm_lc, cm = cm[,-1, drop = FALSE],
                                   step_pattern = step_pattern)
      }else {
         # sakoe chiba warping window
         ret <- cpp_dtw2vec_cm_ws_inc(gcm_lc = init_gcm_lc, cm = cm[,-1, drop = FALSE],
                                      step_pattern = step_pattern, ws = ws, ny = 1)
      }
      ret$gcm_lr_new <- c(tail(init_gcm_lc, 1),  ret$gcm_lr_new)
      
   }else{
      # running case
      if(is.null(ws)){
         ret <- cpp_dtw2vec_cm_inc(gcm_lc = gcm_lc, cm = cm,
                                   step_pattern = step_pattern)
         
      }else if (!is.null(ws) & !is.null(nC)){
         # sakoe chiba warping window
         ret <- cpp_dtw2vec_cm_ws_inc(gcm_lc = gcm_lc, cm = cm,
                                      step_pattern = step_pattern, ws = ws, ny = nC)
         
      }else {
         stop("ERROR: if ws is not NULL, then nC must not be NULL")
         return(NA)
      }   
      
      
      if(!is.null(gcm_lr)){#append former parts of last row of gcm
         ret$gcm_lr_new <- c(gcm_lr,  ret$gcm_lr_new)   
      }
   }
   
   #Normalization
   ret <- c(ret, normalized_distance = NA )
   if(step_pattern == "symmetric2"){
      if(is.null(gcm_lc)){#initial case
         ret$normalized_distance <- ret$distance/(nrow(cm) + ncol(cm)) 
      }else if(!is.null(nC)){
         ret$normalized_distance <-  ret$distance/(nrow(cm) + nC + ncol(cm)) 
      }
   }
return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


idtw2vec_univ <- function(Q, newObs, step_pattern = c("symmetric2", "symmetric1"),
                          gcm_lc = NULL, gcm_lr = NULL, nC = NULL, ws = NULL){
   
   step_pattern <- match.arg(step_pattern)
   initial_dim_check(Q = Q, C = newObs)
   ws <- ws_Inf2NULL(ws)
   
   #--- initial checking
   nQ <- length(Q)
   if(!is.null(ws)){
      if(is.null(gcm_lc)){
         ws <- initial_ws_check2(nQ = nQ, nC = (length(newObs) + 0), ws = ws)
      }else{
         ws <- initial_ws_check2(nQ = nQ, nC = (length(newObs) + nC), ws = ws)
      }
   } 
   
   
   if(is.null(gcm_lc)){
      # initial case
      if(length(newObs) <= 1){
         stop("ERROR: in the initial calculation the time 
              series newObs needs to be at least 2 observations long")
         return(NA)
      }
      
      init_gcm_lc <- cumsum(abs(Q-newObs[1]))
      if(is.null(ws)){
         ret <- cpp_dtw2vec_inc(x=Q, newObs = newObs[-1], 
                                gcm_lc = init_gcm_lc, step_pattern = step_pattern)
         
      }else {
         # sakoe chiba warping window
         if(nQ > ws + 1) init_gcm_lc[(ws + 2):nQ] <- NA
         ret <- cpp_dtw2vec_inc_ws(x=Q, newObs = newObs[-1], 
                                   gcm_lc = init_gcm_lc, 
                                   ws = ws, ny = 1, step_pattern = step_pattern)
      }
      ret$gcm_lr_new <- c(tail(init_gcm_lc, 1),  ret$gcm_lr_new)
      
      
   }else{
      # running case
      if(is.null(ws)){
         ret <- cpp_dtw2vec_inc(x=Q, newObs = newObs, gcm_lc = gcm_lc , step_pattern = step_pattern)
         
      }else if (!is.null(ws) & !is.null(nC)){
         # sakoe chiba warping window
         ret <- cpp_dtw2vec_inc_ws(x=Q, newObs = newObs, gcm_lc = gcm_lc,
                                   ws=ws, ny = nC, step_pattern = step_pattern)
         
      }else {
         stop("ERROR: if ws is not NULL, then nC must not be NULL")
         return(NA)
      }   
      
      
      if(!is.null(gcm_lr)){#append former parts of last row of gcm
         ret$gcm_lr_new <- c(gcm_lr,  ret$gcm_lr_new)   
      }
   }
   
   #Normalization
   ret <- c(ret, normalized_distance = NA )
   if(step_pattern == "symmetric2"){
      if(is.null(gcm_lc)){#initial case
         ret$normalized_distance <- ret$distance/(length(Q) + length(newObs)) 
      }else if(!is.null(nC)){
         ret$normalized_distance <-  ret$distance/(length(Q) + nC + length(newObs)) 
      }
   }
   
   return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


idtw2vec_multiv <- function(Q, newObs, dist_method = c("norm1", "norm2", "norm2_square"),
                            step_pattern = c("symmetric2", "symmetric1"),
                            gcm_lc = NULL, gcm_lr = NULL, nC = NULL, ws = NULL){
   
   step_pattern <- match.arg(step_pattern)
   dist_method <- match.arg(dist_method)
   initial_dim_check(Q = Q, C = newObs)
   ws <- ws_Inf2NULL(ws)
   
   if(!is.null(ws)){
      if(is.null(gcm_lc)){
         ws <- initial_ws_check2(nQ = nrow(Q), nC = (nrow(newObs) + 0), ws = ws)
      }else{
         ws <- initial_ws_check2(nQ = nrow(Q), nC = (nrow(newObs) + nC), ws = ws)
      }
   }
   
   
   if(is.null(gcm_lc)){
      # initial case
      if(length(newObs) <= 1){
         stop("ERROR: in the initial calculation the time 
              series newObs needs to be at least 2 observations long")
         return(NA)
      }
      
      # TODO: also possible to solve with is.NUll in C++ code
      # https://stackoverflow.com/questions/34205925/default-null-parameter-rcpp/34206602
      
      if(is.null(ws)){
         init_gcm_lc <- cpp_cm(Q, newObs[1, , drop = FALSE], dist_method = dist_method, ws = -1, nPrevObs = 0)
      }else{
         init_gcm_lc <- cpp_cm(Q, newObs[1, , drop = FALSE], dist_method = dist_method, ws = ws, nPrevObs = 0)
      }
      init_gcm_lc <- cumsum(init_gcm_lc)   
      
      if(is.null(ws)){
         ret <- cpp_dtw2vec_inc_mv(x=Q, newObs = newObs[-1, , drop = FALSE], 
                                   gcm_lc = init_gcm_lc, dist_method = dist_method, step_pattern = step_pattern)
         
      }else {
         # sakoe chiba warping window
         ret <- cpp_dtw2vec_inc_mv_ws(x=Q, newObs = newObs[-1, , drop = FALSE], 
                                      gcm_lc = init_gcm_lc, dist_method = dist_method,  
                                      ws = ws, ny = 1, step_pattern = step_pattern)
      }
      ret$gcm_lr_new <- c(tail(init_gcm_lc, 1),  ret$gcm_lr_new)
      
      
   }else{
      # running case
      if(is.null(ws)){
         ret <- cpp_dtw2vec_inc_mv(x=Q, newObs = newObs, gcm_lc = gcm_lc,
                                   dist_method = dist_method, step_pattern = step_pattern)
         
      }else if (!is.null(ws) & !is.null(nC)){
         # sakoe chiba warping window
         ret <- cpp_dtw2vec_inc_mv_ws(x=Q, newObs = newObs, gcm_lc = gcm_lc,
                                      dist_method = dist_method, ws=ws, ny = nC, step_pattern = step_pattern)
         
      }else {
         stop("ERROR: if ws is not NULL, then nC must not be NULL")
         return(NA)
      }   
      
      
      if(!is.null(gcm_lr)){#append former parts of last row of gcm
         ret$gcm_lr_new <- c(gcm_lr,  ret$gcm_lr_new)   
      }
   }
   
   #Normalization
   ret <- c(ret, normalized_distance = NA )
   if(step_pattern == "symmetric2"){
      if(is.null(gcm_lc)){#initial case
         ret$normalized_distance <- ret$distance/(nrow(Q) + nrow(newObs)) 
      }else if(!is.null(nC)){
         ret$normalized_distance <-  ret$distance/(nrow(Q) + nC + nrow(newObs)) 
      }
   }
   
   return(ret)
}




#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


switch_names <- function(nms, name1, name2){
   nms <- gsub(x = nms, pattern = name1,replacement = "tmp")
   nms <- gsub(x = nms, pattern = name2, replacement = name1)
   nms <- gsub(x = nms, pattern = "tmp", replacement = name2)
   return(nms)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
