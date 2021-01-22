rundtw <- function(Q, C, dist_method = c("norm1", "norm2", "norm2_square"),
                   step_pattern = c("symmetric1", "symmetric2"), k = NULL, 
                   scale = c("01", "z", "none"),  
                   ws = NULL, threshold = NULL, lower_bound = TRUE,
                   overlap_tol = 0, return_QC = FALSE, normalize = c("01", "z", "none")){
   
   
   if (!missing("normalize")){
      warning("Argument 'normalize' is deprecated, use 'scale' instead. 
              The parameter 'scale' is set equal the parameter 'normalize'.")
      scale <- normalize
   }
   
   if(is.list(C)){
      return(
         rundtw_list(Q, C, dist_method, step_pattern, k, 
                     scale, ws, threshold, lower_bound,
                     overlap_tol, return_QC)
      )
   }
   
   debug <- 0#1 # either 0 or 1, if 1 then the cpp function prints computation details, in that case a sink file is recommended
   dist_method <- match.arg(dist_method)
   step_pattern <- match.arg(step_pattern)
   
   if(is.logical(scale)){
      warning("The values TRUE and FALSE for the parameter 'scale' are deprecated. Use '01' or 'none' instead.")
      if(as.logical(scale)){
         scale <- "01" 
      } else{
         scale <- "none"
      }
   }else{
      scale <- match.arg(scale)   
   }
   
   initial_dim_check(Q = Q, C = C)
   
   if(is.null(k)) {
      k <- 0
   }else{
      k <- as.integer(k)
   }
   if(k < 0) k <- 0L
   # build kNN_Info
   kNN_inf_list <- list(vmax = Inf,
                        which_vmax = 0,
                        which_imax = 0,
                        imax = 0,
                        nr_detected = 1,
                        nr_looking4 = k)
   
   if(kNN_inf_list$nr_detected < kNN_inf_list$nr_looking4){
      kNN_inf_list$which_vmax <- -99;   
   }
   
   if(k > 0) lower_bound <- TRUE # k overrules lowerbound parameter
   
   if(is.null(threshold)){
      threshold <- Inf
      use_ea <- 0
   }else{
      if(threshold < 0){
         stop("threshold needs to be >= 0")
      }
      # threshold is set to any values >= 0 and smaller Inf
      use_ea <- 1
   }
   
   
   if(overlap_tol < 0 || overlap_tol >= length(C)/ifelse(is.null(ncol(C)), 1, ncol(C))){
      stop("The overlap tolerance parameter needs to be between 0 and the length of C")
   }
   
   if(lower_bound){
      # I reconsidered this. Since DTW(x,y,sym1) <= DTW(x,y,sym2) in general holds =>
      #                          lb(x,y) <= DTW(x,y,sym1) <= DTW(x,y,sym2) =>
      #                          lb(x,y) <= DTW(x,y,sym2) 
      #                          => use regular lb() also for sym2 step pattern
      #                          it is a valid lower bound, but not so tight
      #                          => comment the following code block
      
      # if(step_pattern != "symmetric1"){
      #    warning("lower_bound is set to FALSE, currently only implemented for step_pattern symmetric1")
      #    use_lb <- 0
      # }else{
      #    use_ea <- 1
      #    use_lb <- 1
      # }
      
      use_ea <- 1
      use_lb <- 1
   }else{
      use_lb <- 0
   }
   
   if(scale == "01"){
      do_scale <- 1
      Q <- IncDTW::scale(Q, type = "01")
      # Q <- norm(Q, type = "01")
      
   }else if(scale == "z"){
      Q <- IncDTW::scale(Q, type = "z")
      # Q <- norm(Q, type = "z")
      
   }else if (scale == "none"){
      do_scale <- 0
   }
   
   if(is.null(ws)) ws <- ifelse(is.vector(Q), yes = length(Q), no = nrow(Q))
   
   
   if((NCOL(Q) == 1 & NCOL(C) == 1)){
      if(dist_method != "norm1"){
         warning("dist_method is set to 'norm1' for the univariate case")
      }
      
      if(scale == "z"){
         ret <- cpp_rundtw_znorm(h = Q, x = C, 
                                 step_pattern = step_pattern, 
                                 ws = ws, threshold = threshold, overlap_tol = overlap_tol, 
                                 use_ea = use_ea, use_lb = use_lb, debug = debug,
                                 kNN_inf_list = kNN_inf_list) 
         
      }else{
         ret <- cpp_rundtw(h = Q, x = C, step_pattern = step_pattern, 
                           ws = ws, threshold = threshold, overlap_tol = overlap_tol,
                           do_norm = do_scale, use_ea = use_ea, use_lb = use_lb, debug = debug,
                           kNN_inf_list = kNN_inf_list) 
      }
      
   }else{
      
      
      if(scale == "z"){
         ret <- cpp_rundtw_znorm_mv(h = Q, x = C, step_pattern = step_pattern, 
                                    dist_method = dist_method, ws = ws, threshold = threshold, overlap_tol = overlap_tol,
                                    use_ea = use_ea, use_lb = use_lb, debug = debug,
                                    kNN_inf_list = kNN_inf_list)         
      }else{
         ret <- cpp_rundtw_mv(h = Q, x = C, step_pattern = step_pattern, 
                              dist_method = dist_method, ws = ws, threshold = threshold, overlap_tol = overlap_tol,
                              do_norm = do_scale, use_ea = use_ea, use_lb = use_lb, debug = debug,
                              kNN_inf_list = kNN_inf_list)   
      }
   }
   
   names(ret$counter) <- c("scale_reset","scale_new_extreme", "scale_1step",
                           "cm_reset", "cm_1step",
                           "early_abandon", "lower_bound")
   ret$counter["completed"] <- (NROW(C) - NROW(Q) + 1) - (ret$counter["early_abandon"] + ret$counter["lower_bound"]) 
   
   if(scale == "z"){
      ret$counter <- ret$counter[c("early_abandon", "lower_bound", "completed")]
   }else if(scale == "none"){
      ret$counter <- ret$counter[c("cm_reset", "cm_1step","early_abandon", "lower_bound", "completed")]
   }
   
   if(k > 0) {
      # if k > nx/nh then too many NN are initiated and returned, drop those here
      all_best_indices <- ret$all_best_indices + 1
      all_best_dist <- ret$dist[all_best_indices]
      new_order <- order(all_best_dist, decreasing = FALSE)
      
      # k_best_indices
      ret$knn_indices <- all_best_indices[new_order[1:k]]
      # k_best_dist    
      ret$knn_values <- all_best_dist[new_order[1:k]]
      
      # finally
      if(debug == 1){
         print("------all_best indices----------")
         print(ret$all_best_indices +1)
         print("------all_best indices end----------")
      }
      ret$all_best_indices <- NULL
      ix_na <- is.na(ret$knn_values)
      ret$knn_values <- ret$knn_values[!ix_na]
      ret$knn_indices <- ret$knn_indices[!ix_na]
      
   }
   
   if(return_QC) ret <- c(ret, list(Q = Q, C = C))
   class(ret) <- append(class(ret), "rundtw", 0)
   return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


rundtw_list <- function(Q, C, dist_method = c("norm1", "norm2", "norm2_square"),
                        step_pattern = c("symmetric1", "symmetric2"), k = NULL, 
                        scale = c("01", "z", "none"), ws = NULL, threshold = NULL, lower_bound = TRUE,
                        overlap_tol = 0, return_QC = FALSE){
   
   
   
   debug <- 0 # either 0 or 1, if 1 then the cpp function prints computation details, in that case a sink file is recommended
   dist_method <- match.arg(dist_method)
   step_pattern <- match.arg(step_pattern)
   
   if(is.vector(Q)){
      if(any(!sapply(C, is.vector))){
         stop("Q is a vector, but not all elements of C are vectors.")
      }
   }else if(is.matrix(Q)){
      nQ <- ncol(Q)
      if(any(sapply(C, ncol) != nQ)){
         stop("Q is a matrix, but not all elements of C are matrices with equal number of columns.")
      }
   }else{
      stop("Q and all elements of C need to be vectors or matrices of equal number of columns.")
   }
   
   if(is.logical(scale)){
      warning("The values TRUE and FALSE for the parameter 'scale' are deprecated. Use '01' or 'none' instead.")
      if(as.logical(scale)){
         scale <- "01" 
      } else{
         scale <- "none"
      }
   }else{
      scale <- match.arg(scale)   
   }
   
   if(is.null(k)) {
      k <- 0
   }else{
      k <- as.integer(k)
   }
   if(k < 0) k <- 0L
   # build kNN_Info
   kNN_inf_list <- list(vmax = Inf,
                        which_vmax = 0,
                        which_imax = 0,
                        imax = 0,
                        nr_detected = 1,
                        nr_looking4 = k)
   
   if(kNN_inf_list$nr_detected < kNN_inf_list$nr_looking4){
      kNN_inf_list$which_vmax <- -99;   
   }
   
   if(k > 0) lower_bound <- TRUE # k overrules lowerbound parameter
   
   if(is.null(threshold)){
      threshold <- Inf
      use_ea <- 0
   }else{
      if(threshold < 0){
         stop("threshold needs to be >= 0")
      }
      # threshold is set to any values >= 0 and smaller Inf
      use_ea <- 1
   }
   
   
   if(overlap_tol < 0 || overlap_tol >= length(C)/ifelse(is.null(ncol(C)), 1, ncol(C))){
      stop("The overlap tolerance parameter needs to be between 0 and the length of C")
   }
   
   if(lower_bound){
      use_ea <- 1
      use_lb <- 1
   }else{
      use_lb <- 0
   }
   
   if(scale == "01"){
      do_scale <- 1
      Q <- IncDTW::scale(Q, type = "01")
      
   }else if(scale == "z"){
      Q <- IncDTW::scale(Q, type = "z")
      
   }else if (scale == "none"){
      do_scale <- 0
   }
   
   if(is.null(ws)) ws <- ifelse(is.vector(Q), yes = length(Q), no = nrow(Q))
   
   ret <- list()
   kup <- list(kNN_val = rep(Inf , k),
               kNN_ix = rep(-99L, k),
               kNN_lot_ix = rep(0L, k),
               kNN_inf_list = kNN_inf_list)
   if(k > 0) {
      
      knn_bestofall <- data.table::data.table(list_id = integer(),
                                              time_ix = integer(),
                                              distance = numeric())
   }
   if(is.vector(Q) || (ncol(Q) == 1 & ncol(C[[1]]) == 1)){
      if(dist_method != "norm1"){
         warning("dist_method is set to 'norm1' for the univariate case")
      }
      
      if(scale == "z"){
         
         for(i in seq_along(C)){
            ret[[i]] <- cpp_rundtw_znorm_lot(h = Q, x = C[[i]], 
                                             kNN_val_in = kup$kNN_val, kNN_ix_in = kup$kNN_ix, kNN_lot_ix_in = kup$kNN_lot_ix, 
                                             kNN_inf_list_in = kup$kNN_inf_list, lot_ix = i,
                                             step_pattern = step_pattern, ws = ws, threshold = threshold, 
                                             overlap_tol = overlap_tol, use_ea = use_ea, use_lb = use_lb, debug = debug) 
            
            if(k > 0){
               kup <- kNN_update(kup = kup, bix = ret[[i]]$all_best_indices, 
                                 bdi = ret[[i]]$dist[ret[[i]]$all_best_indices + 1], 
                                 ki = ret[[i]]$kNN_inf, i = i)
            }
         }
         
      }else{
         
         for(i in seq_along(C)){
            ret[[i]] <- cpp_rundtw_lot(h = Q, x = C[[i]], 
                                       kNN_val_in = kup$kNN_val, kNN_ix_in = kup$kNN_ix, kNN_lot_ix_in = kup$kNN_lot_ix, 
                                       kNN_inf_list_in = kup$kNN_inf_list, lot_ix = i,
                                       step_pattern = step_pattern, ws = ws, threshold = threshold, 
                                       overlap_tol = overlap_tol, do_norm = do_scale, use_ea = use_ea, use_lb = use_lb, debug = debug) 
            
            if(k > 0){
               kup <- kNN_update(kup = kup, bix = ret[[i]]$all_best_indices, 
                                 bdi = ret[[i]]$dist[ret[[i]]$all_best_indices + 1], 
                                 ki = ret[[i]]$kNN_inf, i = i)
            }
         }
      }
      
   }else{
      
      if(scale == "z"){
         
         for(i in seq_along(C)){
            ret[[i]] <- cpp_rundtw_znorm_mv_lot(h = Q, x = C[[i]], 
                                                kNN_val_in = kup$kNN_val, kNN_ix_in = kup$kNN_ix, kNN_lot_ix_in = kup$kNN_lot_ix, 
                                                kNN_inf_list_in = kup$kNN_inf_list, lot_ix = i,
                                                step_pattern = step_pattern, dist_method = dist_method, ws = ws, threshold = threshold, 
                                                overlap_tol = overlap_tol, use_ea = use_ea, use_lb = use_lb, debug = debug) 
            
            if(k > 0){
               kup <- kNN_update(kup = kup, bix = ret[[i]]$all_best_indices, 
                                 bdi = ret[[i]]$dist[ret[[i]]$all_best_indices + 1], 
                                 ki = ret[[i]]$kNN_inf, i = i)
            }
         }
         
      }else{
         
         for(i in seq_along(C)){
            ret[[i]] <- cpp_rundtw_mv_lot(h = Q, x = C[[i]], 
                                          kNN_val_in = kup$kNN_val, kNN_ix_in = kup$kNN_ix, kNN_lot_ix_in = kup$kNN_lot_ix, 
                                          kNN_inf_list_in = kup$kNN_inf_list, lot_ix = i,
                                          step_pattern = step_pattern, dist_method = dist_method, ws = ws, threshold = threshold, 
                                          overlap_tol = overlap_tol, do_norm = do_scale, use_ea = use_ea, use_lb = use_lb, debug = debug) 
            
            if(k > 0){
               kup <- kNN_update(kup = kup, bix = ret[[i]]$all_best_indices, 
                                 bdi = ret[[i]]$dist[ret[[i]]$all_best_indices + 1], 
                                 ki = ret[[i]]$kNN_inf, i = i)
            }
         }
      }
   }
   
   counter <- Reduce("+", lapply(ret, "[[", "counter"))
   names(counter) <- c("scale_reset","scale_new_extreme", "scale_1step",
                       "cm_reset", "cm_1step",
                       "early_abandon", "lower_bound")
   counter["completed"] <- sum(sapply(C, NROW) - NROW(Q) + 1) - (counter["early_abandon"] + counter["lower_bound"])
   
   if(scale == "z"){
      counter <- counter[c("early_abandon", "lower_bound", "completed")]
   }else if(scale == "none"){
      counter <- counter[c("cm_reset", "cm_1step","early_abandon", "lower_bound", "completed")]
   }
   
   ret <- list(dist = lapply(ret, "[[", "dist"),
               counter = counter)
   if(k > 0){
      ix99 <- which(kup$kNN_ix != -99)
      ret <- c(ret, list(knn_values = kup$kNN_val[ix99],
                         knn_indices = kup$kNN_ix[ix99] + 1,
                         knn_list_indices = kup$kNN_lot_ix[ix99]  ))
   }
   
   if(return_QC) ret <- c(ret, list(Q = Q, C = C))
   class(ret) <- append(class(ret), "rundtw", 0)
   return(ret)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


kNN_update <- function(kup, bix, bdi, ki, i){
   # kup...kNN update
   # bix ... best indices
   # bdi ... best distances
   # ki ... kNN_info
   # i ... runnind index
   
   k <- kup$kNN_inf_list$nr_looking4
   all_best_dist <- c(kup$kNN_val, bdi)
   new_order <- order(all_best_dist, decreasing = FALSE)
   
   kup$kNN_val    <- all_best_dist[new_order[1:k]]
   kup$kNN_ix     <- c(kup$kNN_ix , bix)[new_order[1:k]]
   kup$kNN_lot_ix <- c(kup$kNN_lot_ix , rep(i, length(bix) ))[new_order[1:k]]
   
   kup$kNN_inf_list$imax        <- 0
   #kup$kNN_inf_list$which_imax  <- 0# is set in cpp function
   kup$kNN_inf_list$which_vmax  <- which.max(kup$kNN_val) - 1
   kup$kNN_inf_list$nr_detected <- sum(kup$kNN_val < Inf)
   return(kup)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


find_peaks <- function (x, w, get_min = TRUE, strict = TRUE){
   
   cpp_strict <- strict * 1
   if(get_min){
      return( cpp_local_min(x, w = w, strict = cpp_strict) )
   }else{
      return( cpp_local_min(-x, w = w, strict = cpp_strict) )
   }
}



#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


lowerbound_tube <- function(Q, ws, scale = c("z", "01", "none")){
   
   scale <- match.arg(scale)
   if(scale != "none"){
      Q <- IncDTW::scale(Q, type = scale)
   }
   
   if(is.matrix(Q)){
      # stop("not implemented yet")
      return(cpp_get_tube_mv(Q, ws))
   }else{
      return(cpp_get_tube(Q, ws))   
   }
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


lowerbound <- function(C, ws, scale = c("z", "01", "none"),
                     dist_method = c("norm1", "norm2", "norm2_square"), 
                     Q = NULL, tube = NULL){
   
   if(is.null(Q) & is.null(tube)){
      stop("Either Q or tube must not be NULL")
   }
   
   dist_method <- match.arg(dist_method)
   scale <- match.arg(scale)
   
   # if(scale == "01"){
   #    warning("currently only for scale == 'z' implemented")
   # }
   
   if(is.null(tube)){
      if(scale != "none") Q <- IncDTW::scale(Q, type = scale)
      tube <- lowerbound_tube(Q, ws)
   }
   
   if(scale != "none"){
      C <- IncDTW::scale(C, type = scale)
   }
   
   # adjust the lb functions for other scaling method '01'
   if(is.matrix(C)){
      if(dist_method == "norm1"){
         lb <- get_lb_mv1(tube, C, 0, nrow(C), ncol(C))   
      
      }else if(dist_method == "norm2"){
         lb <- get_lb_mv2(tube, C, 0, nrow(C), ncol(C))
         
      }else if(dist_method == "norm2_square"){
         lb <- get_lb_mv22(tube, C, 0, nrow(C), ncol(C))
      
      }
   }else{
      lb <- get_lb(tube, C, 0, length(C))
   }
   
   return(lb)
}