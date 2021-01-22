is.rundtw <- function(x) {
   inherits(x, "rundtw")
}

print.rundtw <- function(x, ...){

   if(!is.null(x$counter)){
      cat("counter:\n")
      print(x$counter, ...)
      cat("\n")
   }
   
   if(!is.null(x$knn_values)){
      cat("knn_values, distances of k nearest neighbors:\n")
      print(x$knn_values, ...)
      cat("\n")
   }
   
   if(!is.null(x$knn_indices)){
      cat("knn_indices, indices of k nearest neighbors:\n")
      print(x$knn_indices, ...)
      cat("\n")
   }
   
   if(!is.null(x$knn_list_indices)){
      cat("knn_list_indices, indices of list entries of C, where to find the kNN:\n")
      print(x$knn_list_indices, ...)
      cat("\n")
   }
   
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

is.dba <- function(x) {
   inherits(x, "dba")
}

print.dba <- function(x, digits = getOption("digits"), ...){
   
   if(!is.null(x$iterDist_m2m)){
      cat("iterDist_m2m, distances of barycenter_n to barycenter_n-1:\n")
      print(round(x$iterDist_m2m, digits = digits, ...), digits = digits, ...)
      # print(x$iterDist_m2m, ...)
      cat("\n")
   }
   
   if(!is.null(x$m1)){
      cat("m1, barycenter:\n")
      # print(trim_output_matrix(x$m1, ...), ...)
      print(trim_output_matrix(x$m1, digits = digits, ...), digits = digits, ...)
      cat("\n")
   }
   
   if(!is.null(x$input)){
      cat("Input:\n")
      lapply(names(x$input), function(nam){
         if(!is.null(x$input[[nam]])){
            cat(paste0(nam, ": "))
            cat(x$input[[nam]])
            cat("\n")
         }
      })
      cat("\n")
   }
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

is.idtw <- function(x) {
   inherits(x, "idtw")
}

print.idtw <- function(x, digits = getOption("digits"), ...) {
   
   if(!is.null(x$dm)){
      cat("dm, direction matrix:\n")
      print(trim_output_matrix(x$dm, digits = digits, ...), digits = digits, ...)
      cat("\n")
   }
   
   if(!is.null(x$diffM)){
      cat("diffM, matrix of differences:\n")
      print(trim_output_matrix(x$diffM, digits = digits, ...), digits = digits, ...)
      cat("\n")
   }
   
   if(!is.null(x$cm)){
      cat("cm, cost matrix:\n")
      print(trim_output_matrix(x$cm, digits = digits, ...), digits = digits, ...)
      cat("\n")
   }
   
   if(!is.null(x$gcm)){
      cat("gcm, global cost matrix:\n")
      print(trim_output_matrix(x$gcm, digits = digits, ...), digits = digits, ...)
      cat("\n")
   }
   
   
   cat("DTW distance: \n")
   cat(x$distance)
   cat("\n\nNormalized DTW distance: \n")
   cat(x$normalized_distance)
   cat("\n")
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

is.planedtw <- function(x) {
   inherits(x, "planedtw")
}

print.planedtw <- function(x, digits = getOption("digits"), ...) {
   
   cat("control: \n")
   y <- lapply(x$control, function(yy){ifelse(is.null(yy), "NULL", yy)})
   print(as.data.frame(y), row.names = FALSE, ...)
   cat("\n")
   cat("DTW distance: \n")
   cat(x$distance)
   cat("\n\nNormalized DTW distance: \n")
   cat(x$normalized_distance)
   cat("\n")
   
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

trim_output_matrix <- function(x, digits = getOption("digits"), is_dm = FALSE, ...){
   nc <- ncol(x)
   nr <- nrow(x)
   max_cols <- 7
   max_rows <- 7
   
   if(nc > max_cols){
      c1 <- 1
      c2 <- floor(max_cols/2)
      c3 <- nc - c2 + 1
      c4 <- nc
   }else{
      c1 <- 1
      c2 <- 1
      c3 <- 2
      c4 <- nc
   }
   
   if(nr > max_rows){
      r1 <- 1
      r2 <- floor(max_rows/2)
      r3 <- nr - r2 + 1
      r4 <- nr
   }else{
      r1 <- 1
      r2 <- 1
      r3 <- 2
      r4 <- nr
   }
   
   y <- x[c(r1:r2, r3:r4), c(c1:c2, c3:c4)]
   y <- as.data.frame(y)
   y <- round(y, digits = digits)
   y <- format(y, ...)
   
   ncy <- ncol(y)
   nry <- nrow(y)
   if(nc > max_cols){
      y <- cbind(y[, 1:(ncy/2)], "-" ="...", y[, (1 + ncy/2):ncy] )
      colnames(y) <- c(c1:c2, "-", c3:c4)
   }else{
      colnames(y) <- c(c1:c2, c3:c4)
   }
   if(nr > max_rows){
      y <- rbind(y[1:(nry/2), ], "...", y[(1 + nry/2):nry, ])
      rownames(y) <- c(r1:r2, "-" , r3:r4)
   }else{
      rownames(y) <- c(r1:r2, r3:r4)
   }
   return(y)
}
