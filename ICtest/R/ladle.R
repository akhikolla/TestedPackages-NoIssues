# Estimates the number of components after which the spectrum of S stays constant
# x = the data matrix, n x p
# S = a function for the matrix we're interested in
ladle <- function(x, S, n.boots = 200, ...){
  data.name <-  deparse(substitute(x))
  method <- "General"
  
  x <- as.matrix(x)
  
  n <- nrow(x)
  p <- ncol(x)
  eig_x <- eigen(S(x, ...))
  vec_x <- eig_x$vectors
  val_x <- eig_x$values
  
  # The size of the matrix
  q <- length(val_x)
  
  det_matrix <- matrix(0, n.boots, q - 1)
  
  for(b in 1:n.boots){
    x_b <- x[sample(1:n, n, replace = TRUE), ]
    vec_b <- eigen(S(x_b, ...))$vectors 
    
    vec_b_list <- lapply(1:(q - 1), function(i) vec_b[, 1:i, drop = FALSE])
    det_b <- unlist(lapply(vec_b_list, function(x) det(t(vec_x[, 1:(ncol(x))])%*%x)))
    
    det_matrix[b, ] <- det_b
  }
  
  fn0 <- c(0, colMeans(1 - abs(det_matrix)))
  fn <- fn0/(1 + sum(fn0))
  
  phin <- val_x/(1 + sum(val_x))
  gn <- fn + phin
  est.k <- which.min(gn) - 1
  
  # if(est.k != 0){
  #   W <- vec_x[, 1:est.k]
  #   S <- x%*%W
  # }
  # else{
  #   S <- NULL
  # }
  
  res <- list(method = method,
              k = est.k,
              fn = fn,
              phin = phin,
              gn = gn,
              # S = S,
              lambda = val_x,
              data.name = data.name)
  class(res) <- "ladle"
  return(res)
}
