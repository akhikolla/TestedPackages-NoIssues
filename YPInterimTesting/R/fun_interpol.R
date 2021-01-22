fun_interpol <- function(alphas, f) {
  
  y1 <- qnorm(alphas/2, lower.tail = FALSE)
  x1 <- log(2*pnorm(y1, lower.tail = FALSE))
  
  ty1 <- y1
  tx1 <- f(ty1)
  ty2 <- y1*(1.00 - 0.02)
  tx2 <- f(ty2)
  ty3 <- y1*(1.00 + 0.02)
  tx3 <- f(ty3)

  ind_y <- which.min((c(tx2, tx3) - log(alphas))^2)
  
  if (ind_y == 1) {
    y2 <- ty2
    x2 <- tx2
  } else {
    y2 <- ty3
    x2 <- tx3
  }
  
  ny1 <- (y2 - y1)/(x2 - x1)*(log(alphas) - x1) + y1
  nx1 <- f(ny1)
  ncr <- (y2 - ny1)/(x2 - nx1)*(log(alphas) - nx1) + ny1; 
  
  cr2 <- c(y1, ny1, y2)
  pt2 <- c(x1, nx1, x2)
  
  ind_y <- which.min((pt2 - log(alphas))^2)
  crf <- cr2[ind_y]

  return(crf)
}




