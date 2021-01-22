#' @title plot_matrix
#' @description Function to build a heat map displaying the values of given data matrix. This graph is useful for a visual inspection of the 
#' spatial distribution of the observations within the target matrix.
#' @param y data matrix to be evaluated.
#' @param standardize.rows logical argument (default = TRUE) indicating whether to standardize 
#' the rows of y to build the image.
#' @param reorder.rows logical argument (default = TRUE) indicating whether to reorder the rows 
#' of y to highlight a pattern.
#' @param reorder.cols logical argument (default = TRUE) indicating whether to reorder the columns 
#' of y to highlight a pattern.
#' @param high.contrast logical argument (default = TRUE) indicating whether to apply a transformation 
#' to increase contrast in the image of y.
#' @importFrom lattice levelplot
#' @importFrom grDevices colorRampPalette
#' @importFrom stats median
#' @importFrom stats quantile
#' @seealso
#' \code{\link{slfm}}, \code{\link{process_matrix}}, \code{\link{slfm_list}}
#' @export

plot_matrix <- function(y, standardize.rows = TRUE, reorder.rows = TRUE, reorder.cols = TRUE, high.contrast = TRUE) {

  y <- as.matrix(y)
  xl <- "Microarrays" 
  yl <- "Probes"
  
  if(standardize.rows) {
    y <- t(apply(y,1,scale)) # standardize the rows of y
  } 
  
  if(reorder.rows) {
    medrow <- apply(y,1,median) # median of the rows of y
    y <- y[order(medrow),] # order the rows of y 
  } 
  
  if(reorder.cols) {
    medcol <- apply(y,2,median) # median of the columns of y
    y <- y[,order(medcol)]  # order the columns of y
  }
  
  if(high.contrast) {
    y <- (abs(y)^(1/3))*sign(y)
  } # higher contrast
  
  mi <- as.numeric(quantile(y,0.001))
  ma <- as.numeric(quantile(y,0.999))
  nr <- nrow(y)
  nc <- ncol(y)

  if(nr==nc) { # Image of Correlation and Covariance matrices
   xl <- ""
   yl <- ""
   mi <- ifelse(ma==1, -1, ma)
  }

  spr <- ifelse(nr<=5, 1, round(nr/10))
  spc <- ifelse(nc<=5, 1, round(nc/10))
  sc <- list(
    x = list(at=c(seq(1,nc, spc),nc), labels=c(seq(1,nc,spc),nc)), 
    y = list(draw=FALSE)
  )
  col.l <- colorRampPalette(c('blue','white','red'))
  cbar <- seq(mi,ma,length.out=100)
  
  levelplot(t(y[nr:1,]),col.regions=col.l,xlab=xl,ylab=yl,scales=sc,at=cbar,aspect="fill")

}