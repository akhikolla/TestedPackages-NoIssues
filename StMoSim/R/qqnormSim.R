myQQNorm <- function(x, nSim, mOfVar, main, xlab, ylab, qqnormCol, qqnormPch, qqlineCol, qqlineLwd){
  n_size <- length(x)
  
  if(n_size < 5){
    stop("qqnormSim - too small sample")
  }
  if(qqlineLwd > 3){
    stop("qqnormSim - too high qqlineLwd value")
  }
  
  n_mean <- mean(x)
  
  if(mOfVar == "mad"){
    n_sd <- mad(x)    
  }else if(mOfVar == "sd"){
    n_sd <- sd(x)      
  }else{
    stop('qqnormSim - no valid option for argument mOfVar - valid are "mad" or "sd".')
  }

  
  a <- ifelse(n_size <= 10, 3/8, 1/2)
  myseq <- seq(from = (1 - a)/(n_size + (1-a)-a), to = (n_size - a)/(n_size + (1-a)-a), length.out = min(n_size,50))
  
  n_quant <- myQQNormIntern(myseq,n_mean,n_sd,length(x),nSim)
  
  myX <- qnorm((1:n_size - a)/(n_size + (1-a)-a))
  myY <- sort(x)
  
  plot(myX, myY, main = main, xlab = xlab, ylab = ylab, col = qqnormCol, pch = qqnormPch) 
  matlines(qnorm(myseq), n_quant, lty = 1,pch = 1, lwd = qqlineLwd, col = qqlineCol)
  
  points(myX, myY, col = qqnormCol, pch = qqnormPch)
  qqline(x)
  box()
}


#' Quantile-Quantile plot with several Gaussian simulations.
#' 
#' Plots a QQ plot of the variable x with nSim Gaussian simulations.
#' 
#' Two estimators are required for the simulation of the normal distribution. Since the normal distribution is a two-parameter family distribution.
#' Default measure of location is the mean. Default measure of variation is the mad. This gives a robust estimation of the standard deviation even if there are outliers in the sample. 
#' Likewise this can be changed with the parameter \code{mOfVar}.
#' 
#' 
#' @param x a lm-object or a numeric vector. If it's a lm-object its residuals are plotted.
#' @param nSim \emph{[optional]} the number of simulations you like to add to the plot.
#' @param mOfVar \emph{[optinal]} a measure of variation. ("mad" or "sd")
#' @param main \emph{[optional]} an overall title for the plot.
#' @param xlab \emph{[optional]} a title for the x axis.
#' @param ylab \emph{[optional]} a title for the y axis.
#' @param qqnormCol \emph{[optional]} color of the obervations in the plot.
#' @param qqnormPch \emph{[optional]} point character of the observations in the plot.
#' @param qqlineCol \emph{[optional]} color of the simulations in the plot.
#' @param qqlineLwd \emph{[optional]} line width of the simulations. should not be higher than 3.
#'
#'
#' @return invisible(NULL)
#' 
#' @export
#' @rdname qqnormSim
#'
#' @useDynLib StMoSim
#'
#' @importFrom graphics box matlines plot points
#' @importFrom stats sd mad qnorm qqline resid
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @import methods 
#'
#' @examples
#' \dontrun{
#' 
#' ######## qqnorm vs. qqnormSim ########
#' 
#' par(mfrow = c(1,2))
#' x<- rnorm(100)
#' qqnorm(x)
#' qqline(x)
#' qqnormSim(x)
#' par(mfrow = c(1,1))
#' 
#' ######## basic functionality/arguments ########
#' 
#' # The observations should behave like a simulation, 
#' # because the observations are sampled from a Gaussian distribution.
#' qqnormSim(x = rnorm(100))
#' 
#' # If you don't feel comfortable with the mad as 
#' # measure of variation you can change it to the standard deviation.
#' qqnormSim(x = rnorm(100),
#'           mOfVar = "sd")
#' 
#' # On the first glance its obvious that this sample 
#' # doesn't originate from a Gaussian distribution due to the heavy tails.
#' qqnormSim(x = rt(100,df = 4))
#'
#' Reduce the simulation tracks from 500 to 50. (500 is default).
#' Not recommended unless you have not enough computation power.
#' qqnormSim(x = rnorm(100), 
#'           nSim = 50)
#' 
#' ######## graphical arguments ########
#' 
#' # set title and axes labels.
#' qqnormSim(x = rnorm(100), 
#'           main = "main title",
#'           xlab = "x-axis label",
#'           ylab = "y-axis label")
#'           
#' # I don't recommend fancy colors, unless you need it for your corporate identity.
#' qqnormSim(x = rnorm(100), 
#'           qqnormCol = "#ff0000",
#'           qqnormPch = 16,
#'           qqlineCol = "greenyellow",
#'           qqlineLwd = 1)
#' 
#' }
#' @keywords qqnorm
#' 
#' @seealso the basic graph corresponds to \link[stats]{qqnorm}
#' 
#' @author Matthias Salvisberg <matthias.salvisberg@@gmail.com>
#' 
setGeneric("qqnormSim", function(x, 
                                 nSim = 500,
                                 mOfVar = "mad",
                                 main = "Normal Q-Q Plot - SIM",
                                 xlab = "Theoretical Quantiles",
                                 ylab = "Sample Quantiles",
                                 qqnormCol = "black",
                                 qqnormPch = 1,
                                 qqlineCol = "#cdd2d015",
                                 qqlineLwd = 3) 
  standardGeneric("qqnormSim"))


#' @rdname qqnormSim
#' @aliases qqnormSim,lm-method
setMethod("qqnormSim","lm",
          function(x, 
                   nSim = 500,
                   mOfVar = "mad",
                   main = "Normal Q-Q Plot - SIM",
                   xlab = "Theoretical Quantiles",
                   ylab = "Sample Quantiles",
                   qqnormCol = "black",
                   qqnormPch = 1,
                   qqlineCol = "#cdd2d015",
                   qqlineLwd = 3){
            myQQNorm(resid(x), 
                     nSim,
                     mOfVar,
                     main,
                     xlab,
                     ylab,
                     qqnormCol,
                     qqnormPch,
                     qqlineCol,
                     qqlineLwd)
            return(invisible(NULL))
          })

#' @rdname qqnormSim
#' @aliases qqnormSim,numeric-method
setMethod("qqnormSim","numeric",
          function(x, 
                   nSim = 500,
                   mOfVar = "mad",
                   main = "Normal Q-Q Plot - SIM",
                   xlab = "Theoretical Quantiles",
                   ylab = "Sample Quantiles",
                   qqnormCol = "black",
                   qqnormPch = 1,
                   qqlineCol = "#cdd2d015",
                   qqlineLwd = 3){
            myQQNorm(x, 
                     nSim,
                     mOfVar,
                     main,
                     xlab,
                     ylab,
                     qqnormCol,
                     qqnormPch,
                     qqlineCol,
                     qqlineLwd)
            return(invisible(NULL))
          })
