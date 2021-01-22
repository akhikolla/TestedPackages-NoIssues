#' @title Draw Kendall plot and compute AUK.
#' @description This function draws Kendall plot of 2 variables. Also provides an index AUK (area under Kendall plot).
#'
#' @param x a numeric vector stores first variable.
#' @param y a numeric vector stores second variable.
#' @param plot a TRUE/ FALSE flag for generating Kendall plot or not.
#' @param main a character indicating the title of the plot.
#' @param Auxiliary.line a TRUE/ FALSE flag for drawing auxiliary lines or not.
#' @param BS.CI a numeric specifying alpha for Bootstrap confidence interval. When euqal 0, confidence interval won't be computed.
#' @param set.seed a TRUE/ FALSE flag specifying setting seed or not.
#' @author Jeffrey C. Miecznikowski, En-shuo Hsu, Yanhua Chen, Albert Vexler
#' @details
#' AUK is bounded between 0 and 0.75. For positively correlated x and y's, say x = y, AUK = 0.75. And the plot follows the concave auxiliary line. While negatively correlated x and y's, AUK = 0. The plot is horizontal on y = 0. For independent x and y, AUK = 0.5. Kendall plot is on the diagonal.
#' Due to possible variable overflow, this function is only suitable for input size less than 1000. Input size greater than 1000 causes error.
#'
#' @return
#' a list containing a numeric AUK, a numeric vector W.in (x axis of plot), a numeric vector Hi.sort (y axis of plot), and three confidence intervals: normal CI, pivotal CI and percentage CI.
#'
#' @references
#' Vexler, Albert, Xiwei Chen, and Alan D. Hutson. "Dependence and independence: Structure and inference." Statistical methods in medical research (2015): 0962280215594198.
#'
#' R package "VineCopula":
#' Schepsmeier, Ulf, et al. "Package 'VineCopula'." (2015).
#'
#' @export
#' @examples
#' set.seed(123)
#' x = runif(100)
#' y = runif(100)
#'
#' result = AUK(x, y, plot = TRUE)
#' result$AUK
#'
#' #[1] 0.4987523

AUK = function(x, y, plot = F, main = "Kendall plot", Auxiliary.line = T, BS.CI = 0, set.seed = FALSE){
  #check input
  if(!is.numeric(x) || !is.numeric(y) || length(x)<2 || length(y)<2)
    stop("x and y must be numeric vectors.")
  if(length(x) > 1000)
    stop("Input size must be less than 1000.")
  if(length(x) != length(y))
    stop("x and y must be paired (with equal length).")
  if(!is.logical(plot))
    stop("plot must be of type logical")
  if(!is.character(main))
    stop("main must be of type character.")
  if(BS.CI<0 || BS.CI >1)
    stop("BS.CI must be between 0 and 1.")

  #compute x, y's copula

  n = length(x)
  cx = rank(x, ties.method = "min")/n
  cy = rank(y, ties.method = "min")/n

  #CI
  if(BS.CI != 0){
    times = 1000
    BS = function(i){
      if(set.seed){set.seed(i)}
      id = sample(1:n, n, replace = T)
      sx = x[id]
      sy = y[id]
      return (kPlot(sx, sy))
    }

    result = lapply(1:times, BS)
    Ts = c()
    for(i in 1:times){
      Ts[i] = result[[i]]$AUK
    }
    out = kPlot(cx, cy, PLOT = plot, main = main, Auxiliary.line = Auxiliary.line)
    result = c(out, getCI(out$AUK, Ts, BS.CI))
    return (result)
  }

  else{
    return(kPlot(cx, cy, PLOT = plot, main = main, Auxiliary.line = Auxiliary.line))
  }
}

my_rank = function(x){
  uni = unique(sort(x))
  out = c()
  for(i in 1:length(x)){
    out[i] = which(uni == x[i])
  }

  return (out)
}

Hn = function(x, y){
  n = length(x)
  table = emcdf_output(x, y, TRUE)
  rankX = my_rank(x)
  rankY = my_rank(y)
  out = c()
  for(i in 1:n){
    out[i] = table[rankY[i],rankX[i]]-1
  }
  return(out/(n-1))
}

kPlot = function(x, y, PLOT = F, main, Auxiliary.line){

    Hi <- Hn(x, y)
    Hi.sort <- sort(Hi)
    n <- length(x)
    W.in <- rep(NA, n)
    for (i in 1:n) {
        # function to be integrated
        f <- function(w) {
            w * (-log(w)) * (w - w * log(w))^(i - 1) * (1 - w + w * log(w))^(n - i)
        }
        W.in[i] <- n * choose(n - 1, i - 1) * (integrate(f, lower = 0, upper = 1)$value)
    }
    g <- function(w) {
        w - w * log(w)
    }  # K_{0}(w)=P(UV<=w)

    if (PLOT) {
      # should the results be plotted?
        if(Auxiliary.line){
          plot(g, main = main, xlim = c(0, 1), ylim = c(0, 1), pch = "x", xlab = expression(W[1:n]),
               ylab = "H")  #Kurve K_{0}(w)
          points(W.in, Hi.sort, pch = "x", cex = 0.4)
          abline(a = 0, b = 1)  # angle bisector
        }
        else{
          plot(W.in, Hi.sort, cex = 0.4, main = main, xlim = c(0, 1), ylim = c(0, 1), pch = "x", xlab = expression(W[1:n]),
             ylab = "H")  #Kurve K_{0}(w)
        }
    }

      #compute AUK
      W <- c(0,W.in,1)
      H <- c(0,Hi.sort,1)
      idx <- 2:length(W)
      area <- ((W[idx]-W[idx-1])%*%(H[idx]+H[idx-1]))/2
      AUK = min(area,0.75)

      # create output data
      kendall.plot.output <- list(AUK, W.in, Hi.sort)
      names(kendall.plot.output) <- c("AUK","W.in", "Hi.sort")
      return(kendall.plot.output)

}


