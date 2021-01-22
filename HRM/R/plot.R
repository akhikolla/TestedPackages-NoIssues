####################################################################################################################################
### Filename:    plot.R
### Description: Function for plotting the profiles only when one whole- and one subplot factor are used.
###
###
###
####################################################################################################################################

#' Plots profiles of the groups in case of one whole- and one subplot-factor.
#'
#' @param data A data.frame containing the data
#' @param group column name within the data frame data specifying the groups
#' @param factor1 column name within the data frame data specifying the first subplot-factor
#' @param subject column name within the data frame X identifying the subjects
#' @param response column name within the data frame X containing the response variable
#' @param xlab label of the x-axis of the plot
#' @param ylab label of the y-axis of the plot
#' @param legend logical indicating if a legend should be plotted
#' @param legend.title title of the legend
#' @return Plots profiles of the groups.
#' @example R/example_plot.txt
#' @keywords internal
hrm.plot <- function(data, group , factor1, subject, response, xlab="time", ylab="mean", legend = TRUE, legend.title = NULL ){
  X <- as.data.frame(data)
  data<-response
  stopifnot(is.data.frame(X),is.character(subject), is.character(group),is.character(factor1),is.character(data),is.character(xlab),is.character(ylab))
  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0

  group <- as.character(group)
  factor1 <- as.character(factor1)
  subject <- as.character(subject)
  xlab <- as.character(xlab)
  ylab <- as.character(ylab)
  X <- split(X, X[,group], drop=TRUE)
  a <- length(X)
  d <- nlevels(X[[1]][,factor1])
  n <- rep(0,a)

  means <- data.frame(dimension=as.numeric(levels(X[[1]][,factor1])))
  groupnames <- c()


  for(i in 1:a){
    groupnames[i] <- as.character(X[[i]][1,group])
    X[[i]] <- X[[i]][ order(X[[i]][,subject], X[[i]][,factor1]), ]
    X[[i]] <- X[[i]][,data]
    X[[i]] <- matrix(X[[i]], ncol=d,byrow=TRUE)
    n[i] <- dim(X[[i]])[1]
    means[,(i+1)] <- colMeans(X[[i]])
  }

  colnames(means) <- c("dimension",groupnames)

  means <- melt(means, id.vars="dimension")
  colnames(means) <- c("dimension", "group", "value")

  pl <- ggplot() +
        geom_line(data=means, aes(x=means$dimension, y=means$value,group=means$group,colour=means$group)) +
        geom_point(data=means, aes(x=means$dimension, y=means$value,group=means$group,colour=means$group),size=1.5) +
        xlab(xlab) +
        ylab(ylab)

  if(!legend){
    pl <- pl + theme(legend.position = "none")
  } else {
    if(!is.null(legend.title) & is.character(legend.title)){
      pl <- pl + scale_colour_hue(name=legend.title)
    } else {
      pl <- pl + theme(legend.title = element_blank())
    }
    pl <- pl + theme(legend.background = element_rect())
  }

  return(pl)

}

# hrm.plot end ------------------------------------------------------------
