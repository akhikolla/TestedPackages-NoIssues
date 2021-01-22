#' Get matrix of observed transition frequencies
#'
#' Summarizes a sequence of discrete values by the observed transition frequencies.
#'
#' @param z vector of model indices (numerical or character)
#' @param labels fixed labels for models that should be included in transition matrix, e.g., \code{labels=1:20} or \code{c("m1","m2",...)}
#' @param order order of the transition table. If \code{order=1}, a matrix with transition frequencies from \code{z[t+1]} is returned. If \code{order=2}, a 3-dimensional array is returned with transition frequencies for \code{z[t]}, \code{z[t+1]}, and \code{z[t+2]}.
#' @return a square matrix with transition frequencies
#'
#' @examples
#' P <- matrix(c(.9,.1,0,
#'               .1,.6,.3,
#'               .2,.3,.5), 3, byrow=TRUE)
#' z <- rmarkov(2000, P)
#' transitions(z)
#' transitions(z, order = 2)
#'
#' @export
transitions <- function (z, labels, order = 1){
  if (inherits(z, "list") || inherits(z, "mcmc.list")){
    z <- do.call("cbind",z)
  }
  if (missing(labels) || is.null(labels)){
    labels <- sort(unique(as.vector(z)))
  } else if (!all(z %in% labels)){
    stop("'labels' does not include all elements in the sequence 'z'.")
  }

  if (is.matrix(z)){
    ### multiple chains across columns
    tabs <- lapply(split(z, rep(1:ncol(z), each = nrow(z))),
                   transitions, labels=labels, order = order)
    tab <- Reduce("+", tabs)

  } else {
    ### single chain
    n <- length(z)
    if (order == 1){
      tab <- table(factor(z[1:(n-1)], levels=labels),
                   factor(z[2:n], levels=labels))
      class(tab) <- "matrix"
      dimnames(tab) <- list("t" = labels,"t+1" = labels)
    } else {
      tab <- table(factor(z[1:(n-2)], levels=labels),
                   factor(z[2:(n-1)], levels=labels),
                   factor(z[3:n], levels=labels))
      class(tab) <- "array"
      dimnames(tab) <- list("t" = labels,"t+1" = labels,"t+2" = labels)
    }
  }
  tab
}

