#####
#
# general internal functions for ladle
#
########

# Function to compute for a specific iteration fn
# input
# EVboot - eigenvectors for a bootstrapping M matrix
# EVdata - eigenvectors for the M matrix for the observed data
# r - the value for the number of components to consider

fi <- function(EVboot, EVdata, r)
  {
  fni <- numeric(r)
  for (ii in 1:r) {
    fni[ii] <- det(crossprod(EVdata[,1:ii],EVboot[,1:ii]))
  }
  1-abs(fni)
  }

# print method for ladle object.
# prints everything but "S"
print.ladle <- function(x, ...)
  {
  print.listof(x[names(x) != "S"], ...)
  }

# plot method for ladle object.
# ploting function is the same as plot.ictest
plot.ladle <- function (x, which = "all", ...)
  {
    which <- match.arg(which, c("all", "k"))
    if (which == "all") {
      S <- x$S
    }
    else {
      S <- x$S[, 0:x$k, drop = FALSE]
    }
    if (ncol(S) <= 2) {
      plot(S, ...)
    }
    else {
      if (any(class(S) %in% c("mts", "xts", "zoo"))) {plot(S, ...)} else {pairs(S, ...)}
    }
  }

summary.ladle <- function(object, ...)
{
  cat("\n")
  cat(object$method, "Ladle estimator for", object$data.name,"\n")
  cat("k: ", object$k)
  cat("\n")
}


# extracting the components is the same as for ictest
components.ladle <- function (x, which = "all", ...)
{
    which <- match.arg(which, c("all", "k"))
    if (which == "all") {
        S <- x$S
    }
    else {
        S <- x$S[, 0:x$k, drop = FALSE]
    }
    S
}
