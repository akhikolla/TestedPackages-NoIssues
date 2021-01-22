# Wrapped function to estimate alphas
alphaEstFam <- function(gtime, delta) {
  gtimeFactor <- sort(unique(gtime))
  .Call("_groupedSurv_alphaEst1", PACKAGE = "groupedSurv", gtimeFactor, gtime, 
    delta)
}

