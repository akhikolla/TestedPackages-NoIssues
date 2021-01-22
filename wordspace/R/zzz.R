## define .onLoad and .onAttach here if package initialisation functions are needed;
## .Last.lib (has to be exported) or .onUnload for package finalisation


##
## Internal helper functions
##

## extract (score or frequency) matrix from DSM object or pass through matrix; ensures canonical format
find.canonical.matrix <- function (x, triplet=FALSE) {
  if (inherits(x, "dsm")) {
    info <- check.dsm(x)
    M <- if (info$S$ok) x$S else x$M
  } else if (is.matrix(x) || is(x, "dMatrix")) {
    M <- x
  } else {
    stop("first argument must be an object of class 'dsm' or a co-occurrence/score matrix")
  }
  info <- dsm.is.canonical(M)
  if (info$canonical && !triplet) M else dsm.canonical.matrix(M, triplet=triplet)
}
