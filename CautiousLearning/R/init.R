.onLoad <- function(libname, pkgname) {
    if (hasOMP()) setOMPThreads(parallel::detectCores())
    setSITMOSeeds(runif(1))
}

