.onLoad <- function(libname, pkgname) {
	if ("package:parallel" %in% search()) {
		cores <- as.integer(Sys.getenv("OMP_NUM_THREADS"))
		if (is.na(cores)) {
      cores <- detectCores(logical=TRUE)
      if(is.na(cores)) cores <- 1L
			else cores <- 2L
    }
		.Call(`_rpf_setNumberOfCores`, cores)
		#packageStartupMessage(paste("OpenMP will use", cores, "cores"))
	}
}

.onAttach <- function(libname, pkgname) {
	if (! .Call(`_rpf_has_openmp`)) {
		packageStartupMessage("RPF is not compiled to take advantage of computers with multiple cores.")
	}
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call('_rpf_registerCCallable', PACKAGE = 'rpf')
})
