.onAttach <- function(libname, pkgname)
   packageStartupMessage("densEstBayes 1.0 loaded.\nCopyright M.P. Wand 2020.\nFor details on the use of densEstBayes, issue the command:\ndensEstBayesVignette()")

.onUnload <- function(libpath)
    library.dynam.unload("densEstBayes",libpath)
