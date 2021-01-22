.onLoad <- function(libname, pkgname)
{
  library.dynam("fclust", pkgname, libname)
}

.onUnload <- function (lib)
{
  library.dynam.unload("fclust", lib)
}

# fclustStartupMessage <- function()
# {
#   msg <- c(paste0(
# "fclust version"," ",(packageVersion("fclust")),"\n", "A toolbox for fuzzy clustering"), "\n", "Type 'citation(\"fclust\")' for citing this R package in publications.")
#   return(msg)
# }
#
# .onAttach <- function(lib, pkg)
# {
#   # startup message
#   msg <- fclustStartupMessage()
#   if(!interactive())
#     msg[1] <- paste("Package 'fclust' version", packageVersion("fclust"))
#   packageStartupMessage(msg)
#   invisible()
# }


