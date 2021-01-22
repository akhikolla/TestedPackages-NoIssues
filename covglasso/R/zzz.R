covglassoStartupMessage <- function()
{
  msg <- c(paste0(
"                     _
                    | |
  ___ _____   ____ _| | __ _ ___ ___  ___
 / __/ _ \\ \\ / / _` | |/ _` / __/ __|/ _ \\
| (_| (_) \\ V / (_| | | (_| \\__ \\__ \\ (_) |
 \\___\\___/ \\_/ \\__, |_|\\__,_|___/___/\\___/
                __/ |
               |___/
", "Version ", packageVersion("covglasso")),
  "\nType 'citation(\"covglasso\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- covglassoStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'covglasso' version", packageVersion("covglasso"))
  packageStartupMessage(msg)
  invisible()
}
