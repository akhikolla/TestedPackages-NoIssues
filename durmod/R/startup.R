.onLoad <- function(libname, pkgname) {

  if(is.null(cr <- getOption('durmod.threads'))) {
    cr <- as.integer(Sys.getenv('DURMOD_THREADS'))
    if(is.na(cr)) cr <- as.integer(Sys.getenv('OMP_NUM_THREADS'))
    if(is.na(cr)) cr <- as.integer(Sys.getenv('OMP_THREAD_LIMIT'))
    if(is.na(cr)) cr <- as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))
    if(is.na(cr)) cr <- parallel::detectCores(all.tests=TRUE)
    if(is.na(cr) || cr < 1) {
      cr <- 1
    }
    options(durmod.threads=cr)
  }

}

# set initial options, with defaults. Don't set if already set.
setoption <- function(...) {
  arg <- list(...)
  # don't change an existing option
  arg <- arg[sapply(paste('lfe',names(arg),sep='.'), function(n) is.null(getOption(n)))]
  # fill in from environment variable
  if(length(arg) == 0) return()
  nm <- names(arg)
  arg <- lapply(nm, function(n) {
    e <- Sys.getenv(paste('DURMOD',toupper(n),sep='_'))
    if(e != '') {
      val <- try(eval(parse(text=e)))
      if(inherits(val, 'try-error')) val <- arg[[n]]
      val
    } else arg[[n]]
  })
  names(arg) <- paste('lfe',nm,sep='.')
  do.call(options, arg)
}
