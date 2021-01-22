## Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
merge_list = function(x, y) {
  x[names(y)] = y
  x
}

new_defaults = function(value = list()) {
  defaults = value
  
  get = function(name, default = FALSE, drop = TRUE) {
    if (default) defaults = value  # this is only a local version
    if (missing(name)) defaults else {
      if (drop && length(name) == 1) defaults[[name]] else {
        setNames(defaults[name], name)
      }
    }
  }
  resolve = function(...) {
    dots = list(...)
    if (length(dots) == 0) return()
    if (is.null(names(dots)) && length(dots) == 1 && is.list(dots[[1]]))
      if (length(dots <- dots[[1]]) == 0) return()
    dots
  }
  set = function(...) {
    dots = resolve(...)
    if (length(dots)) defaults <<- merge(dots)
    invisible(NULL)
  }
  merge = function(values) merge_list(defaults, values)
  restore = function(target = value) defaults <<- target
  list(get = get, set = set, restore = restore)
}

#' Options for stfit
#'
#' @export
opts_stfit = new_defaults(list(
  temporal_mean_est = smooth_spline
))

#' Remove outlier
#'
#' An outlier is defined as points outside the whiskers of the boxplot
#' over the time domain (DOY).
#'
#' @param rst a *Raster object
#'
#' @return a *Raster object
#' @export
#'
rmOutlier <- function(rst){
  .rmOutlier <- function(x){
    na.idx = is.na(x)
    boxstat = boxplot(x[!na.idx], plot=FALSE)
    x[!na.idx][(x[!na.idx] > boxstat$stats[5,]) | (x[!na.idx] < boxstat$stats[1,])] = NA
    x
  }
  calc(rst, .rmOutlier)
}

#' Missing value percentages
#'
#' @param x A \code{RasterStack} object
#' @param mc.cores Numer of cores to use
#'
#' @return A vector of percent of missing values for each layer
#' @export
#'
pctMissing <- function(x, mc.cores){
  if(missing(mc.cores)) mc.cores = parallel::detectCores()
  doParallel::registerDoParallel(cores=mc.cores)
  nl = nlayers(x)
  i = NULL
  foreach::foreach(i=1:nl, .combine = c) %dopar%{
    val = values(x[[i]])
    sum(is.na(val))/length(val)
  }
}

## test whether all input elemetns are equal
ident <- function(...) {
  args <- c(...)
  if (length(args) > 2L) {
    #  recursively call ident()
    out <- c(identical(args[1] , args[2]) , ident(args[-1]))
  } else{
    out <- identical(args[1] , args[2])
  }
  return(all(out))
}

#' Get missing layer index
#'
#' @param rst.list a RasterStack or RasterBrick object or a list of them
#'
#' @return index of the missing layers
#' @export
#'
getMissingLayers <- function(rst.list){
  if(inherits(rst.list, "RasterStackBrick"))
    return(which(is.infinite(rst.list@data@min) | is.na(rst.list@data@min)))
  if(is.list(rst.list))
    return(lapply(1:length(rst.list), function(i) 
      which(is.infinite(rst.list[[i]]@data@min) | is.na(rst.list[[i]]@data@min))))
}



#' Get image mask
#' 
#' @param object A numeric matrix. Each row is an row stacked image.
#'
#' @param tol If the percentage of missing values for a pixel over time is greater than this
#' value, this pixel is treated as a mask value.
getMask <- function(object, tol = 0.95) {
  return(
    apply(object, 2, function(x) {
      sum(is.na(x))/length(x) >= tol
    }))
}

#' Epanicnicov kernel function
#'
#' @param x numeric vector
#'
#' @return vector
#' @export
#'

epan <- function(x) {
  3 / 4 * pmax(1 - x ^ 2, 0)
}

covInterp = function(tt, phi.fun, omega, nugg.fun, t.grid) {
  m = length(tt)
  J = length(t.grid)
  kpc = length(omega)
  phi = matrix(0, m, kpc)
  nugg = rep(0, m)
  
  for (j in 1:m) {
    s = tt[j]
    if (s == t.grid[1]) {
      phi[j, ] = phi.fun[1, ]
      nugg[j] = nugg.fun[1]
      
    } else if (s == t.grid[J]) {
      phi[j, ] = phi.fun[J, ]
      nugg[j] = nugg.fun[J]
    } else{
      ind = min(which(t.grid >= s))
      p = (s - t.grid[ind - 1]) / (t.grid[ind] - t.grid[ind - 1])
      phi[j, ] = (1 - p) * phi.fun[ind - 1, ] + p * phi.fun[ind, ]
      nugg[j] = (1 - p) * nugg.fun[ind - 1] + p * nugg.fun[ind]
    }
  }
  Vmat=phi%*%diag(omega)%*%t(phi)
  diag(Vmat)=diag(Vmat)+nugg
  return(Vmat)
}

phiInterp = function(tt, phi.fun, t.grid){
  m = length(tt)
  tt = tt[tt < max(t.grid)]
  tt.idx = findInterval(tt, t.grid)
  p = (tt - t.grid[tt.idx])/(t.grid[tt.idx+1] - t.grid[tt.idx])
  phi = (1 - p)*phi.fun[tt.idx,,drop=FALSE] + p*phi.fun[tt.idx+1,,drop=FALSE]
  if(length(tt) < m)
    phi = rbind(phi, phi.fun[rep(length(t.grid), m - length(tt)),])
  return(phi)
}

#' Weight matrix calculation
#' 
#' @param h 'bandwith'
#'
#' @return a weighting matrix
#' @export
#' @example 
#' weightMatrix(3)
weightMatrix <- function(h){
  if(h <=0)
    stop("bandwidth can not be non-positive.")
  if(h == 1)
    warning("No neighborhood information is used")
  aa = seq(-h, h, by = 1)/h
  out = outer(aa, aa, function(x, y) epan(sqrt(x^2+y^2)))
  return(out[2:(2*h), 2:(2*h)])
}

#' Weight vector calculation
#'
#' @param h bandwidth, should be positive numbers
#'
#' @return a vector
#' @export
#' @example 
#' weightVector(5)
weightVector <- function(h){
  if(h <=0)
    stop("bandwidth can not be non-positive.")
  if(h == 1)
    warning("No neighborhood information is used")
  out = epan(seq(-h, h, by = 1)/h)
  return(out[2:(2*h)])
}




