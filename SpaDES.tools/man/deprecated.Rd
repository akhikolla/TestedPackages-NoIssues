% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spades-tools-deprecated.R
\name{fastCrop}
\alias{fastCrop}
\title{\code{fastCrop} is a wrapper around \code{velox::VeloxRaster_crop}, though
\code{raster::crop} is faster under many tests.}
\usage{
fastCrop(x, y, ...)
}
\arguments{
\item{x}{Raster to crop}

\item{y}{Extent object, or any object from which an Extent object can be extracted (see Details)}

\item{...}{Additional arguments as for \code{\link[raster]{writeRaster}}}
}
\description{
\code{fastCrop} is a wrapper around \code{velox::VeloxRaster_crop}, though
\code{raster::crop} is faster under many tests.
}
\seealso{
\code{velox::VeloxRaster_crop}
}
