\name{plot.dist.matrix}
\alias{plot.dist.matrix}
\title{
   Plotting Distance Matrices (wordspace)
}
\description{

  Visualization of a DSM distance matrix as a neighbourhood graph based on multidimensional scaling (MDS).

}
\usage{
  \method{plot}{dist.matrix}(x, y, labels=rownames(x), show.labels=TRUE, label.pos=3,
     selected=attr(x, "selected"), show.selected=TRUE,
     col="black", cex=1, pch=20, pt.cex=1.2, selected.cex=1.2, selected.col="red",
     show.edges=TRUE, edges.lwd=6, edges.col="#AABBFF", edges.threshold=quantile(x, 2/3),
     method=c("isomds", "sammon"), aspect=1, expand=.05, \dots)
}
\arguments{

  \item{x}{a symmetric distance matrix of class \code{dist.matrix}.  NB: similarity values and asymmetric distance measures are not supported.}
  
  \item{y}{unused, must not be specified}

  \item{labels}{a character vector of labels for the DSM vectors (defaults to rownames of \code{x})}

  \item{show.labels}{if \code{TRUE} (default), labels are displayed if available}

  \item{label.pos}{position of labels (default: above points).  Possible values are \code{1} (below), \code{2} (left), \code{3} (above) and \code{4} (right).}
  
  \item{selected}{logical vector of selected points that will be highlighted (defaults to optional \code{selected} attribute of distance matrix)}

  \item{show.selected}{if \code{TRUE} (default), points marked by \code{selected} are highlighted in the plot}

  \item{col}{colour of points and labels}

  \item{cex}{numeric character expansion factor for points and labels}

  \item{pch}{plot symbol for points}

  \item{pt.cex}{character expansion factor for points relative to labels}

  \item{selected.cex}{additional character expansion factor for selected points and labels}

  \item{selected.col}{colour of selected points and labels (if \code{show.selected=TRUE})}

  \item{show.edges}{if \code{TRUE} (default), edges are drawn between points.  The line width of each edge is proportional to the distance between the corresponding points.}

  \item{edges.lwd}{maximal line width of edges (for \eqn{d = 0})}

  \item{edges.col}{colour of edges, usually a light or translucent shade}

  \item{edges.threshold}{maximal distance up to which edges are drawn.  The default is to display two thirds of all edges.}

  \item{method}{whether to perform non-metric (\code{isomds}) or metric (\code{sammon}) multidimensional scaling}

  \item{aspect}{aspect ratio of plot window (e.g. \code{aspect=16/10} for a window that is 8 inches wide and 5 inches high). Setting a correct aspect ratio ensures that the distances between points in the MDS map are correctly represented in the plot.}

  \item{expand}{fraction by which plotting region is extended on each side. Adjust this parameter so that points and labels are completely visible.}

  \item{\dots}{all other arguments are passed to the initial plot function, which sets up the display but does not draw any graphical elements}

}
\value{

  Invisibly returns a two-column matrix with MDS coordinates of all displayed points and \code{labels} as rownames (if available).
  
}
\details{

  For multidimensional scaling (MDS), the functions \code{\link[MASS]{isoMDS}} and \code{\link[MASS]{sammon}} from the \bold{MASS} package are used.
  
}
\seealso{

  \code{\link{nearest.neighbours}}, which produces distance matrices suitable for plotting if the option \code{dist.matrix=TRUE} is specified

}
\examples{
\dontrun{
plot(nearest.neighbours(DSM_Vectors, "walk_V", n=20, dist.matrix=TRUE))
}
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
