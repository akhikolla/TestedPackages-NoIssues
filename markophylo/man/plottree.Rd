\name{plottree}
\alias{plottree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the tree used with different hypothesized branch groupings (or clades) following unique rates coloured differently.
}
\description{
Plotting command for use on an object of class "markophylo".
}
\usage{
plottree(x, colors = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class "markophylo".}
  \item{colors}{Default works well. However, a custom vector of colours can be specified here---should be the same length as length(x$bg). Note that these colours are used to colour the different branch groupings following their own rates.}
  \item{\dots}{Any further commands to ape::plot.phylo.}
}

\author{
Utkarsh J. Dang and G. Brian Golding

\email{udang@mcmaster.ca}}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{estimaterates}}.
}
\examples{
data(simdata2)
model2 <- estimaterates(usertree = simdata2$tree, userphyl = simdata2$data, 
                        alphabet = c(1, 2), bgtype = "ancestornodes", bg = c(7),
                        rootprob = "equal", modelmat = matrix(c(NA, 1, 2, NA), 2, 2))
plottree(model2, colors=c("blue", "darkgreen"), edge.width = 2, show.tip.label = FALSE, 
         no.margin = TRUE)
ape::nodelabels(frame = "circle", cex = 0.7)
ape::tiplabels(frame = "circle", cex = 0.7)
}