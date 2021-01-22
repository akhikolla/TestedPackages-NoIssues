\name{plot.densEstBayes}
\alias{plot.densEstBayes}
\title{Plot the Bayesian density estimate from a \code{densEstBayes()} fit}
\description{The estimated density function obtained via \code{densEstBayes} is plotted.}
\usage{
\method{plot}{densEstBayes}(x,plotIt=TRUE,credLev=0.95,gridSize=1001,varBand=TRUE,
                         shade=TRUE,estCol="darkgreen",varBandCol=NULL,
                         axisCol="slateblue",add=FALSE,lwd=2,xlab=NULL,ylab=NULL,...)
}
\arguments{
\item{x}{A \code{densEstBayes()} fit object.}
\item{plotIt}{Boolean flag:\cr
                TRUE = plot the density estimate and return nothing (the default)\cr
                FALSE = do not plot the density estimate and, instead, return a list containing an equally-spaced grid of abscissae values and three grids of ordinate values corresponding to the estimate and lower and uppwer limits of pointwise (100*credLev)\% sets.}
\item{credLev}{number between 0 and 1 such that the credible interval band has (100*credLev)\% approximate pointwise coverage. The default value is 0.95.}
\item{gridSize}{number of grid points used to display the density estimate curve and the pointwise credible interval band. The default value is 1001.}
\item{varBand}{Boolean flag:\cr
                TRUE = add a pointwise approximate (100*credLev)\% credible set variability band (the default)\cr
                FALSE = only plot the density estimate with a variability band.}
\item{shade}{Boolean flag:\cr
                TRUE = display the variability band using shading (the default)\cr
                FALSE = display the variability band using dashed curves.}
\item{estCol}{colour of the density estimate curve. The default value is "darkgreen".}
\item{varBandCol}{colour of the pointwise credible interval variability band. If \code{shade}=TRUE then the default value is "palegreen". If \code{shade}=FALSE then the default value is \code{estCol}.}
\item{axisCol}{colour of the horizontal axis. The default value is "slateblue".}
\item{add}{Boolean flag:\cr
                TRUE = add the density estimate curve(s) to an existing plot\cr
                FALSE = create a new plot for display of the density estimate (the default).}
\item{lwd}{positive integer indicating the line width of plotted curves. The default value is 2.}
\item{xlab}{character string spedifying the horizontal axis label.}
\item{ylab}{character string specifying the vertical axis label.}
\item{...}{place-holder for other graphic parameters.}
}
\value{If \code{plotIt}=TRUE then a plot is produced on the current device and no numerical values are returned. If \code{plotIt}=FALSE then a list is returned with the following components:
\item{xg}{numerical vector of abscissae values}
\item{densEstg}{numerical vector of ordinate values corresponding to the density estimate}
\item{densLowg}{numerical vector of ordinate values corresponding to the lower limits of the pointwise approximate (100*credLev)\% credible set variability band}
\item{densUppg}{numerical vector of ordinate values corresponding to the upper limits of the pointwise approximate (100*credLev)\% credible set variability band}
}
\author{Matt P. Wand \email{matt.wand@uts.edu.au}}
\examples{
library(densEstBayes) ; data(OldFaithful2011)

# Obtain a density estimate for the `OldFaithful2011' data:

dest <- densEstBayes(OldFaithful2011,method = "SMFVB")

# Plot the density estimate using default settings:

plot(dest,xlab = "time interval between geyser eruptions (minutes)")
rug(jitter(OldFaithful2011,amount=0.2),col = "dodgerblue")

# Plot the density estimate with some user-specified settings:

plot(dest,credLev = 0.999,estCol = "purple",
     varBandCol = "pink",axisCol = "navy",
     xlab = "time interval between geyser eruptions (minutes)")
rug(jitter(OldFaithful2011,amount= 0.2),col = "darkkhaki")

# Plot the density estimate as a black and white line plot:

plot(dest,estCol = "black",shade = FALSE,axisCol = "black",
     xlab = "time interval between geyser eruptions (minutes)")
rug(jitter(OldFaithful2011,amount = 0.2))
}




