\name{predict.densEstBayes}
\alias{predict.densEstBayes}
\title{Obtain the Bayesian density estimate from a \code{densEstBayes()} fit at new abscissae}
\description{Values of the estimated density function, obtained via \code{densEstBayes}, are computed for new abscissae.}
\usage{
\method{predict}{densEstBayes}(object,newdata,cred.fit = FALSE,credLev = 0.95,...)
}
\arguments{
\item{object}{A \code{densEstBayes()} fit object.}
\item{newdata}{numerical vector of abscissae va}
\item{cred.fit}{Boolean flag:\cr
TRUE = compute the lower and upper limits of (100*credLev)\% pointwise credible intervals at 'newdata',\cr
FALSE = do not computer credible interval limits (the default).}
\item{credLev}{number between 0 and 1 such that the credible interval band has (100*credLev)\% approximate coverage. The default value is 0.95.}
\item{...}{place-holder for other graphic parameters.}
}
\value{A list is returned with the following components:
\item{fit}{numerical vector of ordinate values corresponding to the density estimate.}
\item{credLow.fit}{numerical vector of ordinate values corresponding to the lower limits of the pointwise approximate (100*credLev)\% credible set variability band.}
\item{credUpp.fit}{numerical vector of ordinate values corresponding to the upper limits of the pointwise approximate (100*credLev)\% credible set variability band.}
}
\author{Matt P. Wand \email{matt.wand@uts.edu.au}}
\examples{
library(densEstBayes) ; data(OldFaithful2011)

# Obtain a density estimate for the `OldFaithful2011' data:

dest <- densEstBayes(OldFaithful2011,method = "SMFVB")

# Plot the density estimate:

plot(dest,xlab = "time interval between geyser eruptions (minutes)")
rug(jitter(OldFaithful2011,amount = 0.2),col = "dodgerblue")

# Obtain predictions at 60,70,80,90,100 and 110 seconds and
# add to them plot:

newdataVec <- seq(60,110,by = 10)
predictObj <- predict(dest,newdata = newdataVec,cred.fit = TRUE)
print(predictObj$fit)
points(newdataVec,predictObj$fit,col = "blue")

# Print and add to the plot the lower and upper limits of 
# the pointwise 95% credible intervals:

print(predictObj$credLow.fit)
print(predictObj$credUpp.fit)
points(newdataVec,predictObj$credLow.fit,col = "red")
points(newdataVec,predictObj$credUpp.fit,col = "red")
}




