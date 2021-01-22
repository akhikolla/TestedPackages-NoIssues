\name{plot2Densities}
\alias{plot2Densities}
\title{Plotting the propability densitiy functions of two groups}
\description{
plots the pdf of the two groups as Pareto density estimation.
}
\usage{
plot2Densities(Data,Cls,col=c("red","blue"), pde=TRUE, 
	meanLines=FALSE,medianLines=FALSE,...)
}
\arguments{
 \item{Data}{The data of both groups as a vector.}
 \item{Cls}{The class information as a vector of similar length as the data vector.}
 \item{col}{Colors of the two lines to plot.}
 \item{pde}{Plots a standard pdf insted of a non-parametric Pareto density estimation.}
 \item{meanLines}{Adds perpendicular lines at the means of the two groups.}
 \item{medianLines}{Adds perpendicular lines at the medians of the two groups.}
 \item{\dots}{More graphical parameters can be entered as plot arguments if PlotIt = TRUE.}
}
\references{
 Lotsch, J., and Ultsch, A. (2019): 
 Impact -- an R Package for calculation and visualisation of the Impact distance and shape based effect size measure. 
}
\author{
  Jorn Lotsch and Alfred Ultsch
}
\examples{
## example 1
nrPerClass = 100
SampleImpact  = c()
SampleClasses = rep(c(1,2), each=nrPerClass)

SampleData = matrix(ncol=100, nrow=200)

for(i in 1:100){
  SampleVector = c(rnorm(nrPerClass, 1, 5), 
               rnorm(nrPerClass, i, 5)) 
  SampleData[,i] = SampleVector
  SampleImpact[i] = Impact(SampleVector, SampleClasses)$Impact
}
plot(SampleImpact, type='l', main=("Impact: \nIncreasing difference in means"))
}
\keyword{Impact}
