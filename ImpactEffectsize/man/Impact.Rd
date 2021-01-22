\name{Impact}
\alias{Impact}
\title{Impact effect size measure}
\description{
  Calculates the Impact effect size measure that is based on the group distance and the difference in the shape of the data distribution between two groups.
}
\usage{
Impact(Data, Cls, PlotIt = FALSE, pde = TRUE, 
  col = c("red","blue"), meanLines = FALSE, medianLines = FALSE, ...)
}
\arguments{
 \item{Data}{the data of both groups as a vector.}
 \item{Cls}{the class information as a vector of similar length as the data vector.}
 \item{PlotIt}{plots the pdf of the two groups as Pareto density estimation.}
 \item{col}{colors of the two lines to plot.
   Will be ignored if \code{PlotIt = FALSE}.}
 \item{pde}{plots a standard pdf insted of a non-parametric Pareto density estimation.
   Will be ignored if \code{PlotIt = FALSE}.}
 \item{meanLines}{adds perpendicular lines at the means of the two groups.
   Will be ignored if \code{PlotIt = FALSE}.}
 \item{medianLines}{adds perpendicular lines at the medians of the two groups.
   Will be ignored if \code{PlotIt = FALSE}.}
 \item{\dots}{more graphical parameters can be given as plot arguments if PlotIt = TRUE.}
}
\value{
  Returns a list of value incuding the main result and intermediate results.
 \item{Impact}{the main effect size measure.}
 \item{MorphDiff}{the extend of the group-difference in the shapes of the pdf.}
 \item{CTDiff}{the extend of the difference in the group medians.}
 }
\references{
 Lotsch, J., and Ultsch, A. (2019): 
 ImpactEffectsize -- an R Package for calculation and visualisation of the Impact distance and shape based effect size measure.
}
\author{
  Jorn Lotsch and Alfred Ultsch
}
\examples{
## example 1
data("FeatureselectionData")
ImpactSize <- Impact(FeatureselectionData$Var0011, FeatureselectionData$Classes)$Impact
plot2Densities(Data = FeatureselectionData$Var0011, Cls = FeatureselectionData$Classes, 
               xlab = "Variable", ylab = "Density") 
## example 2
data("BcellLymphomaCD79")
data("FeatureselectionData")
data("FlowcytometricData")
data("SameMeansData")
data("StocksFluctuation")
}
\keyword{Impact}
\keyword{effect size}
