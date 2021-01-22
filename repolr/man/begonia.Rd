\name{begonia}
\alias{begonia}
\docType{data}
\title{Begonia Pot Plant Quality Scores}
\description{
Begonia pot plant quality scores (for 2 varieties and 3 transport chains), during 5 weeks in simulated shelf-life conditions (temperature, light and irrigation). Quality scores were originally made on an ordinal scale from 1 to 10 (highest quality). However, only categories 2 to 9 were used, so these were re-coded to scale from 1 to 8. In addition to overall quality scores, a range of plant physiological characteristics were also observed. 
}
\usage{data(begonia)}
\format{
A data frame with 720 observations on the following variables.
  \describe{
    \item{\code{Pot}}{pot plant identifier; 1 to 144.}
    \item{\code{Plot}}{location in growing compartment; 1 to 48.}
    \item{\code{Week}}{week number in simulated shelf-life; 1 to 5.}
    \item{\code{Temp}}{temperature in compartment; \code{16C} or \code{21C}.}
    \item{\code{Light}}{light level in compartment; \code{High} or \code{Low}.}
    \item{\code{Chain}}{transport chain; \code{Cold}, \code{Comm} (commercial) or \code{Optm} (optimum).}
    \item{\code{Irrig}}{irrigation in compartment; \code{Cont} (control) or \code{Fluct} (fluctuating).}    
    \item{\code{Variety}}{variety; \code{Balli} or \code{Batik}.}
    \item{\code{Qual}}{quality score; 1 to 8.}
    \item{\code{FDrop}}{count of dropped flowers.}
    \item{\code{CBDrop}}{count of dropped coloured buds.}
    \item{\code{GBDrop}}{count of dropped green buds.}
    \item{\code{FDam}}{count of damaged flowers.}    
    \item{\code{FSing}}{count of single flowers.}
    \item{\code{FDoub}}{count of double flowers.}
    \item{\code{LYell}}{count of yellow leaves.}
  }
}
\source{
Parsons NR, Edmondson RN, Gilmour SG. A generalized estimating equation method for fitting autocorrelated ordinal score data with an application in horticultural research. \emph{Journal of the Royal Statistical Society C} 2006; 55:507-524.\cr
}
\examples{
data(begonia)
head(begonia)
}
\keyword{datasets}
