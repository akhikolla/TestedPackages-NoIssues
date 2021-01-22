\name{gjamFillMissingTimes}
\alias{gjamFillMissingTimes}
\title{Fill out data for time series (state-space) gjam}
\description{
  Fills in predictor, response, and effort matrics for time series data where there are multiple multivariate time series.  Time series gjam is still under development.
}
\usage{
  gjamFillMissingTimes(xdata, ydata, edata, groupCol, timeCol, groupVars = groupCol,
                                 FILLMEANS = FALSE, typeNames = NULL, missingEffort = .1)
                       
}
\arguments{
  \item{xdata}{\code{n by Q data.frame} holding predictor variables}
  
  \item{ydata}{\code{n by S matrix} holding response variables}
  
  \item{edata}{\code{n by S matrix} holding effort}
  
  \item{groupCol}{column name in \code{xdata} for group variable, i.e., observations part of the same time series}
  
  \item{timeCol}{column name in \code{xdata} for time index}
  
  \item{groupVars}{\code{character vector} of column names in \code{xdata} having values that are fixed for a value of \code{groupCol}, i.e., they do not change with time index in \code{timeCol}}
  
  \item{FILLMEANS}{fill new rows in \code{ydata} with mean for \code{groupCol} times \code{missingEffort}; otherwise NA}
  
  \item{typeNames}{typenames current limited to \code{'DA'} for discrete counts}
  
  \item{missingEffort}{effort assigned to missing values of \code{edata} and \code{ydata}}
}
\details{
Missing times in the data occur where there are gaps in \code{timeCol} column of \code{xdata} and the initial time \code{0} for each sequence.  New versions of the data have \code{NA} (\code{xdata}) or prior values with appropriate weight (\code{ydata}).  Missing times are filled in \code{xdata, ydata, edata}, including a time \code{0} which serves as a prior mean for \code{ydata} for time code{1}.  The group and time indices in columns \code{groupCol} and \code{timeCol} of \code{xdata} reference the time for a given time series.  Missing values in the columns \code{groupVars} of \code{xdata} are filled automatically filled in.  This assumes that values for these variables are fixed for the group.  If \code{FILLMEANS}, the missing  values in \code{ydata} are filled with means for the group and given a low weight specified in \code{missingEffort}.
}
\value{
A list containing the following:

  \item{xdata}{filled version of \code{xdata}}
  
  \item{ydata}{filled version of \code{ydata}}
  
  \item{edata}{filled version of \code{edata}}
  
  \item{timeList}{time indices used for computation, including, \code{timeZero} (row numbers in new data where each time series begins, with times = 0), \code{timeLast} (row numbers in new data where each time series ends), \code{rowInserts} (row numbers for all inserted rows), \code{noEffort} (rows for which effort in \code{edata} is filled with \code{missingEffort})}
 
}

\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2016.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs, 87, 34-56.
}
\seealso{
\code{\link{gjam}} for more on \code{xdata}, \code{ydata}, and effort.

A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

web site 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

