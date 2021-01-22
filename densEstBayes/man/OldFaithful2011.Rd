\name{OldFaithful2011}
\alias{OldFaithful2011}
\usage{data(OldFaithful2011)}
\title{Intervals between geyser eruptions}
\description{
The \code{OldFaithful2011} array has time interval between eruptions (minutes) for all 3,507 eruptions of the Old Faithful Geyser (Yellowstone National Park, U.S.A.) in the year 2011.
}
\format{
Each entry of the array is the time interval in minutes between each geyser eruption in chronological order.
}
\source{The Geyser Observation and Study Association (www.geyserstudy.org).}
\examples{
library(densEstBayes) ; data(OldFaithful2011) 
hist(OldFaithful2011,breaks = 25,col = "gold",
     xlab = "time interval between geyser eruptions (minutes)")
dest <- densEstBayes(OldFaithful2011,method="SMFVB")
plot(dest,xlab = "time interval between geyser eruptions (minutes)")
rug(OldFaithful2011,col = "dodgerblue",quiet = TRUE)
}

