\name{patterns}
\alias{patterns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unique phyletic patterns with counts.
}
\description{
Generate a reduced dataset consisting of unique phyletic patterns
}
\usage{
patterns(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix or data frame of phyletic patterns. Rows represent the discrete character patterns and columns represent the taxa. These data can be numeric or character. }
}
\value{
\item{w        }{Counts of unique phyletic patterns.}
\item{databp_red        }{Reduced data set corresponding to w.}
\item{names        }{Column names of x.}
}
\author{
Utkarsh J. Dang and G. Brian Golding

\email{udang@mcmaster.ca}}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{estimaterates}}.
}
\examples{
data(simdata1)
uniq_patt <- patterns(simdata1$data)
print(uniq_patt)
}