\name{auxKnapsack01bb}
\alias{auxKnapsack01bb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multithreaded binary knapsack problem solver via branch and bound
}
\description{
Given items' weights and values, concurrently solve 0-1 knapsack problems to optimality via branch and bound for multiple knapsacks of different capacities.
}
\usage{
auxKnapsack01bb(
  weight,
  value,
  caps,
  itemNcaps = integer(0),
  maxCore = 7L,
  tlimit = 60,
  ub = "MT",
  simplify = TRUE
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{weight}{
A numeric vector.
}
  \item{value}{
A numeric vector. The size equals that of \code{weight}.
}
  \item{caps}{
A numeric vector of knapsack capacities.
}
  \item{itemNcaps}{
An integer vector of upper bounds on the number of selected items. \code{itemNcaps[i]} corresponds to instance \code{caps[i]}. Empty \code{itemNcaps} implies no size restriction.
}
  \item{maxCore}{
Maximal threads to invoke. No greater than the number of logical CPUs on machine.
}
  \item{tlimit}{
Return the best exsisting solution in \code{tlimit} seconds.
}
  \item{ub}{
Upper bound function.
}
  \item{simplify}{
If \code{length(caps) == 1}, simplify the output.
}
}
\details{
The algorithm takes the Horowitz-Sahni (1974) and the Martello-Toth (1977) upper bound functions and is carefully engineered towards speed. Implementation highlights include (i) an extra option of upper bounding the number of selected items, which only adds trivial overhead; (ii) the return of existing best solutions in time; (iii) the capability of taking numeric weights and values.}
\value{
A list of 2:

\code{maxValue}: a numeric vector. \code{maxValue[i]} equals the sum of values of items selected for capacity \code{caps[i]}.

\code{selection}: a list of integer vectors. \code{selection[i]} indexes the items selected for capacity \code{caps[i]}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
\note{
The function is not to solve the 0-1 multiple knapsack problem. The C++ implementation is fully independent and borrows no code from any open or commercial source.
}


\examples{
set.seed(42)
weight = runif(100, min = 1e3, max = 1e6)
value = weight ^ 0.5 * 100 # Higher correlation between item weights and values
                           # typically implies a harder knapsack problem.
caps = runif(10, min(weight), sum(weight))
rst = FLSSS::auxKnapsack01bb(weight, value, caps, maxCore = 2, tlimit = 2)
str(rst)
}




