\name{auxKnapsack01dp}
\alias{auxKnapsack01dp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multithreaded binary knapsack problem solver via dynamic programming
}
\description{
Given items' weights and values, concurrently solve 0-1 knapsack problems to optimality via dynamic programming for multiple knapsacks of different capacities.
}
\usage{
auxKnapsack01dp(
  weight,
  value,
  caps,
  maxCore = 7L,
  tlimit = 60,
  simplify = TRUE
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{weight}{
An integer vector.
}
  \item{value}{
A numeric vector. The size equals that of \code{weight}.
}
  \item{caps}{
An integer vector of knapsack capacities.
}
  \item{maxCore}{
Maximal threads to invoke. No greater than the number of logical CPUs on machine.
}
   \item{tlimit}{
Return the best exsisting solution in \code{tlimit} seconds.
}
  \item{simplify}{
If \code{length(caps) == 1}, simplify the output.
}
}
\details{
Implementation highlights include (i) lookup matrix is only of space complexity O(N * [max(C) - min(W)]), where N = the number of items, max(C) = maximal knapsack capacity, min(W) = minimum item weight; (ii) threads read and write the same lookup matrix and thus accelerate each other; (iii) the return of existing best solutions in time.
}
\value{
A list of 3:

\code{maxValue}: a numeric vector. \code{maxValue[i]} equals the sum of values of items selected for capacity \code{caps[i]}.

\code{selection}: a list of integer vectors. \code{selection[i]} indexes the items selected for capacity \code{caps[i]}.

\code{lookupTable}: a numeric matrix.
}


\note{
The function is not to solve the 0-1 multiple knapsack problem.
\code{weight} and \code{caps} are integers.
Be cautioned that dynamic programming is not suitable for problems with weights or capacities of high magnitudes due to its space complexity. Otherwise it could outperform branch-and-bound especially for large instances with highly correlated item weights and values.
}


\examples{
set.seed(42)
weight = sample(10L : 100L, 600L, replace = TRUE) # Dynamic programming
                                                  # solution requires integer
                                                  # weights.
value = weight ^ 0.5 * 100 # Higher correlation between item weights and values
                           # typically implies a harder knapsack problem.
caps = as.integer(runif(10, min(weight), 600L))
system.time({rstDp = FLSSS::auxKnapsack01dp(
  weight, value, caps, maxCore = 2, tlimit = 4)})
system.time({rstBb = FLSSS::auxKnapsack01bb(
  weight, value, caps, maxCore = 2, tlimit = 4)})
# Dynamic programming can be faster than branch-and-bound for integer weights
# and capacity of small magnitudes.
}

