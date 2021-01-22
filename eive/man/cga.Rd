\name{cga}
\alias{cga}
\title{
Function for performing compact genetic algorithms
}
\description{
Function performs a compact genetic algorithm search for a given evaluation function.
}
\usage{
cga(chsize, popsize, evalFunc)
}
\arguments{
  \item{chsize}{
Number of bits
}
  \item{popsize}{
Number of population. By default it is 20
}
  \item{evalFunc}{
Function to minimize
}
}
\details{
CGA (Compact genetic algorithms) sample chromosomes using this probability vector.
A probability vector contains [P1,P2,...,PN] and the function generates and returns a chromosome [B1,B2,...,BN]. The probability
of BK having the value of 1 is PK. So, it has more chance to have [1,1,1,0,0] than [0,0,0,1,1] when the probability vector is
[0.9,0.9,0.9,0.1,0.1].
}
\value{
Returns the best chromosome with size of chsize.
}
\author{
Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>
Erkin Diyarbakirlioglu <ediyarbakirlioglu@gmail.com>
}


