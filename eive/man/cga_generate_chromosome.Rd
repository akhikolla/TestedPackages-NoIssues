\name{cga_generate_chromosome}
\alias{cga_generate_chromosome}
\title{
cga_generate_chromosome
}
\description{
Generates vector of zeros and ones for a given probability vector. 
}
\usage{
cga_generate_chromosome(prob_vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prob_vec}{
Vector of probabilities. 
}
}
\details{
This function is not directly called by user. CGA (Compact genetic algorithms) sample chromosomes using this probability vector.
A probability vector contains [P1,P2,...,PN] and the function generates and returns a chromosome [B1,B2,...,BN]. The probability
of BK having the value of 1 is PK. So, it has more chance to have [1,1,1,0,0] than [0,0,0,1,1] when the probability vector is
[0.9,0.9,0.9,0.1,0.1]. 
}
\value{
Returns the generated chromosome for a given probability vector. Return type is vector.
}
\author{
Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>
Erkin Diyarbakirlioglu <ediyarbakirlioglu@gmail.com>
}


