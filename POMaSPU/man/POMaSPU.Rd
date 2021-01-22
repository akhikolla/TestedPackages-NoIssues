\name{POMaSPU}
\alias{POMaSPU}
\title{
An adaptive sum of powered score (SPU) test for single SNP and mutiple phenotypes
}
\description{
POM-aSPU test evaluates an associations between an ordinal response (e.g. SNP) 
and multiple phenotypes;
provides a series of p-values of POM-SPU(pow) and POM-aSPU tests.
}
\usage{
POMaSPU(Y, Y.level, pheno, Z = NULL, pow = c(1:8, Inf), n.perm)
}
\arguments{
\item{Y}{a vector of categorical indicators}
\item{Y.level}{a vector indicating an order of categories either numeric or character e.g. Y.level = c(0, 1, 2) or 
	  Y.level = c("Well", "Mild", "Impaired").}
\item{pheno}{a numeric phenotype matrix with each row as a different individual and each column as a separate phenotype.}
\item{Z}{a numeric covariate matrix with each row as a different individual and each column as a covariate to be adjusted.}
\item{pow}{a vector of the power weight to be used at a trait level (default = c(1:8, Inf)).}
\item{n.perm}{a numeric value of number of null statistics (default = 1000).}
}
\details{
Adaptive association tests for an ordinal variable (e.g. SNP) and multiple phenotypes using POM.
}
\value{a vector of p-values from POM-SPU(pow) tests and POM-aSPU test.
}
\references{
Kim J, and Pan W (2017). "Adaptive testing for multiple traits in a proportional odds model with applications to detect SNP-brain network associations." \emph{Genetic Epidemiology}, 41(3), 259--277.
}
\author{Junghi Kim and Wei Pan
}
\note{
We provide a permutation based POM-aSPU test where null statistics are generated from permutations.
}

\examples{
# Example 1) 
# -- simulating phenotypes
# -- n.subjects: number of subjects
# -- n.traits: number of phenotypes

set.seed(136)
n.subjects <- 100
n.traits <- 10
pheno <- matrix(rnorm(n.subjects*n.traits), n.subjects, n.traits)


# -- simulating genotype (Y) and covariates (Z)
Y <- sample(c(0,1,2), n.subjects, replace = TRUE)
Z <- matrix(rnorm(n.subjects*2, 2, 1), n.subjects, 2)

# -- Computing the p-value of POMaSPU test with the permutation based method
Pvl <- POMaSPU(Y = Y, Y.level = c(0,1,2), pheno = pheno, Z = Z,  pow = c(1,2,4,Inf), n.perm = 1000)


# -- Each element of Pvl is a p value of POM-SPU(pow) in order
# -- The last element of Pvl is a p value of POM-aSPU test
Pvl
Pvl[length(Pvl)]


# > Pvl
#     SPU.1     SPU.2     SPU.4   SPU.Inf      aSPU 
# 0.1800000 0.5010000 0.3220000 0.1720000 0.3506494 


# > Pvl[length(Pvl)]
#      aSPU
# 0.3506494



# Example 2) 
# n.subjects <- 100
# n.traits <- 10
# Y <- sample(c("Well", "Mild", "Sick"), n.subjects, replace = TRUE)
# pheno <- matrix(rnorm(n.subjects*n.traits), n.subjects, n.traits)
# Z <- matrix(rnorm(n.subjects*2, 2, 1), n.subjects, 2)
# Pvl <- POMaSPU(Y = Y, Y.level = c("Sick", "Mild",  "Well"), pheno = pheno, Z = Z,  
# 		 pow = c(1,2,4,Inf), n.perm = 1000)
}
\keyword{ SNP based}
\keyword{ multiple phenotypes}
\keyword{ Proportional odds model}
