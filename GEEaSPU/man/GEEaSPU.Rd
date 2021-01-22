\name{GEEaSPU}
\alias{GEEaSPU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An adaptive sum of powered score (SPU) test for single SNP and mutiple phenotypes
}
\description{
Tests whether a SNP is associated with multivariate outcomes;
provides a series of p-values of GEE-SPU(pow) and GEE-aSPU tests.
}
\usage{
GEEaSPU(pheno, geno, Z = NULL, model = "gaussian", corstr = "independence", 
        pow = c(1:8, Inf), n.perm = 1000, null.type = "perm", score.test = FALSE)
}
\arguments{
\item{pheno}{a numeric phenotype matrix with each row as a different individual and each column as a separate phenotype.}
\item{geno}{a numeric  vector with each element for an individual genotype.}
\item{Z}{a numeric covariate matrix with each row as a different individual and each column as a covariated to be adjested.}
\item{model}{a character string specifying the model of the phenotypes. 
Models supported are "gaussian" for a quantitative trait and "binomial" for a binary trait (default = "gaussian").
}
\item{corstr}{a character string specifying the correlation structure of phenotypes. The following are permitted: "independence", "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable", "AR-M" and "unstructured" (default = "independence").
}
\item{pow}{a vector of the power weight to be used at a trait level (default = c(1:8, Inf)).}
\item{n.perm}{a numeric value of number of null statistics (default = 1000).}
\item{null.type}{a character string specifying how to generate null statistics; "perm" is used when null statistics are generated using permutations and  
"sim" is used when null statistics are generated using simulations (default = "perm").}
\item{score.test}{a logical value indicating whether to include GEEaSPU-Score test along with GEE-Score test (default = FALSE). If TRUE, it computes p-values of GEEaSPU-Score and GEE-Score as well as GEEaSPU test.}
}
\details{
Adaptive association tests for single SNP and multiple phenotypes using GEE.
}
\value{a vector of p-values from GEE-SPU(pow) tests and GEE-aSPU test.
}
\references{
Kim J, Zhang Y, and Pan W (2016). "Powerful and adaptive testing for multi-trait and multi-SNP associations with GWAS and sequencing data." \emph{Genetics} 203(2), 715--731.


Zhang Y, Xu Z, Shen, X and  Pan W (2014). "Testing for association with multiple traits in generalized estimation equations, with application to neuroimaging data." \emph{NeuroImage}, 96, 309--325.
}
\author{Junghi Kim and Wei Pan
}
\note{
When large SNP-set (namely large gene size) or large number of phenotypes are included, the permuation based test (null.type = "perm") is recommended.

An option "binomial" model only supports the option, null.type = "sim".
}

\examples{
# -- simulating phenotypes
# -- n.subjects: number of subjects
# -- n.traits: number of phenotypes
# -- Sigma: covariance matrix of phenotypes (e.g. AR(1))

set.seed(136)
n.subjects <- 100
n.traits <- 3
sigma <- 2; rho <- 0.5
Sigma0 <- diag(n.traits); 
Sigma <- sigma * rho^abs(row(Sigma0) - col(Sigma0))
eS <- eigen(Sigma, symmetric = TRUE)
ev <- eS$values
X <- matrix(rnorm(n.traits * n.subjects), n.subjects)
pheno <- X \%*\% diag(sqrt(pmax(ev, 0)), ncol(Sigma)) \%*\% eS$vectors


# -- simulating genotype
geno <- rbinom(n = n.subjects, size = 2, prob = 0.2)


# -- Computing the p-value of GEEaSPU test with the permutation based method
Pvl <- GEEaSPU(pheno = pheno, geno = geno, Z = NULL, pow = c(1,2,4,Inf), 
        n.perm = 1000, null.type = "perm", score.test = FALSE)


# -- Each element of Pvl is a p value of GEE-SPU(pow) in order
# -- The last element of Pvl is a p value of GEE-aSPU test
Pvl
Pvl[length(Pvl)]


# > Pvl
#     SPU.1     SPU.2     SPU.4   SPU.Inf      aSPU 
# 0.1890000 0.4070000 0.3520000 0.3040000 0.2917083 


# > Pvl[length(Pvl)]
#      aSPU 
# 0.2917083
}
\keyword{ SNP based}
\keyword{ multiple phenotypes}
