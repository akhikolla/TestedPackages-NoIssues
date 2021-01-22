\name{GEEaSPUset}
\alias{GEEaSPUset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An adaptive sum of powered score (SPU) test for SNP-set (gene) and mutiple phenotypes
}
\description{
Tests whether SNP-set (gene) is associated with multivariate outcomes;
provides a series of p-values of GEE-SPU(pow, pow2) and GEEaSPUset tests.
}
\usage{
GEEaSPUset(pheno, geno, Z = NULL, model = "gaussian", corstr = "independence", 
          pow = c(1,2,4,8), pow2 = c(1,2,4,8), n.perm = 1000, 
          null.type = "perm", score.test = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{pheno}{a numeric phenotype matrix with each row as a different individual and each column as a separate phenotype.}
\item{geno}{a numeric genotype matrix with each row as a different individual and each column as a snp.}
\item{Z}{a numeric covariate matrix with each row as a different individual and each column as a covariated to be adjested.}
\item{model}{a character string specifying the model of the phenotypes. 
Models supported are "gaussian" for a quantitative trait and "binomial" for a binary trait (default = "gaussian").
}
\item{corstr}{a character string specifying the correlation structure of phenotypes. The following are permitted: "independence", "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable", "AR-M" and "unstructured" (default = "independence").
}
\item{pow}{a vector of the power weight to be used at a SNP level (default = c(1,2,4,8)).}
\item{pow2}{a vector of the power weight to be used at a trait level (default = c(1,2,4,8)).}
\item{n.perm}{a numeric value of number of null statistics (default = 1000).}
\item{null.type}{a character string specifying how to generate null statistics; "perm" is used when null statistics are generated using permutations and  
"sim" is used when null statistics are generated using simulations (default = "perm").}
\item{score.test}{a logical value indicating whether to include GEEaSPU-Score test along with GEE-Score test (default = FALSE). If TRUE, it computes p-values of GEEaSPU-Score and GEE-Score as well as GEEaSPU test.}
}
\details{
Adaptive association tests for SNP-set (gene) and multiple phenotypes using GEE.
}
\value{a vector of p-values from GEE-SPU(pow, pow2) tests and GEE-aSPUset test.
}
\references{
Kim J, Zhang Y, and Pan W (2016). "Powerful and adaptive testing for multi-trait and multi-SNP associations with GWAS and sequencing data." \emph{Genetics}, 203(2), 715--731.
}
\author{
Junghi Kim and Wei Pan
}
\note{
When large SNP-set (namely large gene size) or large number of phenotypes are included, 
the permuation based test (null.type = "perm") is recommended.

An option "binomial" model only supports the option, null.type="sim".
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
X <- matrix(rnorm(n.subjects * n.traits), n.subjects)
pheno <- X \%*\% diag(sqrt(pmax(ev, 0)), ncol(Sigma)) \%*\% eS$vectors


# -- simulating genotype
# -- n.geno: number of SNPs included in the SNP set/gene

n.geno <- 3
maf <- 0.2
geno <- matrix(rbinom(n = n.subjects * n.geno, size = 2, prob = maf), ncol = n.geno)


# -- Computing the p-value of GEEaSPUpath test with the permutation based method
Pvl <- GEEaSPUset(pheno = pheno, geno = geno, Z = NULL, model = "gaussian", 
        corstr = "independence", pow = c(1,4,Inf), pow2 = c(1,4,Inf), n.perm = 1000,         
        null.type = "perm", score.test = FALSE)


# -- Each element of Pvl is a p value of GEE-SPU(pow,pow2) in order
# -- The last element of Pvl is a p value of GEE-aSPUset test
Pvl
Pvl[length(Pvl)]


# > Pvl
#     SPU.1.1     SPU.1.4   SPU.1.Inf     SPU.4.1     SPU.4.4   SPU.4.Inf   SPU.Inf.1 
#  0.01400000  0.08800000  0.07200000  0.53000000  0.41000000  0.32100000  0.55100000 
#   SPU.Inf.4 SPU.Inf.Inf     aSPUset 
#  0.48700000  0.41000000  0.04095904 
 
 
# > Pvl[length(Pvl)]
#    aSPUset 
# 0.04095904
}
\keyword{ SNP-set/gene based}
\keyword{ multiple phenotypes}
