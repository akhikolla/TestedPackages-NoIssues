\name{GEEaSPUpath}
\alias{GEEaSPUpath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An adaptive sum of powered score (SPU) test for gene-set (pathway) and mutiple phenotypes
}
\description{
Tests whether gene-set (pathway) is associated with multivariate outcomes;
provides a series of p-values of GEE-SPU(pow, pow2, pow3) and GEEaSPUpath tests.
}
\usage{
GEEaSPUpath(pheno, geno, nSNPs, Z = NULL, corstr = "independence", 
           pow = c(1,2,4,8), pow2 = c(1,2,4,8), pow3 = c(1,2,4,8), n.perm = 1000)}

\arguments{
\item{pheno}{a numeric phenotype matrix with each row as a different individual and each column as a separate phenotype.}
\item{geno}{a numeric genotype matrix with each row as a different individual and each column as a snp; the SNPs (with the number stored in nSNPs) from one gene are 
stored consecutively from the first gene.}
\item{nSNPs}{A numeric vector, whose length matches to the total number of genes; 
each element of vector indicate the number of SNPs in each gene.}
\item{Z}{a numeric covariate matrix with each row as a different individual and each column as a covariated to be adjested.}
\item{corstr}{a character string specifying the correlation structure of phenotypes. The following are permitted: "independence", "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable", "AR-M" and "unstructured" (default = "independence").
}
\item{pow}{a vector of the power weight to be used at a SNP level (default = c(1,2,4,8)).}
\item{pow2}{a vector of the power weight to be used at a trait level (default = c(1,2,4,8)).}
\item{pow3}{a vector of the power weight to be used at a gene level (default = c(1,2,4,8)).}
\item{n.perm}{a numeric value of number of null statistics (default = 1000).}
}
\details{
Adaptive association tests for gene-set (pathway) and multiple phenotypes using GEE.
}
\value{a vector of p-values from GEE-SPU(pow, pow2, pow3) tests and GEE-aSPUpath test.
}
\references{
Kim J, Zhang Y, and Pan W (2016). "Powerful and adaptive testing for multi-trait and multi-SNP associations with GWAS and sequencing data." \emph{Genetics}, 203(2), 715--731.
}
\author{Junghi Kim and Wei Pan
}
\note{
GEEaSPUpath function only supports a case for a quantitative trait (model = "gaussian") and a permutation based test (null.type = "perm"). 
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
Sigma0 <- diag(n.traits)
Sigma <- sigma * rho^abs(row(Sigma0) - col(Sigma0))
eS <- eigen(Sigma, symmetric = TRUE)
ev <- eS$values
X <- matrix(rnorm(n.subjects * n.traits), n.subjects)
pheno <- X \%*\% diag(sqrt(pmax(ev, 0)), ncol(Sigma)) \%*\% eS$vectors


# -- simulating genotype
# -- Assume we have two genes each of which has 3 and 5 SNPs respectively.
# -- n.geno1: number of SNPs included in the gene1
# -- n.geno2: number of SNPs included in the gene2
# -- nSNPs <- c(3,5)

n.geno1 <- 3
n.geno2 <- 5
maf1 <- 0.2
maf2 <- 0.4
gene1 <- matrix(rbinom(n = n.subjects*n.geno1, size = 2, prob = maf1), ncol = n.geno1)
gene2 <- matrix(rbinom(n = n.subjects*n.geno2, size = 2, prob = maf2), ncol = n.geno2)
geno <- cbind(gene1, gene2)


# -- Computing the p-value of GEEaSPUpath test
Pvl <- GEEaSPUpath(pheno = pheno, geno = geno, nSNPs = c(3,5), Z = NULL, 
        corstr = "independence", pow = c(1,4,8), pow2 = c(1,4,8), 
        pow3 = c(1,4,8), n.perm = 1000)


# -- Each element of Pvl is a p value of GEE-SPU(pow,pow2,pow3) in order
# -- The last element of Pvl is a p value of GEE-aSPUpath test
Pvl
Pvl[length(Pvl)]


# > Pvl
# SPU.1.1.1  SPU.1.1.4  SPU.1.1.8  SPU.1.4.1  SPU.1.4.4  SPU.1.4.8  SPU.1.8.1  SPU.1.8.4 
# 0.00900000 0.05600000 0.07000000 0.06200000 0.08300000 0.11200000 0.06100000 0.08200000 
#  SPU.1.8.8  SPU.4.1.1  SPU.4.1.4  SPU.4.1.8  SPU.4.4.1  SPU.4.4.4  SPU.4.4.8  SPU.4.8.1 
# 0.10600000 0.58100000 0.54300000 0.49200000 0.62400000 0.64000000 0.62700000 0.64900000 
#  SPU.4.8.4  SPU.4.8.8  SPU.8.1.1  SPU.8.1.4  SPU.8.1.8  SPU.8.4.1  SPU.8.4.4  SPU.8.4.8 
# 0.67100000 0.67500000 0.58300000 0.53700000 0.48100000 0.63400000 0.64600000 0.63800000 
#  SPU.8.8.1  SPU.8.8.4  SPU.8.8.8   aSPUpath 
# 0.66000000 0.68100000 0.67900000 0.04395604 


# > Pvl[length(Pvl)]
#   aSPUpath 
# 0.04395604 
}
\keyword{ gene-set/pathway based}
\keyword{ multiple phenotypes}
