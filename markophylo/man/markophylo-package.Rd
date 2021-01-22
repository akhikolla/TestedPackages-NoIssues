\name{markophylo-package}
\alias{markophylo-package}
\alias{markophylo}
\docType{package}
\title{
Markov chain models for phylogenetic trees.
}
\description{
Allows for fitting of maximum likelihood models using Markov chains on phylogenetic trees for analysis of discrete character data. Examples of such discrete character data include restriction sites, gene family presence/absence, intron presence/absence, and gene family size data. Hypothesis-driven user-specified substitution rate matrices can be estimated. Allows for biologically realistic models combining constrained substitution rate matrices, site rate variation, site partitioning, branch-specific rates, allowing for non-stationary prior root probabilities, correcting for sampling bias, etc. 
}
\details{
\tabular{ll}{
Package: \tab markophylo\cr
Type: \tab Package\cr
Version: \tab 1.0.5\cr
Date: \tab 2016-01-19\cr
License: \tab GPL (>= 2)\cr
}
Function \code{\link{estimaterates}} is the main workhorse of the package, allowing users to specify a substitution rate matrix and estimate it. Partition analysis, gamma rate variation, estimating root frequencies of the discrete characters, and specifying clades (or a set of branches) with their own unique rates are all possible.
}
\author{
Utkarsh J. Dang and G. Brian Golding

\email{udang@mcmaster.ca}
}
\references{
Dirk Eddelbuettel and Conrad Sanderson (2014). RcppArmadillo: Accelerating R with high-performance C++ linear algebra. Computational Statistics and Data Analysis, 71, 1054--1063.
  
Eddelbuettel, Dirk and Romain Francois (2011). Rcpp: Seamless R and C++ Integration. Journal of Statistical Software, 40(8), 1--18.
  
Gilbert, Paul and Ravi Varadhan (2012). numDeriv: Accurate Numerical Derivatives. R package version 2014.2-1.

Harmon Luke J, Jason T Weir, Chad D Brock, Richard E Glor, and Wendell Challenger. 2008. GEIGER: investigating evolutionary radiations. Bioinformatics 24:129-131.

Paradis, E., Claude, J. & Strimmer, K. 2004. APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290. R package version 3.2.
  
Schliep, K.P. 2011. phangorn: phylogenetic analysis in R. Bioinformatics, 27(4) 592-593. R package version 1.99.13.

Yang, Z. (1994). Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods. Journal of Molecular evolution, 39(3), 306--314.
}
\keyword{ package }
