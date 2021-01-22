\name{QTL.gCIMapping.GUI-package}
\alias{QTL.gCIMapping.GUI-package}
\alias{QTL.gCIMapping.GUI}
\docType{package}
\title{QTL Genome-Wide Composite Interval Mapping with Graphical User Interface
}
\description{
Conduct multiple quantitative trait loci (QTL) mapping under the framework of random-QTL-effect mixed linear model. First, each position on the genome is detected in order to construct a negative logarithm P-value curve against genome position. Then, all the peaks on each effect (additive or dominant) curve are viewed as potential QTL, all the effects of the potential QTL are included in a multi-QTL model, their effects are estimated by empirical Bayes in doubled haploid or by adaptive lasso in F2, and true QTL are identified by likelihood radio test.}
\usage{
QTL.gCIMapping.GUI()
}
\details{
\tabular{ll}{
Package: \tab QTL.gCIMapping.GUI\cr
Type: \tab Package\cr
Version: \tab 2.1.1\cr
Date: \tab 2020-10-8\cr
Depends: \tab shiny,MASS,qtl\cr
Imports: \tab methods,openxlsx,stringr,Rcpp\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
}
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\references{
An efficient multi-locus mixed model framework for the detection of small and linked QTLs in F2.Wen Yang-Jun, Zhang Ya-Wen, Zhang Jin, Feng Jian-Ying, Jim M. Dunwell, Zhang Yuan-Ming*
}
\examples{
\dontrun{QTL.gCIMapping.GUI()}
}

