\name{iBATCGH-package}
\alias{iBATCGH-package}
\alias{iBATCGH}
\docType{package}
\title{
Integrative Bayesian Analysis of Transcriptomic and CGH data
}
\description{
Bayesian integrative models of gene expression and comparative genomic hybridization data. The package provides inference on copy number variations and their association with gene expression
}
\details{
\tabular{ll}{
Package: \tab iBATCGH\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2015-07-06\cr
License: \tab GNU\cr
}
The package takes as inputs gene expression, Comparative Genomic Hybridization (CGH) data, and the physical distance between CGH probes. It returns posterior probabilities of inclusion for each potential association gene expression - CGH, and inference on the CGH latent states. Two different models have been considered, see references for details. Available functions are classified in four classes: Preprocessing, Main, Postprocessing, Simulated Data.
}
\author{
Alberto Cassese,
Marina Vannucci,
Michele Guindani,
Mahlet G. Tadesse.

Maintainer: Alberto Cassese <alberto.cassese@maastrichtuniversity.nl>
}
\references{
Cassese A, Guindani M, Tadesse M, Falciani F, Vannucci M. A hierarchical Bayesian model for inference of copy number variants and their association to gene expression. Annals of Applied Statistics, 8(1), 148-175.\cr
Cassese A, Guindani M, Vannucci M. A Bayesian integrative model for genetical genomics with spatially informed variable selection. Cancer Informatics.\cr
Guha S, Li Y, Neuberg D. Bayesian hidden Markov modelling of array cgh data. JASA. 2008;103(482):485-497.
}

\keyword{ Package }
\examples{
\dontrun{
data(NCI_60)

Y <- NCI_60$Affy
X <- NCI_60$aCGH
distance <- NCI_60$distance
disfix <- 146274826
xi <- InitXi(X)
tran <- Tran(xi)
mu <- InitMu()
d=0.2587288

Y <- Center(Y)

res <- iBAT(Y=Y,X=X,distance=distance,disfix=disfix,xi=xi,tran=tran,mu=mu,d=d)

summRes <- Inference(res,G=dim(Y)[[2]],M=dim(X)[[2]],niter=niter,burnin=bi,threshold=0.5)
}
}
