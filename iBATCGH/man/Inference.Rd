\name{Inference}
\alias{Inference}
\title{
Postprocessing - Posterior Inference
}
\description{
Performs posterior inference on the output of the main function
}
\usage{
Inference(listComplete, G, M, niter, burnin, threshold = 0.5)
}
\arguments{
  \item{listComplete}{
Output of the main function
}
  \item{G}{
Number of gene expression probes
}
  \item{M}{
Number of CGH probes
}
  \item{niter}{
Number of Monte Carlo Markov Chain iterations
}
  \item{burnin}{
Burn-in
}
  \item{threshold}{
Threshold on the posterior probabilities of inclusion of the association matrix
}
}
\details{
Wrapper function which calls InferenceR, InferenceXi, InferenceA, InferenceMu, InferenceSd on the output of the main function.
}
\value{
A list made by the following items
\item{R }{Binary matrix of estimated association}
\item{Xi }{Matrix of estimated copy number states}
\item{A }{Estimated transition matrix}
\item{Mu}{Estimated vector of state specific means}
\item{Sd }{Estimated vector of state specific standard deviations}
}
\references{
Cassese A, Guindani M, Tadesse M, Falciani F, Vannucci M. A hierarchical Bayesian model for inference of copy number variants and their association to gene expression. Annals of Applied Statistics, 8(1), 148-175.\cr
Cassese A, Guindani M, Vannucci M. A Bayesian integrative model for genetical genomics with spatially informed variable selection. Cancer Informatics.
}
\author{
Alberto Cassese
}
\seealso{
See Also \code{\link{InferenceR}},\code{\link{InferenceXi}}
}
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
\keyword{Postprocessing}
