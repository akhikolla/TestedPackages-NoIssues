\name{iBAT}
\alias{iBAT}

\title{
Main - Mixture selection prior
}
\description{
Perform MCMC iterations of the model, as described in the reference.
}
\usage{
iBAT(Y, X, distance, disfix, intercept=1, xi,
R=-1, tran, mu, sigma=((rgamma(4,1,1))^(-0.5)),
cmu=1/1000000, c=10, delta=3, d, e=0.001, f=0.999,
alpha=20, deltak=c(-1,0,0.58,1), tauk=c(1,1,1,2),
upp_bounds=c(-0.1, 0.1, 0.73, Inf),
low_bounds=c(-Inf, -0.1, 0.1, 0.73),
alpha_IG=c(1,1,1,1), beta_IG=c(1,1,1,1),
low_IG=c(0.41,0.41,0.41,1), a=c(1,1,1,1),
niter=500000, burnin=200000, Cout=1000,
phi=0.5, pR=0.4, selectioncgh=-1, pXI=0.6, indep=0)
}
\arguments{
  \item{Y}{
Matrix of gene expression data
}
  \item{X}{
Matrix of CGH data
}
  \item{distance}{
Vector of distance between CGH probes
}
  \item{disfix}{
Length of the chromosome under investigation
}
  \item{intercept}{
If set to one an intercept is included in the regression model
}
  \item{xi}{
Initialized matrix of latent states
}
  \item{R}{
Initialized association matrix in a vector form. Default set to -1, that automatically creates a vector with all the positions set to zero
}
  \item{tran}{
Initialized transition matrix
}
  \item{mu}{
Initialized state specific mean vector
}
  \item{sigma}{
Initialized state specific standard deviation vector
}
  \item{cmu}{
Parameter that controls the variance of the prior on the intercept
}
  \item{c}{
Parameter that determines the shrinkage in the model
}
  \item{delta}{
Parameter of the Inverse-Gamma prior on the error variance
}
  \item{d}{
Parameter of the Inverse-Gamma prior on the error variance
}
  \item{e}{
Parameter of the Beta prior on the inclusion probability
}
  \item{f}{
Parameter of the Beta prior on the inclusion probability
}
  \item{alpha}{
Parameter that regulates the strength of the independent part of the mixture
}
  \item{deltak}{
Vector of mean of the prior on the state specific mean
}
  \item{tauk}{
Vector of sd of the prior on the state specific mean
}
  \item{upp_bounds}{
Vector of upper bounds of the prior on the state specific mean
}
  \item{low_bounds}{
Vector of lower bounds of the prior on the state specific mean
}
  \item{alpha_IG}{
Vector of parameters of the prior on the state specific standard deviation
}
  \item{beta_IG}{
Vector of parameters of the prior on the state specific standard deviation
}
  \item{low_IG}{
Truncation of the prior on the state specific standard deviation
}
  \item{a}{
Vector of parameters of the prior on the transition matrix
}
  \item{niter}{
Number of Monte Carlo Markov Chain iterations
}
  \item{burnin}{
Burn-in
}
  \item{Cout}{
Print the number of iterations ran every Cout iterations
}
  \item{phi}{
Probability of an A/D step
}
  \item{pR}{
Parameter of the distribution used to select the rows to be updated at every MCMC iteration
}
  \item{selectioncgh}{
Number of samples not in neutral state in order to consider a CGH as a potential candidate for association with gene expression. Default set to -1 that automatically set it to 10\% of the samples
}
  \item{pXI}{
Parameter of the distribution used to select the rows to be updated at every MCMC iteration
}
  \item{indep}{
If set to an integer different from zero, run the analysis with an independent prior, see reference.
}
}
\value{
The output consists of an R list composed by 4*niter+3 elements objects, where niter is the number of MCMC iterations. The first niter objects of the list are vectors, each containing the positions of the association matrix set to one, at the corresponding MCMC iteration. Each of the following niter objects of the list are the transition matrices at the corresponding MCMC iteration, while the third and the fourth set of niter objects are the vectors of state specific mean and state specific variance, respectively. The last three objects of the list consist of three matrices counting the number of times the corresponding latent state has been set to 1,3 and 4, respectively.
}
\references{
Cassese A, Guindani M, Tadesse M, Falciani F, Vannucci M. A hierarchical Bayesian model for inference of copy number variants and their association to gene expression. Annals of Applied Statistics, 8(1), 148-175.
}
\author{
Alberto Cassese
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
\keyword{Main}
