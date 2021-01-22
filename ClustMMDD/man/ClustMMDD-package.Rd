\name{ClustMMDD-package}
\alias{ClustMMDD-package}
\alias{ClustMMDD}
\docType{package}
\title{
  \code{ClustMMDD} : Clustering by Mixture Models for Discrete Data.
}
\description{
  \code{ClustMMDD} stands for "Clustering by Mixture Models for Discrete Data". This package deals with the two-fold problem of variable selection and model-based unsupervised classification in discrete settings. Variable selection and classification are simultaneously solved via a model selection procedure using penalized criteria: Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), Integrated Completed Log-likelihood (ICL) or a general criterion with penalty function to be data-driven calibrated.
}
\details{
  \tabular{ll}{
  Package: \tab ClustMMDD\cr
  Type: \tab Package\cr
  Version: \tab 1.0.1\cr
  Date: \tab 2015-05-18\cr
  License: \tab GPL (>= 2)\cr
  }

  In this package, \code{K} and \code{S} are respectively the number of clusters and the subset of variables that are relevant for clustering purposes. We assume that a clustering variable has different probability distributions in at least two clusters, and a non-clustering variable has the same distribution in all clusters. We consider a general situation with data described by \eqn{P} random variables \eqn{X^l}, \eqn{l=1,\cdots,P}, where each variable \eqn{X^l} is an unordered set \eqn{\left\{X^{l,1},\cdots,X^{l,ploidy}\right\}} of \eqn{ploidy} categorical variables. For all \eqn{l}, the random variables \eqn{X^{l,1},\cdots,X^{l,ploidy}} take their values in the same set of levels. A typical example of such data comes from population genetics where each genotype of a diploid individual is constituted by \eqn{ploidy = 2} unordered alleles.
  
  The two-fold problem of clustering and variable selection is seen as a model selection problem. A specific collection of competing models associated to different values of \code{(K, S)} is defined, and are compared using penalized criteria. The penalized criteria are of the form \deqn{crit\left(K,S\right)=\gamma_n\left(K,S\right)+pen\left(K,S\right),} where 
  \itemize{
  \item \eqn{\gamma_n\left(K,S\right)} is the maximum log-likelihood,
  \item and \eqn{pen\left(K,S\right)} the penalty function.
  }
  
  The penalty functions used in this package are the following, where \eqn{dim\left(K,S\right)} is the dimension (number of free parameters) of the model defined by \eqn{\left(K,S\right)} :
  \itemize{
  \item{Akaike Information Criterion (AIC) : }{\deqn{pen\left(K,S\right) = dim\left(K,S\right)}}
  \item{Bayesian Information (BIC) : }{\deqn{pen\left(K,S\right) = 0.5*\log (n)*dim\left(K,S\right)}}
  \item{Integrated Complete Likelihood (ICL) : }{\deqn{pen\left(K,S\right) = 0.5*\log (n)*dim\left(K,S\right)+entropy\left(K,S\right),} where \deqn{entropy\left(K,S\right) = -\sum_{i=1}^N\sum_{k=1}^K\tau_{i,k}\log\left(\tau_{i,k}\right)} and \deqn{\tau_{i,k}=P\left(i\in\mathcal{C}_k\right)}.}
  \item{More general penalty function : }{\deqn{pen\left(K,S\right) = \alpha*\lambda*dim\left(K,S\right)} where  
  \itemize{
  \item \eqn{\lambda} is a multiplicative parameter to be calibrated,
  \item \eqn{\alpha} a coefficient in \eqn{[1.5,2]} to be given by the user.
  }
  
  We propose a data driven procedure based the dimension jumb version of the so called "slope heuristics" (see \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} and references therein).}
  }
  
  The maximum log-likelihood is estimated via the Expectation and Maximisation algorithm. The maximum a posteriori classification is derived from the estimated parameters of the selected model.
}

\author{
Wilson Toussile

Maintainer: Wilson Toussile <wilson.toussile@gmail.com>
}
\references{
\itemize{
		\item  \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} : Clustering and 	variable selection for categorical multivariate data. Electronic Journal of Statistics, Volume 7, 2344-2371, ISSN.
		\item 	\href{http://link.springer.com/article/10.1007\%2Fs11634-009-0043-x}{Wilson Toussile and Elisabeth Gassiat (2009)} : Variable 	selection in model-based clustering using multilocus genotype data. Adv Data Anal Classif, Vol 3, number 2, 109-134.
	}
}

\keyword{ Population genetics }
\keyword{ Clustering }
\keyword{ Unsupervized classification }
\keyword{ Mixture of multinomial models }
\keyword{ Expectation and Maximisation }
\keyword{ Model selection}
\keyword{ Slope-heuristics}
\keyword{ Dimension jump}

\seealso{
  The main functions :
  \describe{
	\item{\code{\link{em.cluster.R}}}{Compute an approximation of the maximum likelihood estimates of parameters using Expectation and Maximization algorithm, for a given value of \eqn{(K,S)}. The maximum a posteriori classification is then derived.}
	
	\item{\code{\link{backward.explorer}}}{Gather the most competitive models using a backward-stepwise strategy.}
	
	\item{\code{\link{dimJump.R}}}{Perform the data driven calibration of the penalty function via an estimation of \eqn{\lambda}. Two values are proposed and a graphic is proposed to help user in making a choice.}
	
	\item{\code{\link{selectK.R}}}{Perform the selection of the number \eqn{K} of clusters for a given subset of clustering variables.}
	
	\item{\code{\link{model.selection.R}}}{Perform a model selection from a collection of competing models.}
  }
}
\examples{
data(genotype2)
head(genotype2)
data(genotype2_ExploredModels)
head(genotype2_ExploredModels)

#Calibration of the penalty function
outDimJump = dimJump.R(genotype2_ExploredModels, N = 1000, h = 5, header = TRUE)
cte1 = outDimJump[[1]][1]
outSlection = model.selection.R(genotype2_ExploredModels, cte = cte1, header = TRUE)
outSlection
}
