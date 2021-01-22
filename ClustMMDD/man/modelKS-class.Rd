\name{modelKS-class}
\Rdversion{1.0.00}
\docType{class}

\alias{modelKS-class}
\alias{==,modelKS,ANY-method}
\alias{[<-,modelKS-method}
\alias{[,modelKS-method}
\alias{is.element,modelKS,modelKS-method}
\alias{show,modelKS-method}
\alias{simulData,modelKS,numeric,numeric-method}
\alias{read.modelKS,character-method}
\alias{is.modelKS,ANY-method}
\alias{model,modelKS-method}



\title{\code{modelKS} is a class of parameters of \eqn{(K, S)} model.}

\description{
  \code{modelKS} is a class that can contain the set of parameters associated to a model given by \eqn{(K, S)}.
}

\section{Objects from this class}{
  Objects can be created by calling \code{new("modelKS", ...)}. See \code{\link{new}} for more details. 
}

\section{Slots}{
  \describe{
    \item{\code{N}:}{The number of individuals in the datatset.}
    \item{\code{P}:}{The number of random variables considered in the dataset.}
    \item{\code{N_levels}:}{A vector of the numbers of levels for the considered variables.}
    \item{\code{levels}:}{A \code{"list"} of the observed levels for the variables.}
    \item{\code{K}:}{The number of clusters.}
    \item{\code{S}:}{A vector of \code{"logical"} indicating the clustering variables.}
    \item{\code{dim}:}{The dimension of a model \eqn{(K,S)} defined as the number of free parmaters.}
    \item{\code{mixingProportions}:}{The numeric vector of the mixing proportions.}
    \item{\code{count}:}{A \code{"list"} of the counts of levels for each variable.}
    \item{\code{frequencies}:}{A \code{"list"} of the observed frequencies for each variable.}
    \item{\code{proba}:}{A \code{"list"} of \code{"matrix"} that contains the estimates of the levels probabilities in each clusters.}
    \item{\code{logLik}:}{An approximation of the maximum log-likelihood obtained by the EM algorithm.}
    \item{\code{entropy}:}{The entropy given by \eqn{-\sum_{i=1}^N\sum_{k=1}^K\tau_{i,k}log\left(\tau_{i,k}\right)}, where \eqn{\tau_{i,k}} is the probability that individual \eqn{i} belongs to cluster \eqn{k}}
    \item{\code{membershipProba}:}{The \code{"numeric"} matrix of membership probabilities.}
    \item{\code{mapClassification}:}{The maximum a posteriori classification given by a vector of \code{"integers"}.}
  }
}

\section{Methods}{
  \describe{
    \item{==}{\code{signature(e1 = "modelKS", e2 = "ANY")}: ... }
    \item{[<-}{\code{signature(x = "modelKS")}: ... }
    \item{[}{\code{signature(x = "modelKS")}: ... }
    \item{is.element}{\code{signature(el = "modelKS", set = "modelKS")}: ... }
    \item{show}{\code{signature(object = "modelKS")}: ... }
    \item{simulData}{\code{signature(object = "modelKS", N = "numeric", ploidy = "numeric")}: ... }
    \item{read.modelKS}{\code{signature(file = "character")}: ... }
    \item{is.modelKS}{\code{signature(object = "modelKS")}: ... }
	 }
}

\author{
Wilson Toussile.
}

\references{
\itemize{
		\item  \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} : Clustering and variable selection for categorical multivariate data. Electronic Journal of Statistics, Volume 7, 2344-2371, ISSN.
		
		\item 	\href{http://link.springer.com/article/10.1007\%2Fs11634-009-0043-x}{Wilson Toussile and Elisabeth Gassiat (2009)} : Variable 	selection in model-based clustering using multilocus genotype data. Adv Data Anal Classif, Vol 3, number 2, 109-134.
	}
}

\seealso{
  \code{\link{slotNames}}, \code{\link{new}}, \code{\link{methods}}, \code{\link{showClass}}.
}

\keyword{classes}

\examples{
data(exModelKS)
showClass("modelKS")
slotNames("modelKS")
exModelKS
exModelKS["K"]
exModelKS["S"]
model(exModelKS)
}