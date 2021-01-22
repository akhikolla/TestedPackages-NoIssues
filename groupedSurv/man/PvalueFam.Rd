\name{PvalueFam}
\alias{PvalueFam}
\title{
	Estimate P-values of the Efficient Score Statistics for the Frailty Model
}
\description{
  A method to estimate the p-value for the efficient score statistics for the frailty model.
}

\usage{
  PvalueFam(U)  	
}

\arguments{
  \item{U}{\code{Matrix} of contributions of each family to the efficient score statistics of each variable of interest, (i.e., the output from \code{groupedSurvFam}).}
}

\value{
  A \code{data.frame} with one row for each variables of interest and four columns: \code{stat}, the efficient score statistic, \code{pvalue}, unadjusted p-value, \code{FWER}, the family-wise error rate adjusted p-value (Bonferroni correction), and \code{FDR}, the local false-discovery rate (Storey's Q-value).
}

\references{
  Bonferroni, C.E. (1935). Il calcolo delle assicurazioni su gruppi di teste. Studi in Onore del Professore Salvatore Ortu Carbon, 13-60.\cr
  Storey, J.D. (2003). The positive false discovery rate: a Bayesian interpretation and the q-value. Annals of Statistics, 31:6, 2013-2035.\cr
  Storey, J.D., Taylor, J.E., and Siegmund, D. (2004). Strong control, conservative point estimation and simultaneous conservative consistency of false discovery rates: a unified approach. Journal of the Royal Statistical Society Series B-Statistical Methodology, 66, 187-205\cr
}


\examples{
# Generate dummy data	
U <- matrix(rnorm(10000), ncol=50)

#res   <- PvalueFam(U) 
#res
}

