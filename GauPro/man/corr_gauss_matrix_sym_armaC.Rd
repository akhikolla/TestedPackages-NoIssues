% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{corr_gauss_matrix_sym_armaC}
\alias{corr_gauss_matrix_sym_armaC}
\title{Correlation Gaussian matrix in C using Armadillo (symmetric)}
\usage{
corr_gauss_matrix_sym_armaC(x, theta)
}
\arguments{
\item{x}{Matrix x}

\item{theta}{Theta vector}
}
\value{
Correlation matrix
}
\description{
Correlation Gaussian matrix in C using Armadillo (symmetric)
}
\examples{
corr_gauss_matrix_sym_armaC(matrix(c(1,0,0,1),2,2),c(1,1))
}
