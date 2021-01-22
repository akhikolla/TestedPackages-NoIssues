\name{control}
\alias{control}


\title{Set control parameters}

\description{
Set control parameters of the coordinate descent algorithm for the graphical lasso for sparse covariance matrix estimation.
}

\usage{
control(iter.out = 1000, iter.in = 100, tol.out = 1e-03, tol.in = 1e-02)
}


\arguments{
\item{iter.out}{Maximum number of iterations in the in the outer loop of the coordinate descent algorithm.}

\item{iter.in}{Maximum number of iterations in the in the inner loop of the coordinate descent algorithm.}

\item{tol.out}{Tolerance value for judging when convergence has been reached. Used in the outer loop of the coordinate descent algorithm.}

\item{tol.in}{Tolerance value for judging when convergence has been reached. Used in the inner loop of the coordinate descent algorithm.}

}


\details{
Function \code{control} is used to set control parameters of the coordinate descent algorithm employed for solving the covariance graphical lasso.
}

\value{
A list of parameters values.
}



\references{
Wang, H. (2014). Coordinate descent algorithm for covariance graphical lasso. \emph{Statistics and Computing}, 24:521.
}
