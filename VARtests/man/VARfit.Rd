\name{VARfit}
\alias{VARfit}

\title{
VAR(p) (Vector Autoregression) Model Fitting.
}
\description{
Estimates a VAR(p) model from data. The function's purpose is to be used before the \code{\link{ACtest}} or \code{\link{archBootTest}} functions, and doesn't have the same full functionality as for example the function \code{\link[vars:VAR]{VAR}} in the \code{vars} package.
}
\usage{
VARfit(y, p = 1, const = TRUE, trend = FALSE, exogen = NULL, univariate = FALSE)
}

\arguments{
  \item{y}{
an N x K matrix containing the time series to be modeled.
}
  \item{p}{
the lag length of the VAR(p) process.
}
  \item{const}{
if \code{TRUE} (default) the model will include a constant.
}
  \item{trend}{
if \code{TRUE} the model will include a trend.
}
  \item{exogen}{
a matrix or vector of exogenous variable(s). Must have the same number of rows as 'y' does.
}
  \item{univariate}{
if \code{TRUE} univariate (i.e. an AR(p) model for each timeseries) estimates will also be printed.
}
}



\value{
a list of class \code{"VARfit"} with the following slots: 

\item{y}{
        the data matrix.
    }
\item{p}{
        the lag length \code{p}.
    }
\item{N}{
        the the number of observations/rows in the time series 'y'.
    }
\item{K}{
        the the number of variables in the time series 'y'.
    }
\item{const}{
        \code{TRUE} if the model included a constant.
    }
\item{trend}{
        \code{TRUE} if the model included a trend.
    }
\item{exogen}{
        the exogen variable(s).
    }
\item{Z}{
        the [N - p] x [Kp + numberOf(const, trend, dummy)] design matrix, where the lags of the \code{y} matrix are appended.
    }
\item{call}{
        how the function \code{VARfit()} was called.
    }
\item{coef}{
        a matrix of the estimated parameters.
    }
\item{resid}{
        the [N - p] x K matrix of residuals.
    }
\item{uniCoef}{
        a matrix of the estimated parameter for the univariate case (an AR(p) for each column in 'y').
    }
\item{uniResid}{
        the [N - p] x K matrix of residuals from the univariate (AR(p)) cases.
    }
\item{univariate}{
        the 'univariate' argument.
    }
\item{NnonLagVar}{
        the number of non-endogenous variables (i.e. the number of 'exogen' variables + const + trend).
    }
\item{description}{
        who fitted the model and when.
    }
\item{time}{
        computation time taken for the estimation.
    }

}

\seealso{
\code{\link{ACtest}} to test for error autocorrelations, and \code{\link{archBootTest}} to test for ARCH errors.
}
\examples{

fit <- VARfit(y = VodafoneCDS, p = 3, const = TRUE, trend = FALSE)

}
