\name{VARsim}
\alias{VARsim}

\title{
Simulates vector autoregressive (VAR) series
}
\description{
This function simulates VAR(p) series. For quick use, the function can use the estimated 
model returned by \code{\link{VARfit}} or \code{\link[vars:VAR]{VAR}} as
the DGP (Data Generating Process), if passed to the \code{fittedModel} argument.
}

\usage{
VARsim(N = 200, K = 2, p = 1, const = TRUE, trend = FALSE, exogen = NULL, 
       coef = NULL, dist = "normal", Ystart = NULL, errors = NULL, fittedModel = NULL)
}

\arguments{
  \item{N}{
  The length of the series.
}
  \item{K}{
  The number of series/equations.
}
  \item{p}{
  The lag length of the VAR(p).
}
  \item{const}{
  if \code{TRUE}, a constant will be used.
}
  \item{trend}{
  if \code{TRUE}, a trend will be used.
}
  \item{exogen}{
  a matrix of exogenous variables. Should either have \code{N} or (\code{N} - \code{p} number of rows. If the latter, the last row will be matched with the last row of the simulated series \code{y}.
}
  \item{coef}{
  a matrix of coefficients. 
  E.g. a 2 dimensional VAR(2) with a constant, a trend and one exogenous variable must be entered in the following order (the same as returned by \code{VARfit}):
  
  \ifelse{latex}{
  
\eqn{ 
\begin{array}{cc}
\bold{\mbox{\ \ \ \ \ \ y1}} & \bold{\mbox{\ \ \ \ \ \ y2}} \\ 
\end{array}%
}    

\eqn{
\left[  
\begin{array}{cc}
\mbox{const} & \mbox{const} \\ 
\mbox{trend} & \mbox{trend} \\ 
\mbox{exo1} & \mbox{exo1} \\ 
\mbox{y1[-1]} & \mbox{y1[-1]} \\ 
\mbox{y2[-1]} & \mbox{y2[-1]} \\ 
\mbox{y1[-2]} & \mbox{y1[-2]} \\ 
\mbox{y2[-2]} & \mbox{y2[-2]} \\ 
\end{array}%
\right]
}  


  }{

  \verb{
  
 ___y1_______y2___
| const  | const  |
| trend  | trend  |
| exo1   | exo1   |
| y1[-1] | y1[-1] |
| y2[-1] | y2[-1] |
| y1[-2] | y1[-2] |
| y2[-2] | y2[-2] |
 -----------------
}
}

}
  \item{dist}{
  the distribution of the error terms (currently only \code{"normal"} for i.i.d. standard normal.)
}
  \item{Ystart}{
  a \code{p} x \code{K} matrix of the start values. These will be present in the simulated series. If \code{NULL}, the start values will be set to zero.
}
  \item{errors}{
  (optional) a matrix of error terms. If supplied, \code{VARsim()} will use those instead of sampling the errors. Must have \code{K} columns and either \code{N} or (\code{N} - \code{p}) number of rows. If the latter, the last row will be matched with the last row of the simulated series \code{y}.
}
  \item{fittedModel}{
  (optional) either an object of class \code{"vars"}, as returned by \code{\link[vars:VAR]{VAR}}, or an object of class \code{"VARfit"}, as returned by \code{\link{VARfit}}. If used, the estimated model of the \code{fittedModel} will be used as a DGP (Data Generating Process). Unless the arguments \code{N}, \code{exogen} and/or \code{Ystart} were supplied, they too will be taken from the \code{fittedModel} object.
}
}

\value{
an \code{N} x \code{K} matrix of the simulated VAR(\code{p}).
}

\examples{

# Simulates from a fitted DGP:
fit <- VARfit(y = VodafoneCDS, p = 1, const = TRUE)
y <- VARsim(fittedModel = fit)
VARfit(y = y, p = 1, const = TRUE)

# Simulates from user given parameters. 
#  Includes an exogenous variable:
VARcoef <- matrix(c(1,     2,
                    1,   0.5,
                    0.1,  0.3,
                    0.2,  0.5), 
                  ncol = 2, byrow = TRUE)
exo <- matrix(rnorm(500))
y <- VARsim(N = 500, K = 2, p = 1, const = TRUE, 
            trend = FALSE, exogen = exo,
            coef = VARcoef, dist = "normal")
VARfit(y = y, p = 1, const = TRUE, exogen = exo)

\dontrun{


}}
