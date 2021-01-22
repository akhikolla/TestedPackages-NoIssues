\name{wildBoot}

\alias{wildBoot}
\alias{print.wildBoot}

\title{
Wild Bootstrap Tests for Error Autocorrelation
}
\description{
This function runs Ahlgren and Catani's (2016) wild bootstrap (WB) test for error autocorrelation in VAR models. The function uses the test results from the function \code{\link{ACtest}} for the argument \code{'test'}, so that \code{\link{ACtest}} must be used before running \code{wildBoot}.
}
\usage{
wildBoot(test, WBtype = c("recursive", "fixed"), B = 199, 
         WBdist = c("rademacher", "normal", "mammen"), 
         HCtype = c("LM", "HC0", "HC1", "HC2", "HC3"), univariate = FALSE)
\method{print}{wildBoot}(x, ...)
}


\arguments{
  \item{test}{
An object of class \code{\link{ACtest}}.
}
  \item{WBtype}{
set this to \code{"recursive"} for the recursive version of the test, \code{"fixed"}. for the fixed version, or \code{c("recursive", "fixed")} to run both versions.
}
  \item{B}{
the number of bootstrap simulations (default is 199).
}
  \item{WBdist}{
the distribution used for the wild bootstrap error terms. Either \code{"rademacher"}, \code{"normal"} or \code{"mammen"}.
}
  \item{HCtype}{
a vector of the \code{HCtype}'s used for the WB test. The types must have been used in the \code{test} object.
}
  \item{univariate}{
either \code{FALSE}, \code{TRUE}, or \code{"only"}. If \code{TRUE}, the test will also be performed univariatly, if \code{"only"}, only the univariate version will be used. For \code{TRUE} and \code{"only"}, the univariate version must have been used in the \code{test} object.
}
  \item{x}{Object with class attribute \sQuote{wildBoot}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{

a list of class \code{"wildBoot"}. 

\item{test}{
        the \code{test} argument object.
}\item{WBtype}{
        the \code{WBtype} input.
}
\item{B}{
        the number of bootstrap simulations.
}
\item{WBdist}{
        the distribution used for the wild bootstrap error terms..
}
\item{HCtype}{
        the \code{HCtype}'s used for the WB test.
}
\item{description}{
        who ran the function and when.
}
\item{time}{
        computation time taken to run the WB test.
}
\item{call }{
        how the function \code{wildboot()} was called.
}
\item{numberOfErrors}{
        number of times the bootstrap simulations gave an error (can happen if one of the matrices becomes computational singular). If this happens, the particular bootstrap sample will be resimulated.
}
\item{numberOfNA}{
        number of times the test from the bootstrap simulations returned NA. If this happens, the particular bootstrap sample will be resimulated.
}
\item{WBr.Q}{
        a B x 5 matrix of the \bold{recursive} wild bootstrap (if used) Q statistics. Each row is the Q statistics from a bootstrap sample.
}
\item{WBr.pv}{
        the P. values from the \bold{recursive} design wild bootstrap (if used).
}
\item{WBf.Q}{
        a B x 5 matrix of the \bold{fixed} wild bootstrap (if used) Q statistics. Each row is the Q statistics from a bootstrap sample.
}
\item{WBf.pv}{
        the P. values from the \bold{fixed} design wild bootstrap (if used).
}
\item{uniList}{
        a list of the corresponding above items (\code{WBr.Q}, \code{WBf.Q}, \code{WBr.pv}, \code{WBf.pv}) for each respective univariate test. The elemtents of the list is named "uni1WBr.Q", "uni2WBr.pv", etc.
}

}
\references{
Ahlgren, N. & Catani, P. (2016).
\emph{Wild bootstrap tests for autocorrelation in vector autoregressive models.} Stat Papers, <doi:10.1007/s00362-016-0744-0>.
}

\seealso{
\code{\link{ACtest}} for the test that is required to be run before the \code{wildBoot} function.
}
\examples{

\dontrun{
fit <- VARfit(y = VodafoneCDS, p = 3, const = TRUE, trend = FALSE)
test <- ACtest(fit = fit, h = 1, HCtype = c("LM", "HC0", "HC1", "HC2", "HC3"))
wb <- wildBoot(test = test, WBtype = "recursive", B = 199, 
               WBdist = "rademacher",  HCtype = c("LM", "HC3"))
}
}
