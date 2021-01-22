\name{predict.repolr}
\alias{predict.repolr}
\title{
Predict Method for Fitted repolr Model
}
\description{
Calculates predictions and standard errors of predictions for a fitted \code{repolr} model object.
}
\usage{
\method{predict}{repolr}(object, newdata = NULL, se.fit = FALSE,
           robust.var = TRUE, type = c("link", "response", "terms"), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{is a model fitted using \code{repolr}.}
  \item{newdata}{optionally, a data frame in which to find variables with which to predict; if missing the model fitted values are reported.}
  \item{se.fit}{Logical indicating if standard errors are required.}
  \item{robust.var}{logical; if \code{TRUE}, standard errors are based on the robust variance matrix.}  
  \item{type}{is the type of prediction required. The default \dQuote{\code{link}} is to use the scale of the linear predictors; i.e. the log-odds of cumulative probabilities. The alternative is to report the predicted cumulative probabilities; \dQuote{\code{response}}. The \dQuote{\code{terms}} option returns the matrix of fitted values for each model term on the scale of the linear predictor.}
  \item{\dots}{further arguments passed to or from other methods.}  
}
\details{
If newdata is missing predictions are based on the data used to fit the \code{repolr} model. If newdata are supplied then the format of these data must conform to the same format required for model fitting using \code{repolr}. See \code{\link{repolr}} for details.
}

\value{
 \item{fit}{Predictions.}
 \item{se.fit}{Estimated standard errors.}
}


\examples{
data(HHSpain)
mod.0 <- repolr(HHSpain~Sex*Time, data=HHSpain, categories=4, subjects="Patient",
            times=c(1,2,5), corr.mod="uniform", alpha=0.5)
predict(mod.0, newdata = data.frame(Patient = rep(100, 3), Time = c(1, 2, 5),
       Sex = factor(rep(1, 3), levels=1:2, labels=c("F", "M"))), 
       type="link", se.fit = TRUE)
}