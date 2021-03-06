% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/specs.R
\name{specs_opt}
\alias{specs_opt}
\title{SPECS with data transformation and penalty optimization}
\usage{
specs_opt(
  y,
  x,
  p = 1,
  rule = c("BIC", "AIC", "TSCV"),
  CV_cutoff = 2/3,
  deterministics = c("constant", "trend", "both", "none"),
  ADL = FALSE,
  weights = c("ridge", "ols", "none"),
  k_delta = 1,
  k_pi = 1,
  lambda_g = NULL,
  lambda_i = NULL,
  thresh = 1e-04,
  max_iter_delta = 1e+05,
  max_iter_pi = 1e+05,
  max_iter_gamma = 1e+05
)
}
\arguments{
\item{y}{A vector containing the dependent variable in levels.}

\item{x}{A matrix containing the independent variables in levels.}

\item{p}{Integer indicating the desired number of lagged differences to include. Default is 1.}

\item{rule}{A charcater object indicating which selection rule the optimal choice of the penalty parameters is based on.
Default is "BIC".}

\item{CV_cutoff}{A numeric value between 0 and 1 that decides the proportion of the training sample as a fraction of the complete
sample. Applies only when rule="TSCV". Default is 2/3.}

\item{deterministics}{A character object indicating which deterministic variables should be added
("none","constant","trend","both"). Default is "constant".}

\item{ADL}{Logical object indicating whether an ADL model without error-correction term should be estimated.
Default is FALSE.}

\item{weights}{Choice of penalty weights. The weights can be automatically generated by ridge regression (default)
or ols. Alternatively, a conformable vector of non-negative weights can be supplied.}

\item{k_delta}{The power to which the weights for delta should be raised, if weights are set to "ridge" or "ols".}

\item{k_pi}{The power to which the weights for pi should be raised, if weights are set to "ridge" or "ols".}

\item{lambda_g}{An optional user-specified grid for the group penalty may be supplied. If left empty,
a 10-dimensional grid containing 0 as the minimum value is generated.}

\item{lambda_i}{An optional user-specified grid for the individual penalty may be supplied. If left empty,
a 10-dimensional grid containing 0 as the minimum value is generated.}

\item{thresh}{The treshold for convergence.}

\item{max_iter_delta}{Maximum number of updates for delta. Default is \eqn{10^5}.}

\item{max_iter_pi}{Maximum number of updates for pi. Default is \eqn{10^5}.}

\item{max_iter_gamma}{Maximum number of updates for gamma. Default is \eqn{10^5}.}
}
\value{
\item{D}{A matrix containing the deterministic variables included in the model.}
\item{gammas}{A matrix containing the estimated coefficients of
the stochastic variables in the conditional error-correction model.}
\item{gamma_opt}{A vector containing the estimated coefficients corresponding to the optimal model.}
\item{lambda_g}{The grid of group penalties.}
\item{lambda_i}{The grid of individual penalties.}
\item{Mv}{A matrix containing the independent variables, after regressing out the deterministic components.}
\item{My_d}{A vector containing the dependent variable, after regressing out the deterministic components.}
\item{theta}{The estimated coefficients for the constant and trend. If a deterministic component is excluded,
its coefficient is set to zero.}
\item{theta_opt}{The estimated coefficients for the constant and trend in the optimal model.}
\item{v}{A matrix containing the independent variables (excluding deterministic components).}
\item{weights}{The vector of penalty weights.}
\item{y_d}{A vector containing the dependent variable, i.e. the differences of y.}
}
\description{
This function estimates SPECS and selects the optimal penalty parameter based on a selection rule. All arguments correspond
to those of the function specs(...), but it contains the additional arguments rule and CV_cutoff.
Selection of the penalty parameter can be carried out by BIC or AIC or by time series cross-validation (TSCV). The degrees
of freedom for the information criteria (BIC or AIC) are approximated by the number of non-zero coefficients in the
estimated model. TSCV cuts the sample in two, based on the argument CV_cutoff which determines the proportion of
the training sample. SPECS is estimated on the first part and the estimated model is used to predict the values
in the second part. The selection is then based on the lowest Mean-Squared Forecast Error (MSFE) obtained over the test sample.
}
\examples{

#Estimate an automatically optimized model for unemployment and ten google trends

#Organize data
y <- Unempl_GT[,1]
index_GT <- sample(c(2:ncol(Unempl_GT)),10)
x <- Unempl_GT[,index_GT]

#Estimate a CECM with 1 lagged difference and penalty chosen by the minimum BIC
my_specs <- specs_opt(y,x,p=1,rule="BIC")
coefs <- my_specs$gamma_opt
}
