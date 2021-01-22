\name{predict.ddpcomp}
\alias{predict.ddpcomp}
\docType{package}
\title{
generate predictions for dependent Dirichlet process Weibull model data with competing risks.
}
\description{
generate predictions for dependent Dirichlet process Weibull model data with competing risks.
}
\usage{
\method{predict}{ddpcomp}(object,newdata,alpha=0.05,tpred=NULL,...)
}
\arguments{
\item{object}{Output from dpweib, must be ddpcomp class}
\item{newdata}{The new dataset for predictions}
\item{alpha}{\eqn{1-\alpha} is the probability for constructing credible intervals. The default \eqn{\alpha} is 0.05.}
\item{tpred}{The time points where the predictions are made. If is not given by the user, it will use the time points where the log hazard ratios are calculated in dpweib function.}
\item{...}{Arguments to be passed to method}
}
\value{
\item{tpred}{The time points where the predictions are made.}
\item{alpha}{\eqn{1-\alpha} is the probability for constructing credible intervals.}
\item{Fpred}{A matrix, the estimated cumulative incidence functions of cause 1 for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{Fpredu}{A matrix, the estimated upper pointwise credible interval of the cumulative incidence functions of cause 1 for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{Fpredl}{A matrix, the estimated lower pointwise credible interval of the cumulative incidence functions of cause 1 for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{dpred}{A matrix, the estimated subdistribution density functions for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{dpredu}{A matrix, the estimated upper pointwise credible interval of the subdistribution density functions for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{dpredl}{A matrix, the estimated lower pointwise credible interval of the subdistribution density functions for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{hpred}{A matrix, the estimated subdistribution hazard functions for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{hpredu}{A matrix, the estimated upper pointwise credible interval of the subdistribution hazard functions for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{hpredl}{A matrix, the estimated lower pointwise credible interval of the subdistribution hazard functions for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
}
