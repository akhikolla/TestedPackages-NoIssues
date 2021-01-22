\name{predict.ddp}
\alias{predict.ddp}
\docType{package}
\title{
generate predictions for dependent Dirichlet process Weibull model data without competing risks.
}
\description{
generate predictions for dependent Dirichlet process Weibull model data without competing risks.
}
\usage{
\method{predict}{ddp}(object,newdata,alpha=0.05,tpred=NULL,...)
}
\arguments{
\item{object}{Output from dpweib, must be ddpcomp class}
\item{newdata}{The new dataset for predictions}
\item{tpred}{The time points where the predictions are made. If is not given by the user, it will use the time points where the log hazard ratios are calculated in dpweib function.}
\item{alpha}{\eqn{1-\alpha} is the probability for constructing credible intervals. The default \eqn{\alpha} is 0.05.}
\item{...}{Arguments to be passed to method}
}

\value{
\item{tpred}{The time points where the predictions are made.}
\item{alpha}{\eqn{1-\alpha} is the probability for constructing credible intervals. }
\item{Spred}{A matrix, the estimated survival for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{Spredu}{A matrix, the estimated upper pointwise credible interval of the survival functions for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{Spredl}{A matrix, the estimated lower pointwise credible interval of the survival functions for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{dpred}{A matrix, the estimated density for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{dpredu}{A matrix, the estimated upper pointwise credible interval of the density functions for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{dpredl}{A matrix, the estimated lower pointwise credible interval of the density functions for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{hpred}{A matrix, the estimated hazard for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{hpredu}{A matrix, the estimated upper pointwise credible interval of the hazard functions for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
\item{hpredl}{A matrix, the estimated lower pointwise credible interval of the hazard functions for new covariates. Each row corresponds to a covariate configuration. Each column corresponds to a time point.}
}
