#' Generate example data
#'
#' @description
#' Generate a data table with example data
#' @details
#' The dataset simulates a labour market programme. People entering the dataset are without a job.
#'
#' They experience two hazards, i.e. probabilities per time period. They can either get a job and exit from
#' the dataset, or they can enter a labour market programme, e.g. a subsidised job or similar, and remain
#' in the dataset and possibly get a job later.
#' In the terms of this package, there are two transitions, \code{"job"} and \code{"program"}.
#' 
#' The two hazards are influenced by covariates observed by the researcher, called \code{"x1"} and
#' \code{"x2"}. In addition there are unobserved characteristics influencing the hazards. Being
#' on a programme also influences the hazard to get a job. In the generated dataset, being on
#' a programme is the indicator variable \code{alpha}. While on a programme, the only transition that can
#' be made is \code{"job"}. 
#'
#' The dataset is organized as a series of rows for each individual. Each row is a time period
#' with constant covariates.
#'
#' The length of the time period is in the covariate \code{duration}.
#'
#' The transition being made at the end of the period is coded in the covariate \code{d}. This
#' is an integer which is 0 if no transition occurs (e.g. if a covariate changes), it is 1 for
#' the first transition, 2 for the second transition. It can also be a factor, in which case the
#' level marking no transition must be called \code{"none"}.
#'
#' The covariate \code{alpha} is zero when unemployed, and 1 if on a programme. It is used
#' for two purposes. It is used as an explanatory variable for transition to job, this yields
#' a coefficient which can be interpreted as the effect of being on the programme. It is also
#' used as a "state variable", as an index into a "risk set". I.e. when estimating, the
#' \code{\link{mphcrm}} function must be told which risks/hazards are present. 
#' When on a programme the \code{"toprogram"} transition can not be made. This is implemented
#' by specifying a list of risksets and using \code{alpha+1} as an index into this set.
#'
#' The two hazards are modeled as \eqn{exp(X \beta + \mu)}, where \eqn{X} is a matrix of covariates
#' \eqn{\beta} is a vector of coefficients to be estimated, and \eqn{\mu} is an intercept. All of
#' these quantities are transition specific. This yields an individual likelihood which we call
#' \eqn{M_i(\mu)}. The idea behind the mixed proportional hazard model is to model the
#' individual heterogeneity as a probability distribution of intercepts. We obtain the individual
#' likelihood \eqn{L_i = \sum_j p_j M_i(\mu_j)}, and, thus, the likelihood \eqn{L = \sum_j L_j}.
#'
#' The likelihood is to be maximized over the parameter vectors \eqn{\beta} (one for each transition),
#' the masspoints \eqn{\mu_j}, and probabilites \eqn{p_j}.
#'
#' The probability distribution is built up in steps. We start with a single masspoint, with
#' probability 1. Then we search for another point with a small probability, and maximize the
#' likelihood from there. We continue with adding masspoints until we no longer can improve
#' the likelihood.
#' 
#' @param N integer.
#' The number of individuals in the dataset.
#' @param censor numeric. The total observation period. Individuals are removed
#' from the dataset if they do not exit to \code{"job"} before this time.
#' @note
#' The example illustrates how \code{data(durdata)} was generated.
#' @examples
#' data.table::setDTthreads(1)  # avoid screams from cran-testing
#' dataset <- datagen(5000,80)
#' print(dataset)
#' risksets <- list(unemp=c("job","program"), onprogram="job")
#' # just two iterations to save time
#' Fit <- mphcrm(d ~ x1+x2 + ID(id) + D(duration) + S(alpha+1) + C(job,alpha),
#'           data=dataset, risksets=risksets,
#'           control=mphcrm.control(threads=1,iters=2))
#' best <- Fit[[1]]
#' print(best)
#' summary(best)
#' @export
datagen <- function(N,censor=80) {
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  id <- seq_len(N)
# two transitions, to exit(1) and to program(2)
# exit is absorbing, program is not, but only exit is allowed afterwards
# if on program one can only exit
  means <- c(-2.5,-3)
  cov2 <- matrix(c(1,0.5,0.5,1),2)
  persons <- data.table(id,x1,x2)
# draw correlated unobserved characteristics
  `:=` <- .N <- ve <- vp <- NULL #avoid check NOTEs
  persons[, c('ve','vp') := {
    vv <- mvtnorm::rmvnorm(.N, mean=means, sigma=cov2)
    list(vv[,1],vv[,2])
  }]
  means <- persons[,colMeans(cbind(job=exp(ve),program=exp(vp)))]
  cv <- persons[,cov(cbind(job=exp(ve),program=exp(vp)))]
  # create spells
  spells <- persons[,{
    genspell(x1,x2,ve,vp,censor)
  }, by=id]
  spells$d <- factor(spells$d,levels=0:2,labels=c('none','job','program'))
  spells$state <- factor(spells$alpha+1,levels=1:2,labels=c('unemp','onprogram'))
  setattr(spells,'means',means)
  setattr(spells,'cov',cv)
  spells
}
