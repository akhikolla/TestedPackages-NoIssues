#' Duration data
#'
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
#' a programme is the indicator variable \code{alpha}. While on a programme, the only transition you can
#' make is \code{"job"}. 
#'
#' The dataset is organized as a series of rows for each individual. Each row is a time period
#' with constant covariates.
#'
#' The length of the time period is in the covariate \code{duration}.
#'
#' The transition being made at the end of the period is coded in the covariate \code{d}. This
#' is an integer which is 0 if no transition occurs (e.g. if a covariate changes), it is 1 for
#' the first transition, 2 for the second transition. It can also be a factor, in which case the
#' level marking no transition must be called \code{"none"}. In the test dataset it is a factor
#' with the levels \code{"none"}, \code{"job"}, and \code{"program"}.
#'
#' The covariate \code{alpha} is zero when unemployed, and 1 if on a programme. It is used
#' for two purposes. It is used as an explanatory variable for transition to job, this yields
#' a coefficient which can be interpreted as the effect of being on the programme. It is also
#' used as a "state variable", as an index into a "risk set". I.e. when estimating, the
#' \code{\link{mphcrm}} function must be told which risks/hazards are present. 
#' When on a programme the \code{"toprogram"} transition can not be made. This is implemented
#' by specifying a list of risksets and using \code{alpha+1} as an index into this set.
#'
#' The dataset has already been fitted in the \code{fit} object.
#'
#' @docType data
#' @aliases fit
#' @usage data(durdata)
#'
#' @format A data.frame
#'
#' @keywords datasets
#'
#' @examples
#' data(durdata)
#' print(durdata)
#' print(fit)
#' summary(fit[[1]])
"durdata"
