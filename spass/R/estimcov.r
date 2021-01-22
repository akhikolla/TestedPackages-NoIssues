#' @title Estimation of simulation parameters
#' @description \code{estimcov} estimates the covariance matrix and dropout rates given a dataset and observation-times
#'
#' @param data   matrix with the dataset which is used to estimate the covariance and dropout structure.
#' @param Time    vector with observation-times.
#' @param Startvalues   vector with starting values for variance, \code{rho} and \code{theta} respectively.
#' @param stepwidth   vector describing the step length of previously mentioned values.
#' @param maxiter     maximum amount of iterations
#' @param lower       vector with minimum for the parameters described in Startvalues
#' @param upper    vector with maximum for the parameters described in Startvalues
#'
#' @details
#' This function is designed to estimate the variance, \code{rho} and \code{theta} and a vector with the dropout rate in the data.
#'
#' @return \code{estimcov} returns a list with two entries. In the first the parameters variance, \code{rho} and \code{theta} are returned and in the second a vector with the dropout-rate is returned.
#'
#' @source \code{estimcov} uses code contributed by Roland Gerard Gera.
#'
#'
#' @examples
#' # First generate a dataset with 200 patients, rho =0.25 and tau = 0.5 and
#' # then estimate the parameters using estimcov.
#'
#' set.seed(2015)
#' dataset <- r.gee.1subgroup(n=200, reg=list(c(0,0,0,0.1),c(0,0,0,0.1)), sigma=c(3,2.5),
#'   tau=0.5, rho=0.25, theta=1, k=1.5, Time=c(0:5), OD=0)
#'
#' estimations <- estimcov(data=dataset,Time=c(0:5))
#' estimations[[1]]
#' estimations[[2]]
#' 
#' @export


estimcov <- function(data,
                     Time,
                     Startvalues = c(3,0.5,1),
                     stepwidth=c(1e-3,1e-3,1e-3),
                     maxiter=10000,
                     lower = c(0.0001,0.0001,0.0001),
                     upper = c(Inf,5,3)){

  # starting values for optimasation
  estimation = Startvalues

  # get empiric covariance matrix from data
  CovMat=cov(data$y,use="pairwise.complete.obs")

  estimation = optim(par=estimation,
                     fn = optim.function,
                     CovMat=CovMat,
                     Time=Time,
                     method = "L-BFGS-B",
                     lower = lower,
                     upper=  upper,
                     control = list(maxit=maxiter,ndeps=stepwidth))$par

  if (estimation[2]<0) estimation[2]=0
  if (estimation[3]<0) estimation[3]=0.0001

  #### estimation of dropout

  # How many supervised times exist
  m = dim(data$y)[2]

  # empty vector for dropout
  drop=c()

  # Estimate missing subjects
  for (i in 1:m){
    drop[i] = sum(is.na(data$y[,i]))/length(data$y[,i])
  }

  #integrate estimated parameters and estimated dropout in the list "answer"

  answer=list(estimation=estimation,drop=drop)
  return(answer)
}

# function which is to be optimated
optim.function <- function(para,CovMat,Time){

  if (para[1]<=0) para[1]=0.0001
  if (para[2]<0) para[2]=0.0001
  fit = sum(abs(CovMat-gen_cov_cor(var=para[1],rho=para[2],theta=para[3],Time=Time)))
  return(fit)
}
