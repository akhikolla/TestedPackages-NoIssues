#' @title Sample Size estimation for longitudinal GEE models
#' @description \code{n.gee.1subgroup} calculates the required sample size for proving a desired alternative when testing a regression coefficients in a full and/or a subpopulation. See 'Details' for more information.
#'
#' @param alpha   level (type I error) to which the hypothesis is tested.
#' @param beta    type II error (power=1-beta) to which an alternative should be proven.
#' @param delta   vector of estimated treatment effect in overall and sub population, c(overall population, only subpopulation).
#' @param sigma   vector of estimated standard deviations, c(full population, subpopulation). See 'Details'.
#' @param tau     subgroup prevalence.
#' @param k       sample size allocation factor between control and treatment: see 'Details'.
#' @param nmax    maximum total sample size.
#' @param npow    calculates power of a test if \code{npow} is a sample size.
#' @param tail    which type of test is used, e.g. which quartile und H0 is calculated.
#'
#' @details
#' This function performs a sample size estimation in a design with a nested subgroup within an overall population. To calculate the required sample only the value of tested regressor needs to inserted as \code{delta}. \code{sigma} is the variance of that regressor.
#' The power for the global null hypothesis is given by 1-\code{beta} and \code{alpha} specifies the false positve level for rejecting \eqn{H_0: \Delta_F=\Delta_S=0} to level \code{alpha}.
#'
#' Here argument \code{k} denotes the
#' sample size allocation factor between treatment groups, i.e. \eqn{k = n_T/n_C}.
#'
#' @return \code{n.gee.1subgroup} returns the required sample size within the control group and treatment group.
#'
#' @source \code{n.gee.1subgroup} uses code contributed by Roland Gerard Gera.
#'
#' @seealso \code{\link{bssr.1subgroup}} for blinded sample size re-estimation within a running trial and \code{\link{sandwich}} for estimating asymptotic covarianc mtrices in GEE models.
#'
#' @examples
#' #Calculate required sample size to correctly reject Null with
#' #80% probability when testing global Nullhypothesis H_0: Delta_F=Delta_S = 0, while
#' #assuming the coefficient in and outside of the subgroup is Delta=c(0.1,0,1) with a
#' #subgroup-prevalence of tau=0.4.
#' #The variances of regressors in delta when variances are unequal sigma=c(0.8,0.4).
#'
#' estimate<-n.gee.1subgroup(alpha=0.05,beta=0.2,delta=c(0.1,0.1),sigma=c(0.8,0.4),tau=0.4, k=1)
#' summary(estimate)
#'
#' #Alternatively we can estimate the power our study would have
#' #if we know the effect in and outside our subgroup as
#' #well as the variance of the regressors. Here we
#' #estimate that only 300 Patiens total can be recruited and we are interested
#' #in the power that would give us.
#'
#' n.gee.1subgroup(alpha=0.05,delta=c(0.1,0.1),sigma=c(0.8,0.4),tau=0.4, k=1, npow=300)
#'
#' @import mvtnorm
#' @export

n.gee.1subgroup <- function(alpha,
                            tail="both",
                            beta=NULL,
                            delta,
                            sigma,
                            tau = 0.5 ,
                            k = 1,
                            npow=NULL,
                            nmax=Inf){

  if((is.null(beta) & is.null(npow)) | (!is.null(beta) & !is.null(npow))){
    stop("Either beta OR npow have to be specified")
  }

  Correl = matrix(c(1,sqrt(tau),sqrt(tau),1),nrow=2)

  CritVal = qmvnorm(1-alpha, mean=c(0,0), corr=Correl, tail=tail)$quantile

  n=min(1000,nmax)
  n_s=round(n*tau)

  if (!is.null(npow)) {
    power=1-pmvnorm(mean=c(sqrt(npow)*delta[1]/sigma[1], sqrt(npow*tau)*delta[2]/sigma[2]), corr=Correl,lower=c(-Inf,-Inf),upper=c(CritVal,CritVal))[1]
    names(power)<-c("The power with given sample size is:")
    return(power)
  }

  power=1-pmvnorm(mean=c(sqrt(n)*delta[1]/sigma[1],
                         sqrt(n_s)*delta[2]/sigma[2]),
                  corr=Correl,lower=c(-Inf,-Inf),upper=c(CritVal,CritVal))[1]

  while(power > 0.80){
    n=n-1
    n_s=round(n*tau)
    power=1-pmvnorm(mean=c(sqrt(n)*delta[1]/sigma[1],
                           sqrt(n_s)*delta[2]/sigma[2]),
                    corr=Correl,lower=c(-Inf,-Inf),upper=c(CritVal,CritVal))[1]

  }

  while(power<0.80){
    if (n>nmax){
      print("!!Warning!!")
      cat("calculated total sample size exceeds nmax = ", nmax, ";\n")
      cat("expected power using nmax subjects: ", round(power,2))
      power=0.8
    } else{
      n=n+1
      n_s=round(n*tau)

      power=1-pmvnorm(mean=c(sqrt(n)*delta[1]/sigma[1],
                             sqrt(n_s)*delta[2]/sigma[2]),
                      corr=Correl,lower=c(-Inf,-Inf),upper=c(CritVal,CritVal))[1]
    }
  }
  n_cont=round(n/(k+1))
  n<-c(n_cont,n-n_cont)

  names(n)<-c("Control", "Treatment")
  result<-list(n=n, alpha=alpha, beta=beta, delta=delta, sigma_reg=sigma, tau=tau, k=k, model="normal.gee.1subgroup")
  class(result)<-"ssest"
  return(result)
}
