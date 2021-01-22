#' @title Sample Size Calculation for a One Subgroup Design
#' @description \code{n.1subgroup} calculates the required sample size for proving a desired alternative when testing for
#' an effect in the full or subpopulation. See 'Details' for more information.
#'
#' @param alpha   level (type I error) to which the hypothesis is tested.
#' @param beta    type II error (power=1-beta) to which an alternative should be proven.
#' @param delta   vector of treatment effects to be proven, c(outside subgroup, inside subgroup).
#' @param sigma   vector of standard deviations, c(outside subgroup, inside subgroup).
#' @param tau     subgroup prevalence.
#' @param eps     precision parameter concerning the power calculation in the iterative sample size search algorithm.
#' @param approx  approximation method: Use a conservative multivariate t distribution ("conservative.t"), a liberal multivariate t distribution ("liberal.t") or a multivariate normal distribution ("normal") to approximate the joint distribution of the standardized teststatistics.
#' @param k       sample size allocation factor between groups: see 'Details'.
#' @param nmax    maximum total sample size.
#' @param nmin    minimum total sample size.
#'
#' @details
#' This function performs sample size estimation in a design with a subgroup within a full population where we want to test for treatment effects between a control and a treatment group. Since patients from the subgroup might potentially benefit from the treatment more than patients not included in that subgroup, one might prefer testing hypothesis cercerning the full population and the subpopulation at the same time. Here standardized test statistics are their joined distributions are used to calculate the
#' required sample size for the control and treatment group to prove an existing
#' alternative \code{delta} with a specified power 1-\code{beta} when testing the global null hypothesis \eqn{H_0: \Delta_F=\Delta_S=0} to level \code{alpha}.
#'
#' For sample sizes \eqn{n_C} and \eqn{n_T} of the control and treatment group, respectively, the argument \code{k} is the
#' sample size allocation factor, i.e. \eqn{k = n_T/n_C}.
#'
#'
#' @return \code{n.1subgroup} returns the required sample size within the control group and treatment group.
#'
#' @source \code{n.1subgroup} uses code contributed by Marius Placzek.
#'
#' @seealso #' \code{\link{bssr.1subgroup}} for blinded sample size reestimation within a running trial.
#'
#' @examples
#' #Calculate required sample size to correctly reject with
#' #80% probability when testing the global Nullhypothesis H_0: Delta_F=Delta_S = 0
#' #assuming the true effect Delta_S=1 is in the subgroup (no effect outside of the subgroup)
#' #with subgroup prevalence tau=0.4.
#' #The variances in and outside of the subgroup are unequal, sigma=c(1,1.2).
#'
#' estimate<-n.1subgroup(alpha=0.025,beta=0.1,delta=c(0,1),sigma=c(1,1.2),tau=0.4,eps=0.0001,
#' approx="conservative.t",k=2)
#' summary(estimate)
#' @import mvtnorm
#' @import multcomp
#' @export

n.1subgroup<-function(alpha,beta,delta,sigma,tau,eps=0.001, approx=c("conservative.t","liberal.t","normal"),k=1, nmax=1000, nmin=0){ #delta<-c(deltaF\S,deltaS); sigma<-c(sigmaF\S,sigmaS)

  ##------------------------------------------------------------------------------------##
  ##------recursive search algorithm for determining optimal sample size----------------##
  find<-function(n,alpha,beta,eps,alloc, sigma, delta, tau, approx=c("conservative.t","liberal.t","normal")){
    a<-alloc

    astar<-1+1/a

    n_P<-floor((n[2]+n[1])/2)+1
    ns_P<-floor(tau*n_P)

    S<-cov2cor(matrix(c(1+1/4*delta[1]^2/sigma[1]^2,sqrt(tau)*sigma[2]/sigma[1]*(1+delta[1]*delta[2]/sigma[1]^2/8),sqrt(tau)*sigma[2]/sigma[1]*(1+delta[1]*delta[2]/sigma[1]^2/8),1+1/4*delta[2]^2/sigma[2]^2),byrow=TRUE,ncol=2))

    Z<-c(deltaFiP/sigma[1],delta[2]/sigma[2])

    if(approx=="conservative.t"){

      z<-qmvt(1-alpha, delta=c(0,0), df=(a+1)*ns_P-2, corr=S, tail="lower")$quantile
      pow<-1-pmvt(sqrt(c(n_P/astar,tau*n_P/astar))*Z, df=(a+1)*ns_P-2, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]

    }
    if(approx=="liberal.t"){

      z<-qmvt(1-alpha, delta=c(0,0), df=(a+1)*n_P-4, corr=S, tail="lower")$quantile
      pow<-1-pmvt(sqrt(c(n_P/astar,tau*n_P/astar))*Z, df=(a+1)*n_P-4, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]

    }
    if(approx=="normal"){

      z<-qmvnorm(1-alpha, mean=c(0,0), corr=S, tail="lower")$quantile
      pow<-1-pmvnorm(sqrt(c(n_P/astar,tau*n_P/astar))*Z, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]

    }
    if(n_P!=n[2] && ((a+1)*tau*(floor((n[1]+n_P)/2)+1)-2)>1 ){

      if(pow<(1-beta-eps)){
        n_P<-find(n=c(n_P,n[2]),alpha,beta,eps,alloc, sigma, delta, tau, approx)
      }
      if(pow>(1-beta+eps)){
        n_P<-find(n=c(n[1],n_P),alpha,beta,eps,alloc, sigma, delta, tau, approx)
      }
    }
    return(n_P)
  }
  ##------------------------------------------------------------------------------------##
  ##------------------------------------------------------------------------------------##

  approx<-match.arg(approx)
  deltaFiP<-(1-tau)*delta[1]+tau*delta[2]
  sigmaFiP<-sqrt((1-tau)*(sigma[1]^2)+tau*(sigma[2]^2))
  deltaS<-delta[2];sigmaS<-sigma[2]

  n_P<-find(n=c(0,nmax),alpha=alpha,beta=beta,eps=eps,alloc=k,sigma=c(sigmaFiP,sigmaS),delta=c(deltaFiP,deltaS),tau=tau,approx=approx)

  a<-k
  astar<-1+1/a
  n<-ceiling(c(n_P,a*n_P))

  if (sum(n)>nmax){
    Z<-c(deltaFiP/sigmaFiP,deltaS/sigmaS)
    S<-cov2cor(matrix(c(1+1/4*deltaFiP^2/sigmaFiP^2,sqrt(tau)*sigmaS/sigmaFiP*(1+deltaFiP*deltaS/sigmaFiP^2/8),sqrt(tau)*sigmaS/sigmaFiP*(1+deltaFiP*deltaS/sigmaFiP^2/8),1+1/4*deltaS^2/sigmaS^2),byrow=TRUE,ncol=2))
    n_P<-floor(nmax/(a+1))
    ns_P<-floor(tau*n_P)
    if(approx=="conservative.t"){
      z<-qmvt(1-alpha, delta=c(0,0), df=(a+1)*ns_P-2, corr=S, tail="lower")$quantile
      pow<-1-pmvt(sqrt(c(n_P/astar,tau*n_P/astar))*Z, df=(a+1)*ns_P-2, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
    }
    if(approx=="liberal.t"){
      z<-qmvt(1-alpha, delta=c(0,0), df=(a+1)*n_P-4, corr=S, tail="lower")$quantile
      pow<-1-pmvt(sqrt(c(n_P/astar,tau*n_P/astar))*Z, df=(a+1)*n_P-4, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
    }
    if(approx=="normal"){
      z<-qmvnorm(1-alpha, mean=c(0,0), corr=S, tail="lower")$quantile
      pow<-1-pmvnorm(sqrt(c(n_P/astar,tau*n_P/astar))*Z, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
    }
    print("!!Warning!!")
    cat("calculated total sample size exceeds nmax = ", nmax, ";\n")
    cat("expected power using nmax subjects: ", round(pow,2))
  }
  if (sum(n)<nmin){
    Z<-c(deltaFiP/sigmaFiP,deltaS/sigmaS)
    S<-cov2cor(matrix(c(1+1/4*deltaFiP^2/sigmaFiP^2,sqrt(tau)*sigmaS/sigmaFiP*(1+deltaFiP*deltaS/sigmaFiP^2/8),sqrt(tau)*sigmaS/sigmaFiP*(1+deltaFiP*deltaS/sigmaFiP^2/8),1+1/4*deltaS^2/sigmaS^2),byrow=TRUE,ncol=2))
    n_P<-floor(nmin/(a+1))
    ns_P<-floor(tau*n_P)
    if(approx=="conservative.t"){
      z<-qmvt(1-alpha, delta=c(0,0), df=(a+1)*ns_P-2, corr=S, tail="lower")$quantile
      pow<-1-pmvt(sqrt(c(n_P/astar,tau*n_P/astar))*Z, df=(a+1)*ns_P-2, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
    }
    if(approx=="liberal.t"){
      z<-qmvt(1-alpha, delta=c(0,0), df=(a+1)*n_P-4, corr=S, tail="lower")$quantile
      pow<-1-pmvt(sqrt(c(n_P/astar,tau*n_P/astar))*Z, df=(a+1)*n_P-4, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
    }
    if(approx=="normal"){
      z<-qmvnorm(1-alpha, mean=c(0,0), corr=S, tail="lower")$quantile
      pow<-1-pmvnorm(sqrt(c(n_P/astar,tau*n_P/astar))*Z, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
    }
    print("!!Warning!!")
    cat("calculated total sample size below nmin = ", nmin, ";\n")
    cat("expected power using nmin subjects: ", round(pow,2))
  }
  names(n)<-c("Control", "Treatment")
  result<-list(n=n, alpha=alpha, beta=beta, delta=delta, sigma=sigma, tau=tau, eps=eps, approx=approx, k=a, model="normal1subgroup")
  class(result)<-"ssest"
  return(result)
}

