#' @title Blinded Sample Size Recalculation for a One Subgroup Design
#' @description Given data from an Internal Pilot Study (IPS), \code{bssr.1subgroup} reestimates the nuisance parameters, i.e. variances and prevalence, and recalculates the required sample size for proving a desired alternative when testing for an effect in the full or subpopulation. See 'Details' for more information.
#'
#' @param data   data matrix with data from ongoing trial: see 'Details'.
#' @param alpha  level (type I error) to which the hypothesis is tested.
#' @param beta   type II error (power=1-beta) to which an alternative should be proven.
#' @param delta  vector of treatment effects to be proven, c(outside subgroup, inside subgroup).
#' @param eps    precision parameter concerning the power calculation in the iterative sample size search algorithm.
#' @param approx approximation method: Use a conservative multivariate t distribution ("conservative.t"), a liberal multivariate t distribution ("liberal.t") or a multivariate normal distribution ("normal") to approximate the joint distribution of the standardized test statistics.
#' @param df     in case of a multivariate t distribution approximation, recalculate sample size with degrees of freedom depending on the size of the IPS (df=n1) or depending on the final sample size (df=n).
#' @param adjust adjust blinded estimators for assumed treatment effect ("YES","No").
#' @param k      sample size allocation factor between groups: see 'Details'.
#' @param nmax   maximum total sample size.
#'
#' @details
#' This function performs blinded nuisance parameter reestimation in a design with a subgroup within a full population where we want to test for treatment effects between a control and a treatment group.
#' Then the required sample size for the control and treatment group to prove an existing
#' alternative \code{delta} with a specified power 1-\code{beta} when testing the global null hypothesis \eqn{H_0: \Delta_F=\Delta_S=0} to level \code{alpha} is calculated.
#'
#' The data matrix \code{data} should have three columns: The first column has to be a binary variable (0=treatment group, 1=control group). The second column should also contain a binary variable giving the full population/subgroup differentiation (0=full population, 1=subpopulation). The last column contains the observations.
#'
#' For sample sizes \eqn{n_C} and \eqn{n_T} of the control and treatment group, respectively, the argument \code{k} is the
#' sample size allocation factor, i.e. \eqn{k = n_T/n_C}.
#'
#' The parameter \code{df} provides a difference to the standard sample size calculation procedure implemented in \code{\link{n.1subgroup}}.
#' When applying a multivariate t distribution approximation to approximate the joint distribution of the standardized test statistics it gives the opportunity to use degrees of freedom depending on the number of subjects in the IPS instead of degrees of freedom depending on the projected final sample size.
#' Note that this leads to better performance when dealing with extremely small subgroup sample sizes but significantly increases the calculated final sample size.
#'
#' @return \code{bssr.1subgroup} returns a list containing the recalculated required sample size within the control group and treatment group along with all relevant parameters. Use \code{\link{summary.bssrest}} for a structured overview.
#'
#' @source \code{bssr.1subgroup} uses code contributed by Marius Placzek.
#'
#' @seealso \code{\link{n.1subgroup}} for sample size calculation prior to the trial.
#'
#' @examples
#' #Given data from the Internal Pilot Study, reestimate the nuisance parameters and
#' #recalculate the required sample size to correctly reject with
#' #80% probability when testing the global Nullhypothesis H_0: Delta_F=Delta_S = 0
#' #assuming the true effect Delta_S=1 is in the subgroup (no effect outside of the subgroup).
#'
#' random<-r.1subgroup(n=50, delta=c(0,1), sigma=c(1,1.2), tau=0.4, fix.tau="YES", k=2)
#' reestimate<-bssr.1subgroup(data=random,alpha=0.05,beta=0.1,delta=c(0,1),eps=0.001,
#' approx="conservative.t",df="n1",k=2,adjust="NO")
#' summary(reestimate)
#'
#' @import mvtnorm
#' @export

bssr.1subgroup<-function(data,alpha,beta,delta,eps=0.001, approx=c("conservative.t","liberal.t","normal"),df=c("n","n1"),adjust=c("YES","NO"),k=1,nmax=1000){
  a<-k
  astar<-1+1/a
  tauhat<-sum(data[,2])/length(data[,2])
  deltaFiP<-(1-tauhat)*delta[1]+tauhat*delta[2]
  deltaS<-delta[2]

  #--------------------------------------------------------------------#
  #----------reestimate SDs--------------------------------------------#
  if(adjust=="YES"){
    min.mixed<-function(x, daten, alloc, del){
      sigma<-x
      -sum(log(alloc*dnorm(daten, del, sigma)+(1-alloc)*dnorm(daten, 0, sigma)))
    }
    SigmaShat<-optimize(lower=0,upper=10, f=min.mixed, daten=data[(data$FS==1),3], alloc=1/astar, del=delta[2],tol = .Machine$double.eps*100)$minimum
    SigmaFhat<-sqrt((1-tauhat)*optimize(lower=0,upper=10, f=min.mixed, daten=data[which(data$FS==0),3], alloc=1/astar, del=0,tol = .Machine$double.eps*100)$minimum^2+tauhat*SigmaShat^2)
  }
  if(adjust=="NO"){
    SigmaFhat<-sqrt((1-tauhat)*var(data[which(data$FS==0),3])+tauhat*var(data[which(data$FS==1),3])) #sqrt(1/(2*n1-2)*((2*nfull1-1)*var(data$value[which(data$SF==1)])+(2*nsub1-1)*var(data$value[which(data$SF==0)])))  #sqrt(var(data$value))#
    SigmaShat<-sqrt(var(data[which(data$FS==1),3]))
  }
  sigma<-c(SigmaFhat,SigmaShat)#c(sqrt(tauhat*1^2+(1-tauhat)*1^2), 1)#
  names(sigma)<-c("Full","Sub")
  names(delta)<-c("Full|Sub", "Sub")
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#

  if(df=="n" || approx=="normal"){
    SigmaFffind<-sqrt((SigmaFhat^2-tauhat*SigmaShat^2)/(1-tauhat))
    rec<-n.1subgroup(alpha=alpha,beta=beta, delta=delta, sigma=c(SigmaFffind,SigmaShat), tau=tauhat, eps=eps, approx=approx, k=k,nmax=nmax)
    nrec<-rec$n
  }
  if(df=="n1"){
    #--------------------------------------------------------------------#
    #--------------------------------------------------------------------#
    find<-function(n,alpha,beta,eps,alloc, sigma, delta, tau, approx=c("conservative.t","liberal.t"),n1,z=0){
      a<-alloc
      astar<-1+1/a
      deltaFiP<-(1-tau)*delta[1]+tau*delta[2]
      n_P<-floor((n[2]+n[1])/2)+1
      S<-cov2cor(matrix(c(1+1/4*deltaFiP^2/sigma[1]^2,sqrt(tau)*sigma[2]/sigma[1]*(1+deltaFiP*delta[2]/sigma[1]^2/8),sqrt(tau)*sigma[2]/sigma[1]*(1+deltaFiP*delta[2]/sigma[1]^2/8),1+1/4*delta[2]^2/sigma[2]^2),byrow=TRUE,ncol=2))
      Z<-c(deltaFiP/sigma[1],delta[2]/sigma[2])
      if(z==0){
        if(approx=="conservative.t"){
          z<-qmvt(1-alpha, delta=c(0,0), df=n1[1]-2, corr=S, tail="lower")$quantile
          pow<-1-pmvt(sqrt(c(n_P/astar,tau*n_P/astar))*Z, df=n1[1]-2, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
        }
        if(approx=="liberal.t"){
          z<-qmvt(1-alpha, delta=c(0,0), df=n1[2]-4, corr=S, tail="lower")$quantile
          pow<-1-pmvt(sqrt(c(n_P/astar,tau*n_P/astar))*Z, df=n1[2]-4, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
        }
      }
      if(z!=0){
        if(approx=="conservative.t"){
          pow<-1-pmvt(sqrt(c(n_P/astar,tau*n_P/astar))*Z, df=n1[1]-2, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
        }
        if(approx=="liberal.t"){
          pow<-1-pmvt(sqrt(c(n_P/astar,tau*n_P/astar))*Z, df=n1[2]-4, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
        }
      }

      #print(n_P);print(pow)
      if(n_P!=n[2]){
        if(pow<(1-beta-eps)){
          n_P<-find(n=c(n_P,n[2]),alpha,beta,eps,alloc, sigma, delta, tau, approx,n1,z=z)
        }
        if(pow>(1-beta+eps)){
          n_P<-find(n=c(n[1],n_P),alpha,beta,eps,alloc, sigma, delta, tau, approx,n1,z=z)
        }
      }
      return(n_P)
    }

    #--------------------------------------------------------------------#
    #--------------------------------------------------------------------#
    n1<-c(sum(data[,2]),length(data[,2]))
    n_P<-find(n=c(0,nmax),alpha=alpha,beta=beta,eps=eps,alloc=k, sigma=sigma, delta=delta, tau=tauhat, approx=approx,n1=n1)
    nrec<-c(n_P,a*n_P)
  }

  if (sum(nrec)>nmax){
    Z<-c(deltaFiP/sigma[1],deltaS/sigma[2])
    S<-cov2cor(matrix(c(1+1/4*deltaFiP^2/sigma[1]^2,sqrt(tauhat)*sigma[2]/sigma[1]*(1+deltaFiP*deltaS/sigma[1]^2/8),sqrt(tauhat)*sigma[2]/sigma[1]*(1+deltaFiP*deltaS/sigma[1]^2/8),1+1/4*deltaS^2/sigma[2]^2),byrow=TRUE,ncol=2))
    n_P<-floor(nmax/(a+1))
    ns_P<-floor(tauhat*n_P)
    if(approx=="conservative.t"){
      z<-qmvt(1-alpha, delta=c(0,0), df=(a+1)*ns_P-2, corr=S, tail="lower")$quantile
      pow<-1-pmvt(sqrt(c(n_P/astar,tauhat*n_P/astar))*Z, df=(a+1)*ns_P-2, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
    }
    if(approx=="liberal.t"){
      z<-qmvt(1-alpha, delta=c(0,0), df=(a+1)*n_P-4, corr=S, tail="lower")$quantile
      pow<-1-pmvt(sqrt(c(n_P/astar,tauhat*n_P/astar))*Z, df=(a+1)*n_P-4, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
    }
    if(approx=="normal"){
      z<-qmvnorm(1-alpha, mean=c(0,0), corr=S, tail="lower")$quantile
      pow<-1-pmvnorm(sqrt(c(n_P/astar,tauhat*n_P/astar))*Z, corr=S,lower=c(-Inf,-Inf),upper=c(z,z))[1]
    }
    print("!!Warning!!")
    cat("recalculated total sample size exceeds nmax = ", nmax, ";\n")
    cat("expected power using nmax subjects: ", round(pow,2))
  }

  names(nrec)<-c("Control", "Treatment")
  model<-"normal1subgroup"
  result<-list(n=nrec, alpha=alpha, beta=beta, delta=delta, sigma.est=sigma, tau.est=tauhat, eps=eps, approx=approx, k=a, model=model)
  class(result)<-"bssrest"
  return(result)
}


