#' @title Simulation of a One Subgroup Design with Internal Pilot Study
#' @description Given estimates of the treatment effects to be proven, the variances, and the prevalence,
#' \code{sim.bssr.1subgroup} calculates a initial sample size and performes a blinded sample size recalculation
#' after a prespecified number of subjects have been enrolled. Each oberservation is simulated and a final analysis executed.
#' Several variations are included, such as different approximations or sample size allocation.
#'
#' @param nsim     number of simulation runs.
#' @param alpha    level (type I error) to which the hypothesis is tested.
#' @param beta     type II error (power=1-beta) to which an alternative should be proven.
#' @param delta    vector of true treatment effects, c(outside subgroup, inside subgroup).
#' @param sigma    vector of true standard deviations, c(outside subgroup, inside subgroup).
#' @param tau      subgroup prevalence.
#' @param vdelta   vector of treatment effects to be proven, c(outside subgroup, inside subgroup).
#' @param vsigma   vector of assumed standard deviations, c(outside subgroup, inside subgroup).
#' @param vtau     expected subgroup prevalence.
#' @param rec.at   blinded sample size review is performed after \code{rec.at}*\eqn{100\%} subjects of the initial sample size calculation.
#' @param eps      precision parameter concerning the power calculation in the iterative sample size search algorithm.
#' @param approx   approximation method: Use a conservative multivariate t distribution ("conservative.t"), a liberal multivariate t distribution ("liberal.t") or a multivariate normal distribution ("normal") to approximate the joint distribution of the standardized teststatistics.
#' @param df       in case of a multivariate t distribution approximation, recalculate sample size with degrees of freedom depending on the size of the IPS (df=n1) or depending on the final sample size (df=n).
#' @param fix.tau  subgroup prevalence is fixed by design (e.g. determined by recruitment) or is simulated and has to be reestimated during the blinded review.
#' @param k        sample size allocation factor between groups: see 'Details'.
#' @param adjust   adjust blinded estimators for assumed treatment effect ("YES","No").
#'
#' @details
#' This function combines sample size estimation, blinded sample size reestimation and analysis in a design with a subgroup within a full population where we want to test for treatment effects between a control and a treatment group.
#' The required sample size for the control and treatment group to prove an existing
#' alternative \code{delta} with a specified power 1-\code{beta} when testing the global null hypothesis \eqn{H_0: \Delta_F=\Delta_S=0} to level \code{alpha} is calculated prior to the study and then recalculated in an internal pilot study.
#'
#' For sample sizes \eqn{n_C} and \eqn{n_T} of the control and treatment group, respectively, the argument \code{k} is the
#' sample size allocation factor, i.e. \eqn{k = n_T/n_C}.
#'
#' The parameter \code{df} provides a difference to the standard sample size calculation procedure implemented in \code{\link{n.1subgroup}}.
#' When applying a multivariate t distribution approximation to approximate the joint distribution of the standardized test statistics it gives the opportunity to use degrees of freedom depending on the number of subjects in the IPS instead of degrees of freedom depending on the projected final sample size.
#' Note that this leads to better performance when dealing with extremely small subgroup sample sizes but significantly increases the calculated final sample size.
#'
#' @return \code{sim.bssr.1subgroup} returns a data.frame containing the mean recalculated sample size within the control group and treatment group and the achieved simulated power along with all relevant parameters.
#'
#' @source \code{sim.bssr.1subgroup} uses code contributed by Marius Placzek.
#'
#' @seealso \code{sim.bssr.1subgroup} makes use of \code{\link{n.1subgroup}}, \code{\link{bssr.1subgroup}}, and \code{\link{r.1subgroup}}.
#'
#' @examples
#' sim.bssr.1subgroup(nsim=10,alpha=0.025,beta=0.1,delta=c(0,1),sigma=c(1,1.3),tau=0.2,
#' vdelta=c(0,1),vsigma=c(1,1),vtau=0.3,eps=0.002, approx="conservative.t",df="n",
#' fix.tau="YES",k=1,adjust="NO")
#'
#' @import mvtnorm
#' @import multcomp
#' @export

sim.bssr.1subgroup<-function(nsim=1000,alpha,beta,delta,sigma,tau,vdelta,vsigma,vtau,rec.at=1/2,eps=0.001, approx=c("conservative.t","liberal.t","normal"),df=c("n","n1"), fix.tau=c("YES","NO"),k=1,adjust=c("YES","NO")){
  ntruelist<-n.1subgroup(alpha=alpha,beta=beta,delta=delta,sigma=sigma,tau=tau,eps=eps, approx=approx,k=k)#
  ntrue<-ntruelist$n
  ninit<-n.1subgroup(alpha=alpha,beta=beta,delta=vdelta,sigma=vsigma,tau=vtau,eps=eps, approx="conservative.t",k=k)$n
  nrecs<-c()
  sigmas<-c()
  sigmas2<-c()
  means<-c()
  reject<-0
  for (zz in 1:nsim){
    startn<-floor(rec.at*sum(ninit))
    data1<-r.1subgroup(n=startn, delta=delta, sigma=sigma, tau=tau, fix.tau=fix.tau, k=k)
    nsub1<-sum(data1[,2])
    nrec <-bssr.1subgroup(data=data1,alpha=alpha,beta=beta,delta=vdelta,eps=eps,approx=approx,df=df,k=k,adjust=adjust,nmax=3*sum(ntrue))
    nrecs<-rbind(nrecs,nrec$n)
    sigmas<-rbind(sigmas, nrec$sigma.est)
    nfinal<-max(sum(ninit),sum(nrec$n))
    if(nfinal>startn){
      data2<-r.1subgroup(n=nfinal-startn, delta=delta, sigma=sigma, tau=tau, fix.tau=fix.tau, k=k)
      data<-rbind(data1,data2)
    }
    else{
      data<-data1
    }

    astar=1+1/k
    nrec<-length(data[,1])
    nP=sum(data$TP)
    nT=sum(1-data$TP)
    nSP=sum((data$FS)*(data$TP))
    nFP=sum((1-data$FS)*data$TP)
    nST=sum(data$FS*(1-data$TP))
    nFT=sum((1-data$FS)*(1-data$TP))
    tauhat2<-(nSP+nST)/(nP+nT)
    SigmaFhat2<-sqrt(1/(nrec-4)*((nFT-1)*var(data$value[(data$FS==0)& (data$TP==0)])+(nFP-1)*var(data$value[(data$FS==0)& (data$TP==1)])+(nST-1)*var(data$value[(data$FS==1)& (data$TP==0)])+(nSP-1)*var(data$value[(data$FS==1)& (data$TP==1)])))
    SigmaShat2<-sqrt(1/((nST+nSP)-2)*((nST-1)*var(data$value[(data$FS==1)& (data$TP==0)])+(nSP-1)*var(data$value[(data$FS==1)& (data$TP==1)])))
    sigmas2<-rbind(sigmas2,c(SigmaFhat2,SigmaShat2))

    MeanFhat<-mean( data$value[(data$TP==0)] )-mean( data$value[(data$TP==1)] )
    MeanShat<-mean( data$value[(data$FS==1)& (data$TP==0)] )-mean( data$value[(data$FS==1)& (data$TP==1)] )
    means<-rbind(means,c(MeanFhat,MeanShat,tauhat2))
    Shat2<-cov2cor(matrix(c(1,sqrt(tauhat2)*SigmaShat2/SigmaFhat2,sqrt(tauhat2)*SigmaShat2/SigmaFhat2,1),byrow=TRUE,ncol=2))
    zhat2=qmvt(1-alpha, delta=c(0,0),df=(nSP+nST-2), corr=Shat2, tail="lower")$quantile

    Z<-c(sqrt(nP/astar)*MeanFhat/SigmaFhat2,sqrt(tauhat2*nP/astar)*MeanShat/SigmaShat2)
    reject<-reject+(sum(Z>zhat2)>0)
    print(zz)
  }
  nrecsum<-nrecs[,1]+nrecs[,2]
  quan5<-sort(nrecsum)[floor(0.05*nsim)+1]
  quan95<-sort(nrecsum)[ceiling(0.95*nsim)]
  quan50<-median(nrecsum)
  stdnrec<-sd(nrecsum)
  result<-data.frame(ntrue=t(ntrue),ninit=t(ninit),rec.at=rec.at,nsub1=nsub1,nrec=t(colMeans(nrecs)),pow=reject/nsim,sig.bssr=t(colMeans(sigmas)),sig.fin=t(colMeans(sigmas2)),meansTau=t(colMeans(means)),nsim=nsim,alpha=alpha,beta=beta,delta=t(delta),sigma=t(sigma),tau=tau,vdelta=t(vdelta),vsigma=t(vsigma),vtau=vtau,eps=eps, approx=approx,df=df, fix.tau=fix.tau,k=k,adj=adjust,sd=stdnrec,quan5,quan50,quan95)
  return(result)
}

