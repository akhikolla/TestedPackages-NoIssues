#' @title Generate dataset of normal distributed observations in a one subgroup design
#' @description \code{r.1subgroup} generates data for a design with one subgroup within a full population. Each observation is normal distributed with mean 0 in the placebo group and a potential effect in the treatment group. Whether the effect is solely in the subgroup or additionally a certain amount outside of the subgroup can be specified as well as potentially different variances within the subgroup and outside of the subgroup.
#'
#' @param n        number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param delta    vector of treatment effects in the treatment group, c(outside subgroup, within subgroup).
#' @param sigma    vector of standard deviations, c(outside subgroup, inside subgroup).
#' @param tau      subgroup prevalence.
#' @param fix.tau  subgroup prevalence fix or simulated according to tau, see 'Details'.
#' @param k        sample size allocation factor between groups: see 'Details'.
#'
#' @details
#' For \code{delta}\eqn{=(\Delta_F\S, \Delta_S)'} and \code{sigma}\eqn{=(\sigma_F\S, \sigma_S)'}
#' this function \code{r.1subgroup} generates data as follows:
#'
#' Placebo group outside of subgroup \eqn{~N(0,\sigma^2_F\S)},
#' Placebo group within subgroup \eqn{~N(0,\sigma^2_S)},
#' Treatment group outside of subgroup \eqn{~N(\Delta_F\S,\sigma^2_F\S)},
#' Treatment group within subgroup \eqn{~N(\Delta_S,\sigma^2_S)}.
#'
#' If \code{fix.tau=YES} the subgroup size is generated according to the prevalence \code{tau}, i.e. \eqn{n_S=\tau*n}.
#' If \code{fix.tau=YES}, then each new generated observations probability to belong to the subgroup is \eqn{Ber(\code{tau})}
#' distributed and therefore only \eqn{E(n_s)=\tau*n} holds.
#'
#' The argument \code{k} is the
#' sample size allocation factor, i.e. let \eqn{n_C} and \eqn{n_T} denote the sample sizes of of the control and
#' treatment group, respectively, then \eqn{k = n_T/n_C}.
#'
#' @return \code{r.1subgroup} returns a data matrix of dimension \code{n} x \code{3}. The first column \code{TrPl} defines whether
#' the observation belongs to the treatment group (\code{TrPl=0}) or to the placebo group (\code{TrPl=1}). Second column
#' contains the grouping variable \code{FS}. For \code{FS=1} the observation stems from the subgroup, for \code{FS=0} from
#' the full population without the subgroup. In the last column \code{value} the observation can be found.
#' between time points.
#'
#' @source \code{r.1subgroup} uses code contributed by Marius Placzek.
#' 
#' @examples
#'
#' set.seed(142)
#' random<-r.1subgroup(n=50, delta=c(0,1), sigma=c(1,1), tau=0.4, fix.tau="YES", k=2)
#' random 
#' @export
r.1subgroup<-function(n, delta, sigma, tau, fix.tau=c("YES","NO"), k){

  if(length(n)>1){
    n<-length(n)
  }

  ks<-1/(1+k)
  nP<-round(ks*n)
  nT<-n-nP
  if(fix.tau=="NO"){
    T_FS<-(runif(nT,0,1)<tau)*1
    if(((nT-sum(T_FS))<3)){
      T_FS<-c(rep(1,3),rep(0,nT-3))
    }
    if((sum(T_FS)<3)){
      T_FS<-c(rep(1,nT-3),rep(0,3))
    }
    P_FS<-(runif(nP,0,1)<tau)*1
    if(((nP-sum(P_FS))<3)){
      P_FS<-c(rep(1,3),rep(0,nP-3))
    }
    if((sum(P_FS)<3)){
      P_FS<-c(rep(1,nP-3),rep(0,3))
    }
  }
  if(fix.tau=="YES"){
    nP_S<-round(nP*tau)
    nT_S<-round(n*tau)-nP_S
    T_FS<-c(rep(1,nT_S),rep(0,nT-nT_S))        #sub=1,full=0
    P_FS<-c(rep(1,nP_S),rep(0,nP-nP_S))
  }

  T_FS_means<-(1-T_FS)*delta[1]+T_FS*delta[2]
  P_FS_means<-rep(0,nP)
  TP_vars<-(1-c(T_FS,P_FS))*sigma[1]+c(T_FS,P_FS)*sigma[2]
  data<-data.frame(TP=c(rep(0,nT),rep(1,nP)),FS=c(T_FS,P_FS),value=rnorm(n,mean=c(T_FS_means,P_FS_means),sd=TP_vars))
  return(data)
}

