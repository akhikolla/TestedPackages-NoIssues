#include <iostream>
#include <vector>
#include <Rmath.h>
#include <R.h>
#include <algorithm>
#include <math.h>
#include <Rcpp.h>
#include "arms.h"
#include "commonfunc.h"
#include "auxfuns.h"
#include "noreggilks.h"

using namespace Rcpp;
//calculating the log likelihood
double noreg_loglikelihood(const double tl, const double tr,
const int delta, const int pi,
const double loglambda, const double alpha){
// tl left end point.
// tr right end point.
// delta is the exact observation indicator, 1 if exact observation, 0 if not.
// pi is the right censoring indicator, 1 if right censored, 0 if not.
	if ((pi==0)&(delta == 1)){
//exact obs
	return logdWeibloglambda(tl, alpha, loglambda);
	//return logdWeib(tl,alpha,lambda);
	}else if((pi==1)&(delta == 0)){
// right censored obs
	return logsWeibloglambda(tl,alpha,loglambda);
	}else{
// left censored and interval censored
// to avoid the situation of log(0), we add this cushion thing.
	double cushion=sWeibloglambda(tl,alpha,loglambda)-sWeibloglambda(tr,alpha,loglambda);
//	Rcout<<"leftsurv	"<<sWeibloglambda(tl,alpha,loglambda)<<" rightsurv "<<sWeibloglambda(tr,alpha,loglambda)<<std::endl;
		if((cushion>0.0)&&(testreal(cushion))){
		//	Rcout<<"tl	"<<tl<<" tr "<<tr<<" cushion "<<cushion<<std::endl;
			return log(cushion);
		 }else{
		//	Rcout<<"tl	"<<tl<<" tr "<<tr<<" cushion "<<cushion<<std::endl;
			return -pow(100.0,100.0); 
		 }
	}
}
// Using Neal 8 algorithm to assign cluster to a given observation。
int noreg_group_assign(const double tl, const double tr,
	const int delta,const int pi,
	const int c, const double nu,
	IntegerVector nm,
	NumericVector alpha, NumericVector lambda,
	const double lambda00, const double alpha00, const double alpha0,
	const double alphaalpha, const double alphalambda, const int m,
	IntegerVector allbaskets, std::vector<int> & emptybasket){
// tl left end point.
// tr right end point.
// delta is the exact observation indicator, 1 if exact observation, 0 if not.
// pi is the right censoring indicator, 1 if right censored, 0 if not.
// c is the current cluster that this observation belongs to
// nu concentraton parameter
// nm the number of observatoins each cluster has

	nm(c-1)=nm(c-1)-1;
	int addcluster=m;
	if(nm(c-1)==0){
		emptybasket.insert(emptybasket.begin(),c);
		addcluster=m-1;
	}
	NumericVector alphavector(allbaskets[0]+addcluster);
	NumericVector lambdavector(allbaskets[0]+addcluster);
	NumericVector nmvector(allbaskets[0]+addcluster);
	NumericVector lambda0_new = rgamma(addcluster,alpha00, lambda00);
	for(int i=0;i<(allbaskets[0]+addcluster);i++){
		if(i<allbaskets[0]){
			alphavector(i)=alpha(i);
			lambdavector(i)=lambda(i);
			nmvector(i)=double(nm(i));
		}else{
			lambdavector(i) = R::rgamma(alpha0,1.0/ lambda0_new(i-allbaskets[0]));
			double base_new=findbase(lambdavector(i));
		        if(base_new<80.0){
				double temp1=R::pgamma(base_new,alphaalpha,1.0/alphalambda,1,0);
				alphavector(i)=R::qgamma(temp1+runif(1)[0]*(1.0-temp1),alphaalpha,1.0/alphalambda,1,0);
			}else{
				alphavector(i)=80.0;
			}
			nmvector(i)=nu/double(m);
		}
	}
	if(nm(c-1)==0){
		nmvector(c-1)=nu/double(m);
	}
	NumericVector logprob(allbaskets[0]+addcluster);
	for (int j = 0; j<(allbaskets[0]+addcluster); j++){
			logprob[j] = noreg_loglikelihood(tl,tr,delta,pi, log(lambdavector[j]), alphavector[j]);
	}
    logprob=logprob-max(logprob);
	NumericVector prob=abs(exp(logprob)*nmvector);
	  if(sum(prob)>0.0){
		  prob=prob/sum(prob);
	  }else{
		  prob=nmvector;
		  prob=prob/sum(prob);
	  }
	  //generate a vector from a multivariate distribution
	  IntegerVector newassign=oneMultinomCalt(prob);
	int whichmax=which_max(newassign);
	if((whichmax+1)<=allbaskets[0]){
		nm[whichmax]=nm[whichmax]+1;
		if(whichmax+1==emptybasket[0]){
			emptybasket.erase(emptybasket.begin());
		}
		return whichmax+1;
	}else if(((whichmax+1)>allbaskets[0])&&(emptybasket[0]!=0)){
		nm[emptybasket[0]-1]=1;
		alpha[emptybasket[0]-1]=alphavector(whichmax);
		lambda[emptybasket[0]-1]=lambdavector(whichmax);
		int fill=emptybasket[0];
		nm[emptybasket[0]-1]=1;
		emptybasket.erase(emptybasket.begin());
		return fill;
	}else{
		nm[allbaskets[0]]=1;
		alpha[allbaskets[0]]=alphavector(whichmax);
		lambda[allbaskets[0]]=lambdavector(whichmax);
		allbaskets[0]=allbaskets[0]+1;
		return allbaskets[0];
	}
}

void noreg_update(NumericVector tl, NumericVector tr,
	IntegerVector delta, IntegerVector pi,
	IntegerVector c, IntegerVector nm,
	NumericVector alpha, NumericVector lambda, NumericVector lambda0, 
	const double alpha00,const double alpha0,const double lambda00,
	const double alphaalpha, const double alphalambda,
	double* nextngrp){
		
// tl: the vector of left end point.
// tr: the vector of right end point.
// delta: the vector of exact observation indicator, 1 if exact observation, 0 if not.
// pi: the vector of right censoring indicator, 1 if right censored, 0 if not.
// c is the current cluster that this observation belongs to
// nu concentraton parameter
// c the vector of cluster assignment for all observations

	IntegerVector	uniquec=unique(c);
	*nextngrp=double(uniquec.size());
	for(int k=0;k<*nextngrp;k++){
	NumericVector tltemp=tl[(c==uniquec(k))];
	NumericVector trtemp=tr[(c==uniquec(k))];
	IntegerVector deltatemp=delta[(c==uniquec(k))];
	IntegerVector pitemp=pi[(c==uniquec(k))];
	lambda0[uniquec(k)-1]=R::rgamma(alpha0 + alpha00, 1.0/(lambda00 + lambda[uniquec[k] - 1]));
	noreg_sampleloglambda(&lambda[uniquec[k] - 1], alpha[uniquec[k] - 1],
	alpha0, lambda0[uniquec[k] - 1],
	tltemp.begin(), trtemp.begin(), deltatemp.begin(), pitemp.begin(), tltemp.size());
	noreg_samplealpha(&alpha[uniquec[k] - 1], lambda[uniquec[k] - 1],
	alphaalpha, alphalambda,
	tltemp.begin(), trtemp.begin(), deltatemp.begin(), pitemp.begin(), tltemp.size());
	}
}



List noreg(const int burnin, const int iteration,
	NumericVector tl, NumericVector tr,
	IntegerVector delta,
	IntegerVector pi,
	IntegerVector c,
   	IntegerVector nm, 
	NumericVector alpha,
	NumericVector lambda,
	NumericVector lambda0,
	const double alpha00,
	const double alpha0,
	const double lambda00,
	const double alphaalpha,
	const double alphalambda,
	NumericVector nu,
	NumericVector ngrp,
	const double a, const double b,
	const double ymax, NumericVector t,
	int m, int thin){
	//emptybasket is defined as a series of cluster indicators that have no observations,
	// because it is a cluster indicator, all the elements are positive
	// when you need to create a new cluster, instead of making allbaskets longer,
	// you fill in the emptybasket list first.
    //our algorithm tries to fill in the most recent empty ones first.  
	std::vector<int> emptybasket(1);
	emptybasket[0]=0;
	//allbaskets is not the number of current clusters, 
	//but the number of current clusters plus the emptybasket
	IntegerVector allbaskets(1);
	allbaskets[0]=1;
        int nsave=iteration/thin;
	NumericMatrix alpharec(nsave,tl.size());
	NumericMatrix lambdarec(nsave,tl.size());
	NumericMatrix lambda0rec(nsave,tl.size());
	NumericMatrix lambdarec_orig(nsave,tl.size());
	for (int g = 0; g<(burnin + iteration); g++){
		if((g+1)%100==0){
		Rcout<<g+1<<" iterations out of "<<burnin + iteration<<" iterations done"<<std::endl;
		}
		for (int i = 0; i<tl.size(); i++){
			c[i] = noreg_group_assign(tl[i], tr[i], delta[i], pi[i], c[i], nu[g], nm, alpha, lambda,lambda00, alpha00, alpha0, alphaalpha,alphalambda,m,allbaskets,emptybasket);
		}
		 noreg_update(tl, tr, delta, pi, c, nm, alpha, lambda, lambda0,
		 alpha00, alpha0, lambda00, alphaalpha, alphalambda, &ngrp[g+1]);
		nu[g+1] = nugen(nu[g], tl.size(), ngrp[g+1], a, b);
		if((g>=burnin)&&(g%thin==(thin-1))){
			int index=(g-burnin-(thin-1))/thin;
                for(int z=0;z<tl.size();z++){
				alpharec(index,z) = alpha(c[z] - 1);
				lambda0rec(index,z) = lambda0(c[z] - 1);
				lambdarec(index,z) = lambda(c[z] - 1);
				lambdarec_orig(index,z)=exp(log(lambda(c[z] - 1))-alpha(c[z] - 1)*(log(ymax)-log(10.0)));
			}
		}
	}
	NumericMatrix S(nsave,t.size());
	NumericMatrix d(nsave,t.size());
	NumericMatrix h(nsave,t.size());
	 for(int i=0;i<nsave;i++){
		 for(int j=0;j<t.size();j++){
			 double temp=0.0;
			 double temp1=0.0;
			 for(int k=0;k<tl.size();k++){
				 temp+=sWeib(t(j),alpharec(i,k),lambdarec(i,k));
				 temp1+=dWeib(t(j),alpharec(i,k),lambdarec(i,k));
			 }
			 S(i,j)=temp/tl.size();
			 d(i,j)=temp1/tl.size();
			 h(i,j)=d(i,j)/S(i,j);
		 }
	 }
	 return List::create(
	 Named("c")=c,
	 Named("nm")=nm,
	 Named("S")=S,
	 Named("d")=d,
	 Named("h")=h,
	 Named("ngrp")=ngrp,
    	 Named("alpharec")=alpharec,
         Named("lambda0rec")=lambda0rec,
     	 Named("lambdascaled")=lambdarec,
         Named("lambdarec")=lambdarec_orig,
	 Named("emptybasket")=emptybasket,
	 Named("allbaskets")=allbaskets);
}



List noreg_resume(const int burnin, const int iteration,
	NumericVector tl, NumericVector tr,
	IntegerVector delta,
	IntegerVector pi,
	IntegerVector c,
   	IntegerVector nm, 
	NumericVector alpha,
	NumericVector lambda,
	NumericVector lambda0,
	const double alpha00,
	const double alpha0,
	const double lambda00,
	const double alphaalpha,
	const double alphalambda,
	NumericVector nu,
	NumericVector ngrp,
	const double a, const double b,
	const double ymax, NumericVector t,
	int m, int thin, std::vector<int> emptybasket, IntegerVector allbaskets){
        int nsave=iteration/thin;
	NumericMatrix alpharec(nsave,tl.size());
	NumericMatrix lambda0rec(nsave,tl.size());
	NumericMatrix lambdarec(nsave,tl.size());
	NumericMatrix lambdarec_orig(nsave,tl.size());
	for (int g = 0; g<(burnin + iteration); g++){
		if((g+1)%100==0){
		Rcout<<g+1<<" iterations out of "<<burnin + iteration<<" iterations done"<<std::endl;
		}
		for (int i = 0; i<tl.size(); i++){
			c[i] = noreg_group_assign(tl[i], tr[i], delta[i], pi[i], c[i], nu[g], nm, alpha, lambda,lambda00, alpha00, alpha0, 				alphaalpha,alphalambda,m,allbaskets,emptybasket);
		}
		 noreg_update(tl, tr, delta, pi, c, nm, alpha, lambda, lambda0,
		 alpha00, alpha0, lambda00, alphaalpha, alphalambda, &ngrp[g+1]);
		nu[g+1] = nugen(nu[g], tl.size(), ngrp[g+1], a, b);
		if((g>=burnin)&&(g%thin==(thin-1))){
			int index=(g-burnin-(thin-1))/thin;
                for(int z=0;z<tl.size();z++){
				alpharec(index,z) = alpha(c[z] - 1);
				lambda0rec(index,z) = lambda0(c[z] - 1);
				lambdarec(index,z) = lambda(c[z] - 1);
				lambdarec_orig(index,z)=exp(log(lambda(c[z] - 1))-alpha(c[z] - 1)*(log(ymax)-log(10.0)));
			}
		}
	}
	NumericMatrix S(nsave,t.size());
	NumericMatrix d(nsave,t.size());
	NumericMatrix h(nsave,t.size());
	 for(int i=0;i<nsave;i++){
		 for(int j=0;j<t.size();j++){
			 double temp=0.0;
			 double temp1=0.0;
			 for(int k=0;k<tl.size();k++){
				 temp+=sWeib(t(j),alpharec(i,k),lambdarec(i,k));
				 temp1+=dWeib(t(j),alpharec(i,k),lambdarec(i,k));
			 }
			 S(i,j)=temp/tl.size();
			 d(i,j)=temp1/tl.size();
			 h(i,j)=d(i,j)/S(i,j);
		 }
	 }
	 return List::create(
	 Named("c")=c,
	 Named("nm")=nm,
	 Named("S")=S,
	 Named("d")=d,
	 Named("h")=h,
	 Named("ngrp")=ngrp,
    	 Named("alpharec")=alpharec,
         Named("lambda0rec")=lambda0rec,
     	 Named("lambdascaled")=lambdarec,
         Named("lambdarec")=lambdarec_orig,
	 Named("emptybasket")=emptybasket,
	 Named("allbaskets")=allbaskets);
}
