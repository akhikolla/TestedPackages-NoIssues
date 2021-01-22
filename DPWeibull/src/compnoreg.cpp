#include <iostream>
#include <vector>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <algorithm>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <Rcpp.h>
using namespace Rcpp;
#include "arms.h"
#include "auxfuns.h"
#include "commonfunc.h"
#include "compnoreggilks.h"
double compnoreg_loglikelihood(const double t, const int delta,
const double lambda1, const double alpha1,
const double lambda2, const double alpha2,
const double p){
	if (delta == 1){
		return log(p)+logdWeib(t,alpha1,lambda1);
	}else if(delta==2){
		return log(1.0-p)+logdWeib(t,alpha2,lambda2);
	}else{
		if(p*pWeib(t,alpha1,lambda1)+(1.0-p)*pWeib(t,alpha2,lambda2)<1.0){
		return log(1.0-p*pWeib(t,alpha1,lambda1)-(1.0-p)*pWeib(t,alpha2,lambda2));
		}else{
			return -pow(100.0,100.0);
		}
	}
}



int compnoreg_group_assign(const double t,const int delta,
	const int c, const double nu,
	IntegerVector nm,
	NumericVector alpha1, NumericVector lambda1,
	NumericVector alpha2, NumericVector lambda2,
	NumericVector p,
	const double lambda00, const double alpha00, const double alpha0,
	const double alphaalpha, const double alphalambda,
	const double gamma0, const double gamma1,
	int m,
	IntegerVector allbaskets, std::vector<int> & emptybasket){
	nm(c-1)=nm(c-1)-1;
	int addcluster=m;
	if(nm(c-1)==0){
		emptybasket.insert(emptybasket.begin(),c);
		addcluster=m-1;
	}
	NumericVector alphavector1(allbaskets[0]+addcluster);
	NumericVector lambdavector1(allbaskets[0]+addcluster);
	NumericVector alphavector2(allbaskets[0]+addcluster);
	NumericVector lambdavector2(allbaskets[0]+addcluster);
	NumericVector pvector(allbaskets[0]+addcluster);
	NumericVector nmvector(allbaskets[0]+addcluster);
	for(int i=0;i<(allbaskets[0]+addcluster);i++){
		if(i<allbaskets[0]){
			alphavector1(i)=alpha1(i);
			lambdavector1(i)=lambda1(i);
			alphavector2(i)=alpha2(i);
			lambdavector2(i)=lambda2(i);
			pvector(i)=p(i);
			nmvector(i)=double(nm(i));
		}else{
			pvector(i)=R::rbeta(gamma0,gamma1);
			double lambda0_new1=R::rgamma(alpha00,1.0/lambda00);
			lambdavector1(i) = R::rgamma(alpha0,1.0/ lambda0_new1);
			double lambda0_new2=R::rgamma(alpha00,1.0/lambda00);
			lambdavector2(i) = R::rgamma(alpha0,1.0/ lambda0_new2);
			double base_new=findbase(lambdavector1(i));
			if(base_new<80.0){
				double temp1=R::pgamma(base_new,alphaalpha,1.0/alphalambda,1,0);
				alphavector1(i)=R::qgamma(temp1+runif(1)[0]*(1.0-temp1),alphaalpha,1.0/alphalambda,1,0);
			}else{
				alphavector1(i)=80.0;
			}
			base_new=findbase(lambdavector2(i));
			if(base_new<80.0){
				double temp1=R::pgamma(base_new,alphaalpha,1.0/alphalambda,1,0);
				alphavector2(i)=R::qgamma(temp1+runif(1)[0]*(1.0-temp1),alphaalpha,1.0/alphalambda,1,0);
			}else{
				alphavector2(i)=80.0;
			}
			nmvector(i)=nu/double(m);
		}
	}
	if(nm(c-1)==0){
		nmvector(c-1)=nu/double(m);
	}
	
	NumericVector logprob(allbaskets[0]+addcluster);
	for (int j = 0; j<(allbaskets[0]+addcluster); j++){
			logprob[j] =compnoreg_loglikelihood(t,delta,lambdavector1[j], alphavector1[j],
			 lambdavector2[j], alphavector2[j],pvector[j]);
	}
    logprob=logprob-max(logprob);
	NumericVector prob=exp(logprob)*nmvector;
	if(sum(prob)==0.0){
		prob=nmvector;
	}
	prob=prob/sum(prob);
	IntegerVector newassign=oneMultinomCalt(prob);
	int whichmax=which_max(newassign);
	if((whichmax+1)<=allbaskets[0]){
		nm[whichmax]=nm[whichmax]+1;
		if(whichmax+1==emptybasket[0]){
			emptybasket.erase(emptybasket.begin());
		}
		return whichmax+1;
	}else if(((whichmax+1)>allbaskets[0])&&(emptybasket[0]!=0)){
		alpha1[emptybasket[0]-1]=alphavector1(whichmax);
		lambda1[emptybasket[0]-1]=lambdavector1(whichmax);
		alpha2[emptybasket[0]-1]=alphavector2(whichmax);
		lambda2[emptybasket[0]-1]=lambdavector2(whichmax);
		p[emptybasket[0]-1]=pvector(whichmax);
		int fill=emptybasket[0];
		nm[emptybasket[0]-1]=1;
		emptybasket.erase(emptybasket.begin());
		return fill;
	}else{
		nm[allbaskets[0]]=1;
		alpha1[allbaskets[0]]=alphavector1(whichmax);
		lambda1[allbaskets[0]]=lambdavector1(whichmax);
		alpha2[allbaskets[0]]=alphavector2(whichmax);
		lambda2[allbaskets[0]]=lambdavector2(whichmax);
		p[allbaskets[0]]=pvector(whichmax);
		allbaskets[0]=allbaskets[0]+1;
		return allbaskets[0];
	}
}

void compnoreg_update(NumericVector t,
	IntegerVector delta, 
	IntegerVector c, IntegerVector nm,
	NumericVector alpha1, NumericVector lambda1, NumericVector lambda01, 
	NumericVector alpha2, NumericVector lambda2, NumericVector lambda02, 
	NumericVector p,
	const double alpha00, const double alpha0, const double lambda00,
	const double alphaalpha, const double alphalambda,
	const double gamma0, const double gamma1,
	double* nextngrp){
	IntegerVector	uniquec=unique(c);
	*nextngrp=double(uniquec.size());
	for(int k=0;k<*nextngrp;k++){
	NumericVector ttemp=t[c==uniquec(k)];
	IntegerVector deltatemp=delta[c==uniquec(k)];
	lambda01[uniquec(k)-1]=R::rgamma(alpha0 + alpha00, 1.0/(lambda00 + lambda1[uniquec[k] - 1]));
	compnoreg_sampleloglambda(&lambda1[uniquec[k] - 1], alpha0, lambda01[uniquec(k)-1],
	alpha1[uniquec[k] - 1],alpha2[uniquec[k] - 1],lambda2[uniquec[k] - 1],p[uniquec[k] - 1],
	ttemp.begin(),deltatemp.begin(),ttemp.size());
	compnoreg_samplealpha(&alpha1[uniquec[k] - 1],lambda1[uniquec[k] - 1],alpha2[uniquec[k] - 1],
	lambda2[uniquec[k] - 1],p[uniquec[k] - 1],alphaalpha,alphalambda,
	ttemp.begin(),deltatemp.begin(),ttemp.size());
	lambda02[uniquec(k)-1]=R::rgamma(alpha0 + alpha00, 1.0/(lambda00 + lambda2[uniquec[k] - 1]));
	compnoreg_sampleloglambda(&lambda2[uniquec[k] - 1], alpha0, lambda02[uniquec(k)-1],
	alpha2[uniquec[k] - 1],alpha1[uniquec[k] - 1],lambda1[uniquec[k] - 1],(1.0-p[uniquec[k] - 1]),
	ttemp.begin(),deltatemp.begin(),ttemp.size());
	compnoreg_samplealpha(&alpha2[uniquec[k] - 1],lambda2[uniquec[k] - 1],alpha1[uniquec[k] - 1],
	lambda1[uniquec[k] - 1],(1.0-p[uniquec[k] - 1]),alphaalpha,alphalambda,
	ttemp.begin(),deltatemp.begin(),ttemp.size());
	compnoreg_samplep(&p[uniquec[k] - 1],gamma0,gamma1,ttemp.begin(),deltatemp.begin(),
	alpha1[uniquec[k] - 1],lambda1[uniquec[k] - 1],alpha2[uniquec[k] - 1],lambda2[uniquec[k] - 1],
	ttemp.size());
	}
}



List compnoreg(const int burnin, const int iteration,
	NumericVector t,
	IntegerVector delta,
	IntegerVector c,
    IntegerVector nm, 
	NumericVector alpha1,
	NumericVector lambda1,
	NumericVector lambda01,
	NumericVector alpha2,
	NumericVector lambda2,
	NumericVector lambda02,
	NumericVector p,
	const double alpha00,
	const double alpha0,
	const double lambda00,
	const double alphaalpha,
	const double alphalambda,
	const double gamma0,
	const double gamma1,
	NumericVector nu,
	NumericVector ngrp,
	const double a, const double b,
	const double ymax, NumericVector tplot,
	int m, int thin){
	std::vector<int> emptybasket(1);
	emptybasket[0]=0;
	IntegerVector allbaskets(1);
	allbaskets[0]=1;
	int nsave=iteration/thin;
	NumericMatrix alpharec1(nsave,t.size());
	NumericMatrix lambdarec01(nsave,t.size());
	NumericMatrix lambdarec1(nsave,t.size());
	NumericMatrix lambdarec_orig1(nsave,t.size());
	NumericMatrix alpharec2(nsave,t.size());
	NumericMatrix lambdarec02(nsave,t.size());
	NumericMatrix lambdarec2(nsave,t.size());
	NumericMatrix lambdarec_orig2(nsave,t.size());
	NumericMatrix prec(nsave,t.size());	
	 for (int g = 0; g<(burnin + iteration); g++){
		if((g+1)%100==0){
		Rcout<<g+1<<" iterations out of "<<burnin + iteration<<" iterations done"<<std::endl;
		}
		for (int i = 0; i<t.size(); i++){
			c[i] =compnoreg_group_assign(t[i], delta[i],c[i], nu[g], 
			nm, alpha1, lambda1,alpha2, lambda2, p,
			lambda00, alpha00, alpha0, alphaalpha,alphalambda,
			gamma0,gamma1,
			m,allbaskets,emptybasket);
		}
		 compnoreg_update(t, delta, c, nm, alpha1, lambda1, lambda01,
		 alpha2, lambda2, lambda02,p,
		 alpha00, alpha0, lambda00,
		 alphaalpha, alphalambda,
		 gamma0,gamma1,
		 &ngrp[g+1]);		 
		nu[g+1] = nugen(nu[g], t.size(), ngrp[g+1], a, b);
		if((g>=burnin)&&(g%thin==(thin-1))){
		int iter=(g-burnin-(thin-1))/thin;
			for (int z = 0; z<t.size(); z++){
				alpharec1(iter,z) = alpha1(c[z] - 1);
				lambdarec01(iter,z) = lambda01(c[z] - 1);
				lambdarec1(iter,z) = lambda1(c[z] - 1);
				lambdarec_orig1(iter,z)=exp(log(lambda1(c[z] - 1))-alpha1(c[z] - 1)*(log(ymax)-log(10.0)));
				alpharec2(iter,z) = alpha2(c[z] - 1);
				lambdarec02(iter,z) = lambda02(c[z] - 1);
				lambdarec2(iter,z) = lambda2(c[z] - 1);
				lambdarec_orig2(iter,z)=exp(log(lambda2(c[z] - 1))-alpha2(c[z] - 1)*(log(ymax)-log(10.0)));
				prec(iter,z) = p(c[z] - 1);
			}
		}		
	}
	
	NumericMatrix CIF1(nsave,tplot.size());
	NumericMatrix CIF2(nsave,tplot.size());
	NumericMatrix d1(nsave,tplot.size());
	NumericMatrix d2(nsave,tplot.size());
	NumericMatrix h1(nsave,tplot.size());
	NumericMatrix h2(nsave,tplot.size());
	for(int i=0;i<nsave;i++){
		for(int j=0;j<tplot.size();j++){
		double temp1=0.0;
		double temp2=0.0;
                double temp3=0.0;
		double temp4=0.0;
		for(int k=0;k<t.size();k++){
			temp1+=prec(i,k)*pWeib(tplot(j),alpharec1(i,k),lambdarec1(i,k));
			temp2+=(1.0-prec(i,k))*pWeib(tplot(j),alpharec2(i,k),lambdarec2(i,k));
			temp3+=prec(i,k)*dWeib(tplot(j),alpharec1(i,k),lambdarec1(i,k));
			temp4+=(1.0-prec(i,k))*dWeib(tplot(j),alpharec2(i,k),lambdarec2(i,k));
		}
		CIF1(i,j)=temp1/t.size();
		CIF2(i,j)=temp2/t.size();
		d1(i,j)=temp3/t.size();
		d2(i,j)=temp4/t.size();
		h1(i,j)=d1(i,j)/(1.0-CIF1(i,j));
		h2(i,j)=d2(i,j)/(1.0-CIF2(i,j));
		}
	}

	 return List::create(
	 Named("c")=c,
	 Named("nm")=nm,
	 Named("alpharec1")=alpharec1,
	 Named("lambda0rec1")=lambdarec01,
	 Named("lambdarec1")=lambdarec_orig1,
	 Named("lambdascaled1")=lambdarec1,
	 Named("alpharec2")=alpharec2,
	 Named("lambda0rec2")=lambdarec02,
	 Named("lambdarec2")=lambdarec_orig2,
	 Named("lambdascaled2")=lambdarec2,
	 Named("prec")=prec,
         Named("CIF1")=CIF1,
	 Named("CIF2")=CIF2,
         Named("d1")=d1,
	 Named("d2")=d2,
         Named("h1")=h1,
	 Named("h2")=h2,
	 Named("ngrp")=ngrp,
	 Named("emptybasket")=emptybasket,
	 Named("allbaskets")=allbaskets);
}




List compnoreg_resume(const int burnin, const int iteration,
	NumericVector t,
	IntegerVector delta,
	IntegerVector c,
    	IntegerVector nm, 
	NumericVector alpha1,
	NumericVector lambda1,
	NumericVector lambda01,
	NumericVector alpha2,
	NumericVector lambda2,
	NumericVector lambda02,
	NumericVector p,
	const double alpha00,
	const double alpha0,
	const double lambda00,
	const double alphaalpha,
	const double alphalambda,
	const double gamma0,
	const double gamma1,
	NumericVector nu,
	NumericVector ngrp,
	const double a, const double b,
	const double ymax, NumericVector tplot,
	int m, int thin, std::vector<int> emptybasket, IntegerVector allbaskets){
	int nsave=iteration/thin;
	NumericMatrix alpharec1(nsave,t.size());
	NumericMatrix lambdarec01(nsave,t.size());
	NumericMatrix lambdarec1(nsave,t.size());
	NumericMatrix lambdarec_orig1(nsave,t.size());
	NumericMatrix alpharec2(nsave,t.size());
	NumericMatrix lambdarec02(nsave,t.size());
	NumericMatrix lambdarec2(nsave,t.size());
	NumericMatrix lambdarec_orig2(nsave,t.size());
	NumericMatrix prec(nsave,t.size());	
	 for (int g = 0; g<(burnin + iteration); g++){
		if((g+1)%100==0){
		Rcout<<g+1<<" iterations out of "<<burnin + iteration<<" iterations done"<<std::endl;
		}
		for (int i = 0; i<t.size(); i++){
			c[i] =compnoreg_group_assign(t[i], delta[i],c[i], nu[g], 
			nm, alpha1, lambda1,alpha2, lambda2, p,
			lambda00, alpha00, alpha0, alphaalpha,alphalambda,
			gamma0,gamma1,
			m,allbaskets,emptybasket);
		}
		 compnoreg_update(t, delta, c, nm, alpha1, lambda1, lambda01,
		 alpha2, lambda2, lambda02,p,
		 alpha00, alpha0, lambda00,
		 alphaalpha, alphalambda,
		 gamma0,gamma1,
		 &ngrp[g+1]);		 
		nu[g+1] = nugen(nu[g], t.size(), ngrp[g+1], a, b);
		if((g>=burnin)&&(g%thin==(thin-1))){
		int iter=(g-burnin-(thin-1))/thin;
			for (int z = 0; z<t.size(); z++){
				alpharec1(iter,z) = alpha1(c[z] - 1);
				lambdarec01(iter,z) = lambda01(c[z] - 1);
				lambdarec1(iter,z) = lambda1(c[z] - 1);
				lambdarec_orig1(iter,z)=exp(log(lambda1(c[z] - 1))-alpha1(c[z] - 1)*(log(ymax)-log(10.0)));
				alpharec2(iter,z) = alpha2(c[z] - 1);
				lambdarec02(iter,z) = lambda02(c[z] - 1);
				lambdarec2(iter,z) = lambda2(c[z] - 1);
				lambdarec_orig2(iter,z)=exp(log(lambda2(c[z] - 1))-alpha2(c[z] - 1)*(log(ymax)-log(10.0)));
				prec(iter,z) = p(c[z] - 1);
			}
		}		
	}
	
	NumericMatrix CIF1(nsave,tplot.size());
	NumericMatrix CIF2(nsave,tplot.size());
	NumericMatrix d1(nsave,tplot.size());
	NumericMatrix d2(nsave,tplot.size());
	NumericMatrix h1(nsave,tplot.size());
	NumericMatrix h2(nsave,tplot.size());
	for(int i=0;i<nsave;i++){
		for(int j=0;j<tplot.size();j++){
		double temp1=0.0;
		double temp2=0.0;
                double temp3=0.0;
		double temp4=0.0;
		for(int k=0;k<t.size();k++){
			temp1+=prec(i,k)*pWeib(tplot(j),alpharec1(i,k),lambdarec1(i,k));
			temp2+=(1.0-prec(i,k))*pWeib(tplot(j),alpharec2(i,k),lambdarec2(i,k));
			temp3+=prec(i,k)*dWeib(tplot(j),alpharec1(i,k),lambdarec1(i,k));
			temp4+=(1.0-prec(i,k))*dWeib(tplot(j),alpharec2(i,k),lambdarec2(i,k));
		}
		CIF1(i,j)=temp1/t.size();
		CIF2(i,j)=temp2/t.size();
		d1(i,j)=temp3/t.size();
		d2(i,j)=temp4/t.size();
		h1(i,j)=d1(i,j)/(1.0-CIF1(i,j));
		h2(i,j)=d2(i,j)/(1.0-CIF2(i,j));
		}
	}

	 return List::create(
	 Named("c")=c,
	 Named("nm")=nm,
	 Named("alpharec1")=alpharec1,
	 Named("lambda0rec1")=lambdarec01,
	 Named("lambdarec1")=lambdarec_orig1,
	 Named("lambdascaled1")=lambdarec1,
	 Named("alpharec2")=alpharec2,
	 Named("lambda0rec2")=lambdarec02,
	 Named("lambdarec2")=lambdarec_orig2,
	 Named("lambdascaled2")=lambdarec2,
	 Named("prec")=prec,
         Named("CIF1")=CIF1,
	 Named("CIF2")=CIF2,
         Named("d1")=d1,
	 Named("d2")=d2,
         Named("h1")=h1,
	 Named("h2")=h2,
	 Named("ngrp")=ngrp,
	 Named("emptybasket")=emptybasket,
	 Named("allbaskets")=allbaskets);
}
