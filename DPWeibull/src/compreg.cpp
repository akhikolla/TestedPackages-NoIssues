#include <RcppArmadillo.h>
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

using namespace Rcpp ;
using namespace arma;

#include "arms.h"
#include "auxfuns.h"
#include "commonfunc.h"
#include "compreggilks.h"

double compreg_loglikelihood(const double t, const int delta,
const NumericVector x,
const double lambda1, const double alpha1,
const double lambda2, const double alpha2,
const NumericVector beta1, const NumericVector beta2,
const double p){
	if (delta == 1){
		return logf1(t, p,alpha1, lambda1, x, beta1);
	}else if(delta==2){
		return logf2(t, p, alpha2, lambda2, x, beta1, beta2);
	}else{
		double xbeta1=sum(x*beta1);
		double xbeta2=sum(x*beta2);
		return logScomp(t,alpha1,lambda1,alpha2,lambda2,xbeta1,xbeta2,p);
	}
}

int compreg_group_assign(const double t,const int delta, const NumericVector x,
	const int c, const double nu,
	IntegerVector nm,
	NumericVector alpha1, NumericVector lambda1,
	NumericVector alpha2, NumericVector lambda2,
	NumericMatrix beta1, NumericMatrix beta2,
	NumericVector p,
	const double lambda00, const double alpha00, const double alpha0,
	const double alphaalpha, const double alphalambda,
	const double gamma0, const double gamma1,
	const double betasl,
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
	NumericMatrix betamatrix1(allbaskets[0]+addcluster,x.size());
	NumericMatrix betamatrix2(allbaskets[0]+addcluster,x.size());
	NumericVector pvector(allbaskets[0]+addcluster);
	NumericVector nmvector(allbaskets[0]+addcluster);

	for(int i=0;i<(allbaskets[0]+addcluster);i++){
		if(i<allbaskets[0]){
			alphavector1(i)=alpha1(i);
			lambdavector1(i)=lambda1(i);
			alphavector2(i)=alpha2(i);
			lambdavector2(i)=lambda2(i);
			betamatrix1(i,_)=beta1(i,_);
			betamatrix2(i,_)=beta2(i,_);
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
			for(int j=0;j<x.size();j++){
				samptruncauchy(&betamatrix1(i,j),betasl);
				samptruncauchy(&betamatrix2(i,j),betasl);
			}
			nmvector(i)=nu/double(m);
		}
	}
	if(nm(c-1)==0){
		nmvector(c-1)=nu/double(m);
	}
	NumericVector logprob(allbaskets[0]+addcluster);
	for (int j = 0; j<(allbaskets[0]+addcluster); j++){
			logprob[j] = compreg_loglikelihood(t,delta,x,lambdavector1[j], alphavector1[j],
			 lambdavector2[j], alphavector2[j],betamatrix1(j,_),betamatrix2(j,_),pvector[j]);
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
		beta1((emptybasket[0]-1),_)=betamatrix1(whichmax,_);
		beta2((emptybasket[0]-1),_)=betamatrix2(whichmax,_);
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
		beta1(allbaskets[0],_)=betamatrix1(whichmax,_);
		beta2(allbaskets[0],_)=betamatrix2(whichmax,_);
		allbaskets[0]=allbaskets[0]+1;
		return allbaskets[0];
	}
}

void compreg_update(NumericVector t,
	IntegerVector delta, 
	NumericMatrix x,
	IntegerVector c, IntegerVector nm,
	NumericVector alpha1, NumericVector lambda1, NumericVector lambda01, 
	NumericVector alpha2, NumericVector lambda2, NumericVector lambda02, 
	NumericMatrix beta1, NumericMatrix beta2, const double betasl,
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
	NumericMatrix xtemp=submat_rcpp(x,c==uniquec(k));
	NumericVector xbeta1=matrixtimesvector(xtemp,beta1(uniquec[k]-1,_));
	NumericVector xbeta2=matrixtimesvector(xtemp,beta2(uniquec[k]-1,_));
	
	lambda01[uniquec(k)-1]=R::rgamma(alpha0 + alpha00, 1.0/(lambda00 + lambda1[uniquec[k] - 1]));
	
	compreg_sampleloglambda1(&lambda1[uniquec[k] - 1], alpha0, lambda01[uniquec(k)-1],
	alpha1[uniquec[k] - 1],alpha2[uniquec[k] - 1],lambda2[uniquec[k] - 1],p[uniquec[k] - 1],
	ttemp.begin(),deltatemp.begin(),ttemp.size(),xbeta1.begin(),xbeta2.begin());
	
	compreg_samplealpha1(&alpha1[uniquec[k] - 1],lambda1[uniquec[k] - 1],alpha2[uniquec[k] - 1],
	lambda2[uniquec[k] - 1],p[uniquec[k] - 1],alphaalpha,alphalambda,
	ttemp.begin(),deltatemp.begin(),ttemp.size(),xbeta1.begin(),xbeta2.begin());
	
	for(int i=0;i<xtemp.ncol();i++){
			  NumericVector xsample=xtemp(_,i);
			  NumericMatrix index(clone(xtemp));
			  index(_,i)=rbinom(ttemp.size(),1,0.0);
			  NumericVector xbeta1left=matrixtimesvector(index,beta1(uniquec[k]-1,_));
	  		  compreg_samplebeta1(&beta1(uniquec[k] - 1,i), betasl,
			  alpha1[uniquec[k] - 1],lambda1[uniquec[k] - 1],
			  alpha2[uniquec[k] - 1],lambda2[uniquec[k] - 1],
			  ttemp.begin(),deltatemp.begin(),ttemp.size(),
			  xsample.begin(),
			  xbeta1left.begin(),xbeta2.begin(),p[uniquec[k] - 1]);
  	}
	xbeta1=matrixtimesvector(xtemp,beta1(uniquec[k]-1,_));
	
	lambda02[uniquec(k)-1]=R::rgamma(alpha0 + alpha00, 1.0/(lambda00 + lambda2[uniquec[k] - 1]));
	
	compreg_sampleloglambda2(&lambda2[uniquec[k] - 1], alpha0, lambda02[uniquec(k)-1],
	alpha2[uniquec[k] - 1],alpha1[uniquec[k] - 1],lambda1[uniquec[k] - 1],p[uniquec[k] - 1],
	ttemp.begin(),deltatemp.begin(),ttemp.size(),xbeta1.begin(),xbeta2.begin());
	
	compreg_samplealpha2(&alpha2[uniquec[k] - 1],lambda2[uniquec[k] - 1],alpha1[uniquec[k] - 1],
	lambda1[uniquec[k] - 1],p[uniquec[k] - 1],alphaalpha,alphalambda,
	ttemp.begin(),deltatemp.begin(),ttemp.size(),xbeta1.begin(),xbeta2.begin());
	
	for(int i=0;i<xtemp.ncol();i++){
			  NumericVector xsample=xtemp(_,i);
			  NumericMatrix index(clone(xtemp));
			  index(_,i)=rbinom(ttemp.size(),1,0.0);
			  NumericVector xbeta2left=matrixtimesvector(index,beta2(uniquec[k]-1,_));
	  		  compreg_samplebeta2(&beta2(uniquec[k] - 1,i), betasl,
			  alpha1[uniquec[k] - 1],lambda1[uniquec[k] - 1],
			  alpha2[uniquec[k] - 1],lambda2[uniquec[k] - 1],
			  ttemp.begin(),deltatemp.begin(),ttemp.size(),
			  xsample.begin(),
			  xbeta2left.begin(),xbeta1.begin(),p[uniquec[k] - 1]);
  	}
	
	compreg_samplep(&p[uniquec[k] - 1],gamma0,gamma1,ttemp.begin(),deltatemp.begin(),
	alpha1[uniquec[k] - 1],lambda1[uniquec[k] - 1],alpha2[uniquec[k] - 1],lambda2[uniquec[k] - 1],
	ttemp.size(),xbeta1.begin(),xbeta2.begin());
	}
}



List compreg(const int burnin, const int iteration,
	NumericVector t, IntegerVector delta, NumericMatrix x, IntegerVector c, IntegerVector nm, 
	NumericVector alpha1, NumericVector lambda1, NumericVector lambda01,
	NumericVector alpha2, NumericVector lambda2, NumericVector lambda02,
	NumericMatrix beta1, NumericMatrix beta2, NumericVector p,
	const double alpha00, const double alpha0, const double lambda00,
	const double alphaalpha, const double alphalambda,
	const double gamma0, const double gamma1,
	NumericVector nu, NumericVector ngrp,
	const double a, const double b,
	const double ymax, NumericVector tplot,
	int m, int thin, double betasl,NumericVector xplot1, NumericVector xplot2){
	std::vector<int> emptybasket(1);
	emptybasket[0]=0;
	IntegerVector allbaskets(1);
	allbaskets[0]=1;
	int nsave=iteration/thin;
	NumericMatrix alpharec1(nsave,t.size());
	NumericMatrix lambdarec01(nsave,t.size());
	NumericMatrix lambdarec1(nsave,t.size());
	NumericMatrix alpharec2(nsave,t.size());
	NumericMatrix lambdarec02(nsave,t.size());
	NumericMatrix lambdarec2(nsave,t.size());
	NumericMatrix lambdarec_orig1(nsave,t.size());
	NumericMatrix lambdarec_orig2(nsave,t.size());
	NumericMatrix betarec1(nsave,t.size()*x.ncol());
	NumericMatrix betarec2(nsave,t.size()*x.ncol());
	NumericMatrix prec(nsave,t.size());	
	 for (int g = 0; g<(burnin + iteration); g++){
		if((g+1)%100==0){
		Rcout<<g+1<<" iterations out of "<<burnin + iteration<<" iterations done"<<std::endl;
		}
		for (int i = 0; i<t.size(); i++){
			c[i] = compreg_group_assign(t[i], delta[i],x(i,_),c[i], nu[g], 
			nm, alpha1, lambda1,alpha2, lambda2,beta1,beta2, p,
			lambda00, alpha00, alpha0, alphaalpha,alphalambda,
			gamma0,gamma1,betasl,
			m,allbaskets,emptybasket);
		}
		
		compreg_update(t, delta,x, c, nm, alpha1, lambda1, lambda01,
		alpha2, lambda2, lambda02,
		beta1,beta2,betasl,p,
		alpha00, alpha0, lambda00,
		alphaalpha, alphalambda,
		gamma0,gamma1,
		&ngrp[g+1]);		 
		nu[g+1] = nugen(nu[g], t.size(), ngrp[g+1], a, b);

		if((g>=burnin)&&(g%thin==(thin-1))){
			int comp=(g-burnin-thin+1)/thin;
			for(int z=0;z<t.size();z++){
				alpharec1(comp,z) = alpha1(c[z] - 1);
				lambdarec01(comp,z) = lambda01(c[z] - 1);
				lambdarec1(comp,z) = lambda1(c[z] - 1);
				lambdarec_orig1(comp,z)=exp(log(lambda1(c[z] - 1))-alpha1(c[z] - 1)*(log(ymax)-log(10.0)));
				alpharec2(comp,z) = alpha2(c[z] - 1);
				lambdarec02(comp,z) = lambda02(c[z] - 1);
				lambdarec2(comp,z) = lambda2(c[z] - 1);
				lambdarec_orig2(comp,z)=exp(log(lambda2(c[z] - 1))-alpha2(c[z] - 1)*(log(ymax)-log(10.0)));
				prec(comp,z) = p(c[z] - 1);
				for(int i=0;i<x.ncol();i++){
					betarec1(comp,z*x.ncol()+i) = beta1(c[z] - 1,i);
					betarec2(comp,z*x.ncol()+i) = beta2(c[z] - 1,i);
				}
			}
		}
		
	}
	NumericMatrix loghr(nsave,tplot.size()*x.ncol());
    for(int covnum=0; covnum<x.ncol(); covnum++){   
	  NumericMatrix S1(nsave,tplot.size());
	  NumericMatrix d1(nsave,tplot.size());
	  NumericMatrix h1(nsave,tplot.size());
	  NumericMatrix S2(nsave,tplot.size());
	  NumericMatrix d2(nsave,tplot.size());
	  NumericMatrix h2(nsave,tplot.size());
 	       for(int i=0;i<nsave;i++){
		 for(int j=0;j<tplot.size();j++){
		double temps1=0.0;
		double temps2=0.0;
		double tempd1=0.0;
		double tempd2=0.0;
			for(int k=0;k<t.size();k++){
			double temptemp1=xplot1(covnum)*betarec1(i,k*x.ncol()+covnum);
			double temptemp2=xplot2(covnum)*betarec1(i,k*x.ncol()+covnum);
			temps1+=1.0-F1v2(tplot(j),prec(i,k),alpharec1(i,k),lambdarec1(i,k),temptemp1);
			temps2+=1.0-F1v2(tplot(j),prec(i,k),alpharec1(i,k),lambdarec1(i,k),temptemp2);
			tempd1+=f1v2(tplot(j),prec(i,k),alpharec1(i,k),lambdarec1(i,k),temptemp1);
			tempd2+=f1v2(tplot(j),prec(i,k),alpharec1(i,k),lambdarec1(i,k),temptemp2);
			}
		S1(i,j)=temps1/t.size();
		d1(i,j)=tempd1/t.size();
		h1(i,j)=d1(i,j)/S1(i,j);
		S2(i,j)=temps2/t.size();
		d2(i,j)=tempd2/t.size();
		h2(i,j)=d2(i,j)/S2(i,j);
	    loghr(i,j+covnum*tplot.size())=log(d1(i,j))-log(S1(i,j))-log(d2(i,j))+log(S2(i,j));	
		}
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
	 Named("betarec1")=betarec1,
	 Named("betarec2")=betarec2,
	 Named("prec")=prec,
     Named("loghr")=loghr,
	 Named("ngrp")=ngrp,
	 Named("emptybasket")=emptybasket,
	 Named("allbaskets")=allbaskets);
}


List compreg_resume(const int burnin, const int iteration,
	NumericVector t, IntegerVector delta, NumericMatrix x,	IntegerVector c, IntegerVector nm, 
	NumericVector alpha1, NumericVector lambda1, NumericVector lambda01,
	NumericVector alpha2, NumericVector lambda2, NumericVector lambda02,
	NumericMatrix beta1, NumericMatrix beta2, NumericVector p,
	const double alpha00, const double alpha0, const double lambda00,
	const double alphaalpha, const double alphalambda,
	const double gamma0, const double gamma1,
	NumericVector nu, NumericVector ngrp,
	const double a, const double b,
	const double ymax, NumericVector tplot,
	int m, int thin, double betasl, NumericVector xplot1, NumericVector xplot2,
	std::vector<int> emptybasket, IntegerVector allbaskets){
	int nsave=iteration/thin;
	NumericMatrix alpharec1(nsave,t.size());
	NumericMatrix lambdarec01(nsave,t.size());
	NumericMatrix lambdarec1(nsave,t.size());
	NumericMatrix alpharec2(nsave,t.size());
	NumericMatrix lambdarec02(nsave,t.size());
	NumericMatrix lambdarec2(nsave,t.size());
	NumericMatrix lambdarec_orig1(nsave,t.size());
	NumericMatrix lambdarec_orig2(nsave,t.size());
	NumericMatrix betarec1(nsave,t.size()*x.ncol());
	NumericMatrix betarec2(nsave,t.size()*x.ncol());
	NumericMatrix prec(nsave,t.size());	
	 for (int g = 0; g<(burnin + iteration); g++){
		if((g+1)%100==0){
		Rcout<<g+1<<" iterations out of "<<burnin + iteration<<" iterations done"<<std::endl;
		}
		for (int i = 0; i<t.size(); i++){
			c[i] = compreg_group_assign(t[i], delta[i],x(i,_),c[i], nu[g], 
			nm, alpha1, lambda1,alpha2, lambda2,beta1,beta2, p,
			lambda00, alpha00, alpha0, alphaalpha,alphalambda,
			gamma0,gamma1,betasl,
			m,allbaskets,emptybasket);
		}
		
		 compreg_update(t, delta,x, c, nm, alpha1, lambda1, lambda01,
		 alpha2, lambda2, lambda02,
		 beta1,beta2,betasl,p,
		 alpha00, alpha0, lambda00,
		 alphaalpha, alphalambda,
		 gamma0,gamma1,
		 &ngrp[g+1]);		 
		nu[g+1] = nugen(nu[g], t.size(), ngrp[g+1], a, b);

		if((g>=burnin)&&(g%thin==(thin-1))){
			int comp=(g-burnin-thin+1)/thin;
			for(int z=0;z<t.size();z++){
				alpharec1(comp,z) = alpha1(c[z] - 1);
				lambdarec01(comp,z) = lambda01(c[z] - 1);
				lambdarec1(comp,z) = lambda1(c[z] - 1);
				lambdarec_orig1(comp,z)=exp(log(lambda1(c[z] - 1))-alpha1(c[z] - 1)*(log(ymax)-log(10.0)));
				alpharec2(comp,z) = alpha2(c[z] - 1);
				lambdarec02(comp,z) = lambda02(c[z] - 1);
				lambdarec2(comp,z) = lambda2(c[z] - 1);
				lambdarec_orig2(comp,z)=exp(log(lambda2(c[z] - 1))-alpha2(c[z] - 1)*(log(ymax)-log(10.0)));
				prec(comp,z) = p(c[z] - 1);
				for(int i=0;i<x.ncol();i++){
					betarec1(comp,z*x.ncol()+i) = beta1(c[z] - 1,i);
					betarec2(comp,z*x.ncol()+i) = beta2(c[z] - 1,i);
				}
			}
		}	
	}
	
	NumericMatrix loghr(nsave,tplot.size()*x.ncol());
    for(int covnum=0; covnum<x.ncol(); covnum++){   
		NumericMatrix S1(nsave,tplot.size());
		NumericMatrix d1(nsave,tplot.size());
		NumericMatrix h1(nsave,tplot.size());
		NumericMatrix S2(nsave,tplot.size());
		NumericMatrix d2(nsave,tplot.size());
		NumericMatrix h2(nsave,tplot.size());
 	    for(int i=0;i<nsave;i++){
			for(int j=0;j<tplot.size();j++){
			double temps1=0.0;
			double temps2=0.0;
			double tempd1=0.0;
			double tempd2=0.0;
			for(int k=0;k<t.size();k++){
				double temptemp1=xplot1(covnum)*betarec1(i,k*x.ncol()+covnum);
				double temptemp2=xplot2(covnum)*betarec1(i,k*x.ncol()+covnum);
				temps1+=1.0-F1v2(tplot(j),prec(i,k),alpharec1(i,k),lambdarec1(i,k),temptemp1);
				temps2+=1.0-F1v2(tplot(j),prec(i,k),alpharec1(i,k),lambdarec1(i,k),temptemp2);
				tempd1+=f1v2(tplot(j),prec(i,k),alpharec1(i,k),lambdarec1(i,k),temptemp1);
				tempd2+=f1v2(tplot(j),prec(i,k),alpharec1(i,k),lambdarec1(i,k),temptemp2);
			}
			S1(i,j)=temps1/t.size();
			d1(i,j)=tempd1/t.size();
			h1(i,j)=d1(i,j)/S1(i,j);
			S2(i,j)=temps2/t.size();
			d2(i,j)=tempd2/t.size();
			h2(i,j)=d2(i,j)/S2(i,j);
	    	loghr(i,j+covnum*tplot.size())=log(d1(i,j))-log(S1(i,j))-log(d2(i,j))+log(S2(i,j));	
			}
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
	 Named("betarec1")=betarec1,
	 Named("betarec2")=betarec2,
	 Named("prec")=prec,
     Named("loghr")=loghr,
	 Named("ngrp")=ngrp,
	 Named("emptybasket")=emptybasket,
	 Named("allbaskets")=allbaskets);
}



List predcompreg(NumericMatrix alpharec1,NumericMatrix lambdarec1,NumericMatrix betarec1,
NumericMatrix alpharec2,NumericMatrix lambdarec2,NumericMatrix betarec2,
NumericMatrix prec, NumericMatrix xplot, NumericVector tplot, double alpha){
    int nsave=alpharec1.nrow();
    NumericMatrix Fpred(xplot.nrow(),tplot.size());    
    NumericMatrix Fpredl(xplot.nrow(),tplot.size());  
    NumericMatrix Fpredu(xplot.nrow(),tplot.size());  
    NumericMatrix dpred(xplot.nrow(),tplot.size());  
    NumericMatrix dpredl(xplot.nrow(),tplot.size());  
    NumericMatrix dpredu(xplot.nrow(),tplot.size());  
    NumericMatrix hpred(xplot.nrow(),tplot.size());  
    NumericMatrix hpredl(xplot.nrow(),tplot.size());  
    NumericMatrix hpredu(xplot.nrow(),tplot.size());  
	cube F(nsave,tplot.size(),xplot.nrow());
	cube d(nsave,tplot.size(),xplot.nrow());
	cube h(nsave,tplot.size(),xplot.nrow());
	for(int covnum=0; covnum<xplot.nrow(); covnum++){   
 	    for(int i=0;i<nsave;i++){
			for(int j=0;j<tplot.size();j++){
			double tempF=0.0;
			double tempd=0.0;
				for(int k=0;k<alpharec1.ncol();k++){
		            double temptemp=0.0;
		               	for(int l=0; l<xplot.ncol();l++){
							temptemp+=xplot(covnum,l)*betarec1(i,k*xplot.ncol()+l);
						}
					tempF+=F1v2(tplot(j),prec(i,k),alpharec1(i,k),lambdarec1(i,k),temptemp);
					tempd+=f1v2(tplot(j),prec(i,k),alpharec1(i,k),lambdarec1(i,k),temptemp);
				}
			F(i,j,covnum)=tempF/alpharec1.ncol();
			d(i,j,covnum)=tempd/alpharec1.ncol();
			h(i,j,covnum)=d(i,j,covnum)/(1.0-F(i,j,covnum));
		}
	   }
	Fpred(covnum,_)=colpercentileRcpp(wrap(F.slice(covnum)),0.5);
	Fpredu(covnum,_)=colpercentileRcpp(wrap(F.slice(covnum)),1.0-alpha/2.0);
	Fpredl(covnum,_)=colpercentileRcpp(wrap(F.slice(covnum)),alpha/2.0);
	dpred(covnum,_)=colpercentileRcpp(wrap(d.slice(covnum)),0.5);
	dpredu(covnum,_)=colpercentileRcpp(wrap(d.slice(covnum)),1.0-alpha/2.0);
	dpredl(covnum,_)=colpercentileRcpp(wrap(d.slice(covnum)),alpha/2.0);
	hpred(covnum,_)=colpercentileRcpp(wrap(h.slice(covnum)),0.5);
	hpredu(covnum,_)=colpercentileRcpp(wrap(h.slice(covnum)),1.0-alpha/2.0);
	hpredl(covnum,_)=colpercentileRcpp(wrap(h.slice(covnum)),alpha/2.0);
	}
	return List::create(
	Named("F")=F,
	Named("d")=d,
	Named("h")=h,
	Named("Fpred")=Fpred,
	Named("Fpredl")=Fpredl,
	Named("Fpredu")=Fpredu,
	Named("dpred")=dpred,
	Named("dpredl")=dpredl,
	Named("dpredu")=dpredu,
	Named("hpred")=hpred,
	Named("hpredl")=hpredl,
	Named("hpredu")=hpredu);
}
