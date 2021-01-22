#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <Rmath.h>
#include <R.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp ;
using namespace arma;

#include "arms.h"
#include "auxfuns.h"
#include "commonfunc.h"
#include "noreggilks.h"
#include "noreg.h"
#include "reggilks.h"

double reg_loglikelihood(const double tl, const double tr,
const int delta, const int pi,
const double lambda, const double alpha,
NumericVector x, NumericVector beta){
	double loglambda_xbeta=log(lambda)+sum(x*beta);
	return noreg_loglikelihood(tl,tr,delta,pi, loglambda_xbeta, alpha);
}

int reg_group_assign(const double tl, const double tr,
	const int delta,const int pi,
	NumericVector x,
	const int c, const double nu,
	IntegerVector nm,
	NumericVector alpha, NumericVector lambda,
	NumericMatrix beta,
	const double lambda00, const double alpha00, const double alpha0,
	const double alphaalpha, const double alphalambda,
	const double betasl,
	const int m,
	IntegerVector allbaskets, std::vector<int> & emptybasket){
	nm(c-1)=nm(c-1)-1;
	int addcluster=m;
	if(nm(c-1)==0){
		emptybasket.insert(emptybasket.begin(),c);
		addcluster=m-1;
	}
	NumericVector alphavector(allbaskets[0]+addcluster);
	NumericVector lambdavector(allbaskets[0]+addcluster);
	NumericMatrix betamatrix(allbaskets[0]+addcluster,x.size());
	NumericVector nmvector(allbaskets[0]+addcluster);
	double base_new=0.0;
	double temp1=0.0;
	for(int i=0;i<(allbaskets[0]+addcluster);i++){
		if(i<allbaskets[0]){
			alphavector(i)=alpha(i);
			lambdavector(i)=lambda(i);
			betamatrix(i,_)=beta(i,_);
			nmvector(i)=double(nm(i));
		}else{
			double lambda0_new=R::rgamma(alpha00,1.0/lambda00);
			lambdavector(i) = R::rgamma(alpha0,1.0/ lambda0_new);
			for(int j=0;j<x.size();j++){
					samptruncauchy(&betamatrix(i,j),betasl);
				}
			base_new=findbase(lambdavector(i));
			if(base_new<80.0){
				temp1=R::pgamma(base_new,alphaalpha,1.0/alphalambda,1,0);
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
			logprob[j] = reg_loglikelihood(tl,tr,delta,pi,lambdavector[j], alphavector[j],x,betamatrix(j,_));
		//if(!testreal(logprob[j])){logprob[j]=-pow(100.0,100.0);}
	}
        logprob=logprob-max(logprob);
		//		Rcout<<"logprob	"<<logprob<<std::endl;
	NumericVector prob=exp(logprob)*nmvector;
		//		Rcout<<"prob	"<<prob<<std::endl;
	  if(sum(prob)>0.0){
		  prob=prob/sum(prob);
	  }else{
		  prob=nmvector;
		  prob=prob/sum(prob);
	  }
    IntegerVector newassign=oneMultinomCalt(prob);
	int whichmax=which_max(newassign);

	if((whichmax+1)<=allbaskets[0]){
		nm[whichmax]=nm[whichmax]+1;
		if(whichmax+1==emptybasket[0]){
			emptybasket.erase(emptybasket.begin());
		}
		return whichmax+1;
	}else if(((whichmax+1)>allbaskets[0])&&(emptybasket[0]!=0)){
		alpha[emptybasket[0]-1]=alphavector(whichmax);
		lambda[emptybasket[0]-1]=lambdavector(whichmax);
		beta((emptybasket[0]-1),_)=betamatrix(whichmax,_);
		int fill=emptybasket[0];
		nm[emptybasket[0]-1]=1;
		emptybasket.erase(emptybasket.begin());
		return fill;
	}else{
		nm[allbaskets[0]]=1;
		alpha[allbaskets[0]]=alphavector(whichmax);
		lambda[allbaskets[0]]=lambdavector(whichmax);
		beta(allbaskets[0],_)=betamatrix(whichmax,_);
		allbaskets[0]=allbaskets[0]+1;
		return allbaskets[0];
	}
}

void reg_update(NumericVector tl, NumericVector tr,
	IntegerVector delta, IntegerVector pi,  NumericMatrix x,
	IntegerVector c, IntegerVector nm,
	NumericVector alpha, NumericVector lambda, 
	NumericMatrix beta,
	NumericVector lambda0, 
	const double alpha00, const double alpha0, const double lambda00,
	const double alphaalpha, const double alphalambda,
    const double betasl,double* nextngrp){
	IntegerVector	uniquec=unique(c);
	*nextngrp=double(uniquec.size());
	for(int k=0;k<*nextngrp;k++){
	NumericVector tltemp=tl[c==uniquec(k)];
	NumericVector trtemp=tr[c==uniquec(k)];
	IntegerVector deltatemp=delta[c==uniquec(k)];
	IntegerVector pitemp=pi[c==uniquec(k)];
	NumericMatrix xtemp=submat_rcpp(x,c==uniquec(k));
	lambda0[uniquec(k)-1]=R::rgamma(alpha0 + alpha00, 1.0/(lambda00 + lambda[uniquec[k] - 1]));
	NumericVector xbeta=matrixtimesvector(xtemp,beta(uniquec[k]-1,_));
  	reg_sampleloglambda(&lambda[uniquec[k] - 1], alpha[uniquec[k] - 1],
	alpha0, lambda0[uniquec(k)-1], tltemp.begin(),trtemp.begin(),
	deltatemp.begin(),pitemp.begin(),
	tltemp.size(),xbeta.begin());
	reg_samplealpha(&alpha[uniquec[k] - 1], log(lambda[uniquec[k] - 1]),
	alphaalpha, alphalambda,tltemp.begin(),trtemp.begin(),deltatemp.begin(),pitemp.begin(),
	tltemp.size(), xbeta.begin());
	for(int i=0;i<xtemp.ncol();i++){
			  NumericVector xsample=xtemp(_,i);
			  NumericMatrix emptycol(clone(xtemp));
			  emptycol(_,i)=rbinom(tltemp.size(),1,0.0);
			  NumericVector partloglambda=rep(log(lambda[uniquec[k] - 1]),tltemp.size())
			  +matrixtimesvector(emptycol,beta(uniquec[k]-1,_));
	  		  reg_samplebeta(&beta(uniquec[k] - 1,i), betasl,
			  alpha[uniquec[k] - 1],tltemp.begin(),trtemp.begin(), 
			  deltatemp.begin(),pitemp.begin(),tltemp.size(),
			  xsample.begin(), partloglambda.begin());
  	}
 }
}	


List reg(const int burnin, const int iteration,
	NumericVector tl, NumericVector tr,
	IntegerVector delta,
	IntegerVector pi,
	NumericMatrix x, 
	IntegerVector c,
  	IntegerVector nm, 
	NumericVector alpha,
	NumericVector lambda,
	NumericMatrix beta,
	NumericVector lambda0,
	const double alpha00,
	const double alpha0,
	const double lambda00,
	const double alphaalpha,
	const double alphalambda,
	NumericVector nu,
	NumericVector ngrp,
	const double a, const double b,
	const double ymax, int m,
	double betasl, NumericVector time, NumericVector xpred1, NumericVector xpred2, int thin){
	std::vector<int> emptybasket(1);
	emptybasket[0]=0;
	IntegerVector allbaskets(1);
	allbaskets[0]=1;
	int nsave=iteration/thin;
	NumericMatrix alpharec(nsave,tl.size());
	NumericMatrix lambda0rec(nsave,tl.size());
	NumericMatrix lambdarec(nsave,tl.size());
	NumericMatrix lambdaorigrec(nsave,tl.size());
	NumericMatrix betarec(nsave,tl.size()*x.ncol());
	for (int g = 0; g<(burnin + iteration); g++){
		if((g+1)%100==0){
		Rcout<<g+1<<" iterations out of "<<burnin + iteration<<" iterations done"<<std::endl;
		}
		for (int i = 0; i<tl.size(); i++){
			c[i] =reg_group_assign(tl[i], tr[i], delta[i], pi[i], x(i,_), c[i], nu[g], nm, alpha, lambda,beta,lambda00, alpha00, alpha0, alphaalpha,alphalambda,betasl,m,allbaskets,emptybasket);
		}
		reg_update(tl, tr, delta, pi, x, c, nm, alpha, lambda, beta, lambda0,
		alpha00, alpha0, lambda00, alphaalpha, alphalambda, betasl, &ngrp[g+1]);
 		nu[g+1] = nugen(nu[g], tl.size(), ngrp[g+1], a, b);
	if((g>=burnin)&&(g%thin==(thin-1))){
			int index=(g-burnin-(thin-1))/thin;
                for(int z=0;z<tl.size();z++){
				alpharec(index,z) = alpha(c[z] - 1);
				lambda0rec(index,z) = lambda0(c[z] - 1);
				lambdarec(index,z) = lambda(c[z] - 1);
				lambdaorigrec(index,z)=exp(log(lambda(c[z] - 1))-alpha(c[z] - 1)*(log(ymax)-log(10.0)));
				for(int i=0;i<x.ncol();i++){
					betarec(index,z*x.ncol()+i) = beta(c[z] - 1,i);
				}
			}
		}
	}
          NumericMatrix loghr(nsave,time.size()*x.ncol());

for(int covnum=0; covnum<x.ncol();covnum++){
	  NumericMatrix S1(nsave,time.size());
	  NumericMatrix d1(nsave,time.size());
	  NumericMatrix h1(nsave,time.size());
	  NumericMatrix S2(nsave,time.size());
	  NumericMatrix d2(nsave,time.size());
	  NumericMatrix h2(nsave,time.size());
       	  for(int i=0;i<nsave;i++){
		for(int j=0;j<time.size();j++){
		double temps1=0.0;
		double temps2=0.0;
		double tempd1=0.0;
		double tempd2=0.0;
		for(int k=0;k<tl.size();k++){
		double temptemp1=xpred1(covnum)*betarec(i,k*x.ncol()+covnum);
		double temptemp2=xpred2(covnum)*betarec(i,k*x.ncol()+covnum);
		temps1+=sWeibloglambda(time[j],alpharec(i,k),log(lambdarec(i,k))+temptemp1);
		tempd1+=dWeibloglambda(time[j],alpharec(i,k),log(lambdarec(i,k))+temptemp1);
		temps2+=sWeibloglambda(time[j],alpharec(i,k),log(lambdarec(i,k))+temptemp2);
		tempd2+=dWeibloglambda(time[j],alpharec(i,k),log(lambdarec(i,k))+temptemp2);
		}
		S1(i,j)=temps1/tl.size();
		d1(i,j)=tempd1/tl.size();
		h1(i,j)=d1(i,j)/S1(i,j);
		S2(i,j)=temps2/tl.size();
		d2(i,j)=tempd2/tl.size();
		h2(i,j)=d2(i,j)/S2(i,j);
	    loghr(i,j+covnum*time.size())=log(d1(i,j))-log(S1(i,j))
			-log(d2(i,j))+log(S2(i,j));				
		}
	  }
}
	return List::create(
			 Named("nu")=nu,
	 Named("c")=c,
	 Named("nm")=nm,
	Named("loghr")=loghr,
    	Named("betarec")=betarec,
	Named("alpharec")=alpharec,
	Named("lambda0rec")=lambda0rec,
	Named("lambdascaled")=lambdarec,
	Named("lambdarec")=lambdaorigrec,
	Named("ngrp")=ngrp,
        Named("allbaskets")=allbaskets,
        Named("emptybasket")=emptybasket);
}


List reg_resume(const int burnin, const int iteration,
	NumericVector tl, NumericVector tr,
	IntegerVector delta,
	IntegerVector pi,
	NumericMatrix x, 
	IntegerVector c,
  	IntegerVector nm, 
	NumericVector alpha,
	NumericVector lambda,
	NumericMatrix beta,
	NumericVector lambda0,
	const double alpha00,
	const double alpha0,
	const double lambda00,
	const double alphaalpha,
	const double alphalambda,
	NumericVector nu,
	NumericVector ngrp,
	const double a, const double b,
	const double ymax, int m,
	double betasl, NumericVector time, NumericVector xpred1, NumericVector xpred2, int thin,
        std::vector<int> emptybasket, IntegerVector allbaskets){

	int nsave=iteration/thin;
	NumericMatrix alpharec(nsave,tl.size());
	NumericMatrix lambda0rec(nsave,tl.size());
	NumericMatrix lambdarec(nsave,tl.size());
	NumericMatrix lambdaorigrec(nsave,tl.size());
	NumericMatrix betarec(nsave,tl.size()*x.ncol());
	for (int g = 0; g<(burnin + iteration); g++){
		if((g+1)%100==0){
		Rcout<<g+1<<" iterations out of "<<burnin + iteration<<" iterations done"<<std::endl;
		}
		for (int i = 0; i<tl.size(); i++){
			c[i] =reg_group_assign(tl[i], tr[i], delta[i], pi[i], x(i,_), c[i], nu[g], nm, alpha, lambda,beta,
     			 lambda00, alpha00, alpha0, alphaalpha,alphalambda,betasl,m,allbaskets,emptybasket);
		}
		 reg_update(tl, tr, delta, pi, x, c, nm, alpha, lambda, beta, lambda0,
		 alpha00, alpha0, lambda00, alphaalpha, alphalambda, betasl, &ngrp[g+1]);
 		nu[g+1] = nugen(nu[g], tl.size(), ngrp[g+1], a, b);
	if((g>=burnin)&&(g%thin==(thin-1))){
			int index=(g-burnin-(thin-1))/thin;
                for(int z=0;z<tl.size();z++){
				alpharec(index,z) = alpha(c[z] - 1);
				lambda0rec(index,z) = lambda0(c[z] - 1);
				lambdarec(index,z) = lambda(c[z] - 1);
				lambdaorigrec(index,z)=exp(log(lambda(c[z] - 1))-alpha(c[z] - 1)*(log(ymax)-log(10.0)));
				for(int i=0;i<x.ncol();i++){
					betarec(index,z*x.ncol()+i) = beta(c[z] - 1,i);
				}
			}
		}
	}
          NumericMatrix loghr(nsave,time.size()*x.ncol());

for(int covnum=0; covnum<x.ncol();covnum++){
	  NumericMatrix S1(nsave,time.size());
	  NumericMatrix d1(nsave,time.size());
	  NumericMatrix h1(nsave,time.size());
	  NumericMatrix S2(nsave,time.size());
	  NumericMatrix d2(nsave,time.size());
	  NumericMatrix h2(nsave,time.size());
       	  for(int i=0;i<nsave;i++){
		for(int j=0;j<time.size();j++){
		double temps1=0.0;
		double temps2=0.0;
		double tempd1=0.0;
		double tempd2=0.0;
		for(int k=0;k<tl.size();k++){
		double temptemp1=xpred1(covnum)*betarec(i,k*x.ncol()+covnum);
		double temptemp2=xpred2(covnum)*betarec(i,k*x.ncol()+covnum);
		temps1+=sWeibloglambda(time[j],alpharec(i,k),log(lambdarec(i,k))+temptemp1);
		tempd1+=dWeibloglambda(time[j],alpharec(i,k),log(lambdarec(i,k))+temptemp1);
		temps2+=sWeibloglambda(time[j],alpharec(i,k),log(lambdarec(i,k))+temptemp2);
		tempd2+=dWeibloglambda(time[j],alpharec(i,k),log(lambdarec(i,k))+temptemp2);
		}
		S1(i,j)=temps1/tl.size();
		d1(i,j)=tempd1/tl.size();
		h1(i,j)=d1(i,j)/S1(i,j);
		S2(i,j)=temps2/tl.size();
		d2(i,j)=tempd2/tl.size();
		h2(i,j)=d2(i,j)/S2(i,j);
	    loghr(i,j+covnum*time.size())=log(d1(i,j))-log(S1(i,j))
			-log(d2(i,j))+log(S2(i,j));		
		}
	  }
}
	return List::create(
	 Named("nu")=nu,
	 Named("c")=c,
	 Named("nm")=nm,
	Named("loghr")=loghr,
    	Named("betarec")=betarec,
	Named("alpharec")=alpharec,
	Named("lambda0rec")=lambda0rec,
	Named("lambdascaled")=lambdarec,
	Named("lambdarec")=lambdaorigrec,
	Named("ngrp")=ngrp,
        Named("allbaskets")=allbaskets,
        Named("emptybasket")=emptybasket);
}

List predreg(NumericMatrix alpharec,
       NumericMatrix lambdarec,
       NumericMatrix betarec,
       NumericMatrix xplot, NumericVector tplot, double alpha){
       int nsave=alpharec.nrow();
    NumericMatrix Spred(xplot.nrow(),tplot.size());    
    NumericMatrix Spredl(xplot.nrow(),tplot.size());  
    NumericMatrix Spredu(xplot.nrow(),tplot.size());  
    NumericMatrix dpred(xplot.nrow(),tplot.size());  
    NumericMatrix dpredl(xplot.nrow(),tplot.size());  
    NumericMatrix dpredu(xplot.nrow(),tplot.size());  
    NumericMatrix hpred(xplot.nrow(),tplot.size());  
    NumericMatrix hpredl(xplot.nrow(),tplot.size());  
    NumericMatrix hpredu(xplot.nrow(),tplot.size());  
	cube S(nsave,tplot.size(),xplot.nrow());
	cube d(nsave,tplot.size(),xplot.nrow());
	cube h(nsave,tplot.size(),xplot.nrow());
	for(int covnum=0; covnum<xplot.nrow(); covnum++){   
 	    for(int i=0;i<nsave;i++){
			for(int j=0;j<tplot.size();j++){
			double tempS=0.0;
			double tempd=0.0;
				for(int k=0;k<alpharec.ncol();k++){
		            double temptemp=0.0;
		                for(int l=0; l<xplot.ncol();l++){
							temptemp+=xplot(covnum,l)*betarec(i,k*xplot.ncol()+l);
					}
					tempS+=sWeibloglambda(tplot(j),alpharec(i,k),log(lambdarec(i,k))+temptemp);
					tempd+=dWeibloglambda(tplot(j),alpharec(i,k),log(lambdarec(i,k))+temptemp);
				}
			S(i,j,covnum)=tempS/alpharec.ncol();
			d(i,j,covnum)=tempd/alpharec.ncol();
			h(i,j,covnum)=d(i,j,covnum)/S(i,j,covnum);
		}
	   }
	Spred(covnum,_)=colpercentileRcpp(wrap(S.slice(covnum)),0.5);
	Spredu(covnum,_)=colpercentileRcpp(wrap(S.slice(covnum)),1.0-alpha/2.0);
	Spredl(covnum,_)=colpercentileRcpp(wrap(S.slice(covnum)),alpha/2.0);
	dpred(covnum,_)=colpercentileRcpp(wrap(d.slice(covnum)),0.5);
	dpredu(covnum,_)=colpercentileRcpp(wrap(d.slice(covnum)),1.0-alpha/2.0);
	dpredl(covnum,_)=colpercentileRcpp(wrap(d.slice(covnum)),alpha/2.0);
	hpred(covnum,_)=colpercentileRcpp(wrap(h.slice(covnum)),0.5);
	hpredu(covnum,_)=colpercentileRcpp(wrap(h.slice(covnum)),1.0-alpha/2.0);
	hpredl(covnum,_)=colpercentileRcpp(wrap(h.slice(covnum)),alpha/2.0);
	}
	 return List::create(
	Named("S")=S,
	Named("d")=d,
	Named("h")=h,	 
	Named("Spred")=Spred,
	Named("Spredl")=Spredl,
	Named("Spredu")=Spredu,
	Named("dpred")=dpred,
	Named("dpredl")=dpredl,
	Named("dpredu")=dpredu,
	Named("hpred")=hpred,
	Named("hpredl")=hpredl,
	Named("hpredu")=hpredu);
}


