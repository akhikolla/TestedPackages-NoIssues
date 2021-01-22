#include <iostream>
#include <vector>
#include <Rmath.h>
#include <R.h>
#include <algorithm>
#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;
#include "arms.h"
#include "commonfunc.h"
#include "auxfuns.h"
#include "compnoreg.h"

struct compreg_loglambda1_parm{
	double alpha1,lambda01, alpha2,  lambda2,alpha0, p;
	double *t;
	int *delta;
	double *xbeta1;
	double *xbeta2;
	int size;
};

double compreg_logdloglambda1(double loglambda, void* loglambdacomp1_data){
struct compreg_loglambda1_parm *d = (struct compreg_loglambda1_parm*)loglambdacomp1_data;
double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
		if((testreal(temp))&&(abs(temp)>=0.0)){
			 if(d->delta[i]==1){
				 temp+=logf1v2(d->t[i],d->p,d->alpha1, exp(loglambda), d->xbeta1[i]);
			 }else if(d->delta[i]==2){
				 temp=temp*1.0;
			 }else{
				 temp+=logScomp(d->t[i],d->alpha1,exp(loglambda),
				 d->alpha2,d->lambda2,d->xbeta1[i],d->xbeta2[i],d->p);
			 }
		}
	}
	return (d->alpha0)*loglambda-(d->lambda01)*exp(loglambda)+temp;
}

void compreg_sampleloglambda1(double* lambda,
double alpha0, double lambda01,
 double alpha1, 
 double alpha2,double lambda2,  double p,
 double* t,
 int* delta, int size, double* xbeta1, double* xbeta2)
{ 
  double loglambda=log(*lambda);
  double xl=log(inversebase(alpha1));
  double xr=1.0;
  int ninit=5;
  int dometrop = 1;
  double xprev = loglambda;
  struct compreg_loglambda1_parm loglambdacomp1_data;
  loglambdacomp1_data.alpha0=alpha0;
  loglambdacomp1_data.lambda01=lambda01;
  loglambdacomp1_data.alpha1=alpha1;
  loglambdacomp1_data.alpha2=alpha2;
  loglambdacomp1_data.lambda2=lambda2;
  loglambdacomp1_data.p=p;
  loglambdacomp1_data.size=size;
  loglambdacomp1_data.t=t;
  loglambdacomp1_data.delta=delta;
  loglambdacomp1_data.xbeta1=xbeta1;
  loglambdacomp1_data.xbeta2=xbeta2;
  int err=arms_simple(ninit, &xl, &xr,compreg_logdloglambda1,&loglambdacomp1_data,dometrop,&xprev, &loglambda);
  if(testreal(loglambda)|(err!=0)){
	  *lambda=exp(loglambda);
  }
}

struct compreg_loglambda2_parm{
	double alpha1,lambda1,lambda02, alpha2, alpha0, p;
	double *t;
	int *delta;
	double *xbeta1;
	double *xbeta2;
	int size;
};

double compreg_logdloglambda2(double loglambda, void* loglambdacomp2_data){
struct compreg_loglambda2_parm *d = (struct compreg_loglambda2_parm*)loglambdacomp2_data;
double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
		if((testreal(temp))&&(abs(temp)>=0.0)){
			 if(d->delta[i]==1){
				 temp=temp*1.0;
			 }else if(d->delta[i]==2){
				 temp+=logf2v2(d->t[i],d->p,d->alpha2,exp(loglambda), d->xbeta1[i],d->xbeta2[i]);
			 }else{
				 temp+=logScomp(d->t[i],d->alpha1,d->lambda1,d->alpha2,
				 exp(loglambda),d->xbeta1[i],d->xbeta2[i],d->p);
			 }
		  }
		}

	return (d->alpha0)*loglambda-(d->lambda02)*exp(loglambda)+temp;
}

void compreg_sampleloglambda2(double* lambda,
double alpha0, double lambda02,
 double alpha2, 
 double alpha1,double lambda1,  double p,
 double* t,
 int* delta, int size, double* xbeta1, double* xbeta2)
{ 
  double loglambda=log(*lambda);
  double xl=log(inversebase(alpha2));
  double xr=1.0;
  int ninit=5;
  int dometrop = 1;
  double xprev = loglambda;
  struct compreg_loglambda2_parm loglambdacomp2_data;
  loglambdacomp2_data.alpha0=alpha0;
  loglambdacomp2_data.lambda02=lambda02;
  loglambdacomp2_data.alpha1=alpha1;
  loglambdacomp2_data.alpha2=alpha2;
  loglambdacomp2_data.lambda1=lambda1;
  loglambdacomp2_data.p=p;
  loglambdacomp2_data.size=size;
  loglambdacomp2_data.t=t;
  loglambdacomp2_data.delta=delta;
  loglambdacomp2_data.xbeta1=xbeta1;
  loglambdacomp2_data.xbeta2=xbeta2;
  int err=arms_simple(ninit, &xl, &xr,compreg_logdloglambda2,&loglambdacomp2_data,dometrop,&xprev, &loglambda);
  if(testreal(loglambda)|(err!=0)){
	  *lambda=exp(loglambda);
 }
}


struct compreg_alpha1_parm{
	double lambda1,alpha2,lambda2, alphaalpha, alphalambda, p;
	double *t;
	int *delta;
	double *xbeta1;
	double *xbeta2;
	int size;
};

double compreg_logdalpha1(double alpha, void* alpha1_data){
	struct compreg_alpha1_parm *d = (struct compreg_alpha1_parm*)alpha1_data;
	double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
		if((testreal(temp))&&(abs(temp)>=0.0)){
			 if(d->delta[i]==1){
				 temp+=logf1v2(d->t[i],d->p,alpha,d->lambda1,d->xbeta1[i]);
			 }else if(d->delta[i]==2){
				 temp=temp*1.0;
			 }else{
				 temp+=logScomp(d->t[i],alpha,d->lambda1,
				 d->alpha2,d->lambda2,d->xbeta1[i],d->xbeta2[i],d->p);
			 }
		}
	}
	return (d->alphaalpha-1.0)*log(alpha)-(d->alphalambda)*alpha+temp;
}



void compreg_samplealpha1(double* alpha, double lambda1,
double alpha2, double lambda2,
double p,
double alphaalpha, double alphalambda,
 double* t,
 int* delta, int size, double* xbeta1, double* xbeta2)
{
  double xl=findbase(lambda1);
  double xr=80.0;
  int ninit=5;
  int dometrop = 1;
  double xprev = *alpha;
  struct compreg_alpha1_parm alpha1_data;
  alpha1_data.alphaalpha=alphaalpha;
  alpha1_data.alphalambda=alphalambda;
  alpha1_data.alpha2=alpha2;
  alpha1_data.lambda1=lambda1;
  alpha1_data.lambda2=lambda2;
  alpha1_data.size=size;
  alpha1_data.t=t;
  alpha1_data.delta=delta;
  alpha1_data.p=p;
  alpha1_data.xbeta1=xbeta1;
  alpha1_data.xbeta2=xbeta2;
  int err=arms_simple(ninit, &xl, &xr,compreg_logdalpha1,&alpha1_data,dometrop,&xprev, alpha);
  if(err!=0){
  *alpha=xprev;	
	}
}

struct compreg_alpha2_parm{
	double lambda2,alpha1,lambda1, alphaalpha, alphalambda, p;
	double *t;
	int *delta;
	double *xbeta1;
	double *xbeta2;
	int size;
};

double compreg_logdalpha2(double alpha, void* alpha2_data){
	struct compreg_alpha2_parm *d = (struct compreg_alpha2_parm*)alpha2_data;
	double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
		if((testreal(temp))&&(abs(temp)>=0.0)){
			 if(d->delta[i]==1){
				 temp=temp*1.0;
			 }else if(d->delta[i]==2){
				 temp+=logf2v2(d->t[i],d->p,alpha,d->lambda2, d->xbeta1[i],d->xbeta2[i]);
			 }else{
				 temp+=logScomp(d->t[i],d->alpha1,d->lambda1,alpha,
				 d->lambda2,d->xbeta1[i],d->xbeta2[i],d->p);
			 }
		}
	}
	return (d->alphaalpha-1.0)*log(alpha)-(d->alphalambda)*alpha+temp;
}

void compreg_samplealpha2(double* alpha, double lambda2,
double alpha1, double lambda1,
double p,
double alphaalpha, double alphalambda,
 double* t,
 int* delta, int size, double* xbeta1, double *xbeta2)
{
  double xl=findbase(lambda1);
  double xr=80.0;
  int ninit=5;
  int dometrop = 1;
  double xprev = *alpha;
  struct compreg_alpha2_parm alpha2_data;
  alpha2_data.alphaalpha=alphaalpha;
  alpha2_data.alphalambda=alphalambda;
  alpha2_data.alpha1=alpha1;
  alpha2_data.lambda1=lambda1;
  alpha2_data.lambda2=lambda2;
  alpha2_data.size=size;
  alpha2_data.t=t;
  alpha2_data.delta=delta;
  alpha2_data.p=p;
  alpha2_data.xbeta1=xbeta1;
  alpha2_data.xbeta2=xbeta2;
  int err=arms_simple(ninit, &xl, &xr,compreg_logdalpha2,&alpha2_data,dometrop,&xprev, alpha);
  if(err!=0){
  *alpha=xprev;	
	}
}

struct compreg_p_parm{
	double alpha1,lambda1,alpha2, lambda2, gamma0, gamma1;
	int size;
	double* t;
	int* delta;
	double *xbeta1;
	double *xbeta2;
};

double compreg_logdp(double p, void* p_data){
	struct compreg_p_parm *d = (struct compreg_p_parm*)p_data;
	double temp=0.0;
    for(int i=0;i<(d->size);i++){
		if(d->delta[i]==1){
			temp+=logf1v2(d->t[i],p,d->alpha1, d->lambda1, d->xbeta1[i]);
		}else if(d->delta[i]==2){
			temp+=logf2v2(d->t[i],p,d->alpha2,d->lambda2, d->xbeta1[i],d->xbeta2[i]);
		}
		else{
			temp+=logScomp(d->t[i],d->alpha1,d->lambda1,d->alpha2,d->lambda2,
			d->xbeta1[i],d->xbeta2[i],p);
		}
	}
	double logd=(d->gamma0-1.0)*log(p)+(d->gamma1-1.0)*log(1.0-p)+temp;
 return logd;
}

void compreg_samplep(double* p, double gamma0, double gamma1,double *t, int* delta,
 double alpha1, double lambda1,double alpha2, double lambda2, int size,
 double* xbeta1, double *xbeta2)
{
    double xl=0.0;
	double xr=1.0;
	int  ninit = 5;
	int dometrop = 1;
	double xprev = *p;
    double sample=*p;
	struct compreg_p_parm p_data;
	p_data.gamma0=gamma0;
	p_data.gamma1=gamma1;
	p_data.alpha1=alpha1;
	p_data.lambda1=lambda1;
	p_data.alpha2=alpha2;
	p_data.lambda2=lambda2;
    p_data.size = size;
	p_data.t = t;
    p_data.delta=delta;
	p_data.xbeta1=xbeta1;
	p_data.xbeta2=xbeta2;
	int err = arms_simple(ninit, &xl, &xr, compreg_logdp, &p_data, dometrop, &xprev, &sample);
	if(testreal(sample)|(err!=0)){
			*p=sample;
	}
}

struct compreg_beta1_parm{
	double betasl, alpha1, lambda1, alpha2, lambda2,p;
	int size;
    double *t;
	int *delta;
	double *x;
	double *xbeta1left;
	double *xbeta2;
};

double compreg_logdbeta1(double beta, void* beta1_data){
	struct compreg_beta1_parm *d = (struct compreg_beta1_parm*)beta1_data;
	double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
	if((testreal(temp))&&(abs(temp)>=0.0)){
		if(d->delta[i]==1){
			temp+=logf1v2(d->t[i],d->p,d->alpha1, d->lambda1, d->xbeta1left[i]+d->x[i]*beta);
		}else if(d->delta[i]==2){
			temp+=logf2v2(d->t[i],d->p,d->alpha2,d->lambda2, d->xbeta1left[i]+d->x[i]*beta,d->xbeta2[i]);
		}
		else{
			temp+=logScomp(d->t[i],d->alpha1,d->lambda1,d->alpha2,d->lambda2,
			d->xbeta1left[i]+d->x[i]*beta,d->xbeta2[i],d->p);
		}
	 }
	}
	return -log(1.0+pow(beta,2.0)/pow(d->betasl,2.0))+temp;
}

void compreg_samplebeta1(double* beta, double betasl,
 double alpha1,double lambda1,
 double alpha2, double lambda2,
 double* t,
 int* delta, int size, double*x, double *xbeta1left, double* xbeta2, double p) 
{
	double xl = -10.0;
	double xr = 10.0;
	int  ninit = 4;
	int dometrop = 1;
	double xprev = *beta;
	struct compreg_beta1_parm beta1_data;
	beta1_data.betasl=betasl;
	beta1_data.lambda1 = lambda1;
	beta1_data.alpha1 = alpha1;
	beta1_data.lambda2 = lambda2;
	beta1_data.alpha2 = alpha2;
	beta1_data.size = size;
	beta1_data.t=t;
	beta1_data.delta=delta;
	beta1_data.x = x;
	beta1_data.xbeta1left = xbeta1left;
	beta1_data.xbeta2 = xbeta2;
	beta1_data.p=p;
	int err= arms_simple(ninit, &xl, &xr, compreg_logdbeta1, &beta1_data, dometrop, &xprev, beta);
        if(err!=0){
	 *beta=xprev;	
	}
}


struct compreg_beta2_parm{
	double betasl, alpha1, lambda1, alpha2, lambda2,p;
	int size;
    double *t;
	int *delta;
	double *x;
	double *xbeta2left;
	double *xbeta1;
};

double compreg_logdbeta2(double beta, void* beta2_data){
	struct compreg_beta2_parm *d = (struct compreg_beta2_parm*)beta2_data;
	double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
	if((testreal(temp))&&(abs(temp)>=0.0)){
		if(d->delta[i]==1){
			temp=temp*1.0;
		}else if(d->delta[i]==2){
			temp+=logf2v2(d->t[i],d->p,d->alpha2,d->lambda2, 
			d->xbeta1[i],d->xbeta2left[i]+d->x[i]*beta);
		}
		else{
			temp+=logScomp(d->t[i],d->alpha1,d->lambda1,d->alpha2,d->lambda2,d->xbeta1[i],
			d->xbeta2left[i]+d->x[i]*beta,d->p);
		}
	 }

	}
	return -log(1.0+pow(beta,2.0)/pow(d->betasl,2.0))+temp;
}

void compreg_samplebeta2(double* beta, double betasl, double alpha1,double lambda1,
 double alpha2, double lambda2,
 double* t,
 int* delta, int size, double*x, double *xbeta2left, double* xbeta1, double p) 
{
	double xl = -10.0;
	double xr = 10.0;
	int  ninit = 4;
	int dometrop = 1;
	double xprev = *beta;
	struct compreg_beta2_parm beta2_data;
	beta2_data.betasl=betasl;
	beta2_data.lambda1 = lambda1;
	beta2_data.alpha1 = alpha1;
	beta2_data.lambda2 = lambda2;
	beta2_data.alpha2 = alpha2;
	beta2_data.size = size;
	beta2_data.t=t;
	beta2_data.delta=delta;
	beta2_data.x = x;
	beta2_data.xbeta2left = xbeta2left;
	beta2_data.xbeta1 = xbeta1;
	beta2_data.p=p;
	int err= arms_simple(ninit, &xl, &xr, compreg_logdbeta2, &beta2_data, dometrop, &xprev, beta);
        if(err!=0){
	 *beta=xprev;	
	}
	}
