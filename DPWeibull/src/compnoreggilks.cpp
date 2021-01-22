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


struct compnoreg_loglambda_parm{
	double alpha1,  lambda01, alpha2,  lambda2,alpha0, p;
	double *t;
	int *delta;
	int size;
};

double compnoreg_logdloglambda(double loglambda, void* loglambdacomp_data){
struct compnoreg_loglambda_parm *d = (struct compnoreg_loglambda_parm*)loglambdacomp_data;
double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
		if((testreal(temp))&&(abs(temp)>=0.0)){
			 if(d->delta[i]==1){
				 temp+=logdWeibloglambda(d->t[i],d->alpha1,loglambda);
			 }else if(d->delta[i]==2){
				 temp=temp*1.0;
			 }else{
				 if((d->p)*pWeib(d->t[i],d->alpha1,exp(loglambda))
					 +(1.0-d->p)*pWeib(d->t[i],d->alpha2,d->lambda2)<1.0){
						 temp+=log(1.0-(d->p)*pWeib(d->t[i],d->alpha1,exp(loglambda))
					 -(1.0-d->p)*pWeib(d->t[i],d->alpha2,d->lambda2));
				}else{
					temp+=-pow(100.0,100.0);
				}
			 }
		}
	}
	return (d->alpha0)*loglambda-(d->lambda01)*exp(loglambda)+temp;
}

void compnoreg_sampleloglambda(double* lambda,
double alpha0, double lambda01,
 double alpha1, 
 double alpha2,double lambda2,  double p,
 double* t,
 int* delta, int size)
{ 
  double loglambda=log(*lambda);
  double xl=log(inversebase(alpha1));
  double xr=1.0;
  int ninit=5;
  int dometrop = 1;
  double xprev = loglambda;
  struct compnoreg_loglambda_parm loglambdacomp_data;
  loglambdacomp_data.alpha0=alpha0;
  loglambdacomp_data.lambda01=lambda01;
  loglambdacomp_data.alpha1=alpha1;
  loglambdacomp_data.alpha2=alpha2;
  loglambdacomp_data.lambda2=lambda2;
  loglambdacomp_data.p=p;
  loglambdacomp_data.size=size;
  loglambdacomp_data.t=t;
  loglambdacomp_data.delta=delta;
  int err=arms_simple(ninit, &xl, &xr,compnoreg_logdloglambda,&loglambdacomp_data,dometrop,&xprev, &loglambda);
  if(testreal(loglambda)|(err!=0)){
	  *lambda=exp(loglambda);
  }
}

struct compnoreg_alpha_parm{
	double lambda1,alpha2,lambda2, alphaalpha, alphalambda, p;
	double *t;
	int *delta;
	int size;
};

double compnoreg_logdalpha(double alpha, void* alpha_data){
	struct compnoreg_alpha_parm *d = (struct compnoreg_alpha_parm*)alpha_data;
	double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
		if((testreal(temp))&&(abs(temp)>=0.0)){
			 if(d->delta[i]==1){
				 temp+=logdWeib(d->t[i],alpha,d->lambda1);
			 }else if(d->delta[i]==2){
				 temp=temp*1.0;
			 }else{
				 if((d->p)*pWeib(d->t[i],alpha,d->lambda1)
					 +(1.0-d->p)*pWeib(d->t[i],d->alpha2,d->lambda2)<1.0){
						 temp+=log(1.0-(d->p)*pWeib(d->t[i],alpha,d->lambda1)
					 -(1.0-d->p)*pWeib(d->t[i],d->alpha2,d->lambda2));
				}else{
					temp+=-pow(100.0,100.0);
				}
			 }
		}
	}
	return (d->alphaalpha-1.0)*log(alpha)-(d->alphalambda)*alpha+temp;
}

void compnoreg_samplealpha(double* alpha, double lambda1,
double alpha2, double lambda2,
double p,
double alphaalpha, double alphalambda,
 double* t,
 int* delta, int size)
{
  double xl=findbase(lambda1);
  double xr=80.0;
  int ninit=5;
  int dometrop = 1;
  double xprev = *alpha;
  struct compnoreg_alpha_parm alpha_data;
  alpha_data.alphaalpha=alphaalpha;
  alpha_data.alphalambda=alphalambda;
  alpha_data.alpha2=alpha2;
  alpha_data.lambda1=lambda1;
  alpha_data.lambda2=lambda2;
  alpha_data.size=size;
  alpha_data.t=t;
  alpha_data.delta=delta;
  alpha_data.p=p;
  int err=arms_simple(ninit, &xl, &xr,compnoreg_logdalpha,&alpha_data,dometrop,&xprev, alpha);
  if(err!=0){
	  *alpha=xprev;
  }
}

struct compnoreg_p_parm{
	double alpha1,lambda1,alpha2, lambda2, gamma0, gamma1;
	int size;
	double* t;
	int* delta;
};

double compnoreg_logdp(double p, void* p_data){
	struct compnoreg_p_parm *d = (struct compnoreg_p_parm*)p_data;
	double temp=0.0;
    for(int i=0;i<(d->size);i++){
		if(d->delta[i]==1){
			temp+=log(p);
		}else if(d->delta[i]==2){
			temp+=log(1.0-p);
		}
		else{
			 if(p*pWeib(d->t[i],d->alpha1,d->lambda1)
			+(1.0-p)*pWeib(d->t[i],d->alpha2,d->lambda2)<1.0){
			 temp+=log(1.0-p*pWeib(d->t[i],d->alpha1,d->lambda1)
			-(1.0-p)*pWeib(d->t[i],d->alpha2,d->lambda2));
			}else{
					temp+=-pow(100.0,100.0);
				}
		}
	}
	double logd=(d->gamma0-1.0)*log(p)+(d->gamma1-1.0)*log(1.0-p)+temp;
 return logd;
}

void compnoreg_samplep(double* p, double gamma0, double gamma1,double *t, int* delta,
 double alpha1, double lambda1,double alpha2, double lambda2, int size)
{
    double xl=0.0;
	double xr=1.0;
	int  ninit = 5;
	int dometrop = 1;
	double xprev = 0.5;
    double sample=*p;
	struct compnoreg_p_parm p_data;
	p_data.gamma0=gamma0;
	p_data.gamma1=gamma1;
	p_data.alpha1=alpha1;
	p_data.lambda1=lambda1;
	p_data.alpha2=alpha2;
	p_data.lambda2=lambda2;
    p_data.size = size;
	p_data.t = t;
    p_data.delta=delta;
	int err = arms_simple(ninit, &xl, &xr, compnoreg_logdp, &p_data, dometrop, &xprev, &sample);
	if(testreal(sample)|(err!=0)){
			*p=sample;
	}
}
