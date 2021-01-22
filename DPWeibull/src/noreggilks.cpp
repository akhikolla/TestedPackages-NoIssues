#include <iostream>
#include <vector>
#include <Rmath.h>
#include <R.h>
#include <algorithm>
#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;
#include "arms.h"
#include "auxfuns.h"
#include "commonfunc.h"
struct noreg_loglambda_parm{
	double alpha, alpha0, lambda0;
	double *tl;
	double *tr;
	int *delta;
	int *pi;
	int size;
};

struct noreg_alpha_parm{
	double lambda, alphaalpha, alphalambda;
	double *tl;
	double *tr;
	int *delta;
	int *pi;
	int size;
};


double noreg_logdloglambda(double loglambda, void* loglambda_data){
	struct noreg_loglambda_parm *d = (struct noreg_loglambda_parm*)loglambda_data;
	double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
		if((testreal(temp))&&(abs(temp)>=0.0)){
			 if((d->delta[i]==1)&&(d->pi[i]==0)){
				 temp+=logdWeibloglambda(d->tl[i],d->alpha,loglambda);
			 }else if((d->delta[i]==0)&&(d->pi[i]==1)){
				 temp+=logsWeibloglambda(d->tl[i],d->alpha,loglambda);
			 }else{
				 double cushion=sWeibloglambda(d->tl[i],d->alpha,loglambda)-
				 sWeibloglambda(d->tr[i],d->alpha,loglambda);
				 if((cushion>0.0)&&(testreal(cushion))){
					 temp+=log(cushion);
				 }
			 }
		}
	}
	return (d->alpha0)*loglambda-(d->lambda0)*exp(loglambda)+temp;
}

void noreg_sampleloglambda(double* lambda, double alpha,
double alpha0, double lambda0,
 double* tl, double* tr,
 int* delta, int* pi, int size)
{ 
  double loglambda=log(*lambda);
  double xl=log(inversebase(alpha));
  double xr=10.0;
  int ninit=5;
  int dometrop = 1;
  double xprev = loglambda;
  struct noreg_loglambda_parm loglambda_data;
  loglambda_data.alpha0=alpha0;
  loglambda_data.lambda0=lambda0;
  loglambda_data.alpha=alpha;
  loglambda_data.size=size;
  loglambda_data.tl=tl;
  loglambda_data.tr=tr;
  loglambda_data.delta=delta;
  loglambda_data.pi=pi;
  int err=arms_simple(ninit, &xl, &xr,noreg_logdloglambda,(void*)&loglambda_data,dometrop,&xprev, &loglambda);
  if(testreal(loglambda)|(err!=0)){
	  *lambda=exp(loglambda);
  }
}


double noreg_logdalpha(double alpha, void* alpha_data){
	struct noreg_alpha_parm *d = (struct noreg_alpha_parm*)alpha_data;
	double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
		if((testreal(temp))&&(abs(temp)>=0.0)){
			if((d->delta[i]==1)&&(d->pi[i]==0)){
			 temp+=logdWeib(d->tl[i],alpha,d->lambda);
			}else if((d->delta[i]==0)&&(d->pi[i]==1)){
			 temp+=logsWeib(d->tl[i],alpha,d->lambda);
			}else{
			 double cushion=sWeib(d->tl[i],alpha,d->lambda)-sWeib(d->tr[i],alpha,d->lambda);
			 if((cushion>0.0)&&(testreal(cushion))){
				temp+=log(cushion);
			 }
			}
		 }

	}
	return (d->alphaalpha-1.0)*log(alpha)-(d->alphalambda)*alpha+temp;
}

void noreg_samplealpha(double* alpha, double lambda,
double alphaalpha, double alphalambda,
 double* tl, double* tr,
 int* delta, int* pi, int size)
{
  double xl=findbase(lambda);
  double xr=80.0;
  int ninit=5;
  int dometrop = 1;
  double xprev = *alpha;
  struct noreg_alpha_parm alpha_data;
  alpha_data.alphaalpha=alphaalpha;
  alpha_data.alphalambda=alphalambda;
  alpha_data.lambda=lambda;
  alpha_data.size=size;
  alpha_data.tl=tl;
  alpha_data.tr=tr;
  alpha_data.delta=delta;
  alpha_data.pi=pi;
  int err=arms_simple(ninit, &xl, &xr,noreg_logdalpha,(void*)&alpha_data,dometrop,&xprev, alpha);
if(err!=0){
*alpha=xprev;
}
}
