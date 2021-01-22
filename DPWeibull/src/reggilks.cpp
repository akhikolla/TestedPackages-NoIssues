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
#include "noreg.h"

struct reg_loglambda_parm{
	double alpha, alpha0, lambda0;
	double *tl;
	double *tr;
	int *delta;
	int *pi;
	int size;
	double *xbeta;
};

struct reg_alpha_parm{
	double loglambda, alphaalpha, alphalambda;
	double *tl;
	double *tr;
	int *delta;
	int *pi;
	int size;
	double* xbeta;
};

struct reg_beta_parm{
	double betasl, alpha;
	int size;
	double *tl;
	double *tr;
	int *delta;
	int *pi;
	double* x;
	double *loglambda;
};

struct reg_beta2_parm{
	double betasl, alpha;
	int size;
	double *tl;
	double *tr;
	int *delta;
	int *pi;
	double *loglambda;
};

double reg_logdloglambda(double loglambda, void* loglambdareg_data){
	struct reg_loglambda_parm *d = (struct reg_loglambda_parm*)loglambdareg_data;
	double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
		if((testreal(temp))&&(abs(temp)>=0.0)){
			 if((d->delta[i]==1)&&(d->pi[i]==0)){
				 temp+=logdWeibloglambda(d->tl[i],d->alpha,loglambda+d->xbeta[i]);
			 }else if((d->delta[i]==0)&&(d->pi[i]==1)){
				 temp+=logsWeibloglambda(d->tl[i],d->alpha,loglambda+d->xbeta[i]);
			 }else{
				 double cushion=sWeibloglambda(d->tl[i],d->alpha,loglambda+d->xbeta[i])-
				 sWeibloglambda(d->tr[i],d->alpha,loglambda+d->xbeta[i]);
				 if((cushion>0.0)&&(testreal(cushion))){
					 temp+=log(cushion);
				 }else{
				temp+=-pow(100.0,100.0);
				}
			 }
		}
	}
	return (d->alpha0)*loglambda-(d->lambda0)*exp(loglambda)+temp;
}

void reg_sampleloglambda(double* lambda, double alpha,
double alpha0, double lambda0,
 double* tl, double* tr,
 int* delta, int* pi, int size, double* xbeta)
{ 
  double loglambda=log(*lambda);
  double xl=log(inversebase(alpha));
  double xr=10.0;
  int ninit=5;
  int dometrop = 1;
  double xprev = loglambda;
  struct reg_loglambda_parm loglambdareg_data;
  loglambdareg_data.alpha0=alpha0;
  loglambdareg_data.lambda0=lambda0;
  loglambdareg_data.alpha=alpha;
  loglambdareg_data.size=size;
  loglambdareg_data.tl=tl;
  loglambdareg_data.tr=tr;
  loglambdareg_data.delta=delta;
  loglambdareg_data.pi=pi;
  loglambdareg_data.xbeta=xbeta;  
  int err=arms_simple(ninit, &xl, &xr,reg_logdloglambda,(void*)&loglambdareg_data,dometrop,&xprev, &loglambda);
    if(testreal(loglambda)|(err!=0)){
	  *lambda=exp(loglambda);
  }
}

double reg_logdalpha(double alpha, void* alpha_data){
	struct reg_alpha_parm *d = (struct reg_alpha_parm*)alpha_data;
	double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
		if((testreal(temp))&&(abs(temp)>=0.0)){
			if((d->delta[i]==1)&&(d->pi[i]==0)){
			 temp+=logdWeibloglambda(d->tl[i],alpha,d->loglambda+d->xbeta[i]);
			}else if((d->delta[i]==0)&&(d->pi[i]==1)){
			 temp+=logsWeibloglambda(d->tl[i],alpha,d->loglambda+d->xbeta[i]);
			}else{
			 double cushion=sWeibloglambda(d->tl[i],alpha,d->loglambda+d->xbeta[i])-
			 sWeibloglambda(d->tr[i],alpha,d->loglambda+d->xbeta[i]);
			 if((cushion>0.0)&&(testreal(cushion))){
				temp+=log(cushion);
			 }else{
				temp+=-pow(100.0,100.0);
			}
		  }
		 }
	}
	return (d->alphaalpha-1.0)*log(alpha)-(d->alphalambda)*alpha+temp;
}

void reg_samplealpha(double* alpha, double loglambda,
double alphaalpha, double alphalambda,
 double* tl, double* tr,
 int* delta, int* pi, int size, double* xbeta)
{
  double xl=findbase(exp(loglambda));
  double xr=80.0;
  int ninit=5;
  int dometrop = 1;
  double xprev = *alpha;
  struct reg_alpha_parm alpha_data;
  alpha_data.alphaalpha=alphaalpha;
  alpha_data.alphalambda=alphalambda;
  alpha_data.loglambda=loglambda;
  alpha_data.size=size;
  alpha_data.tl=tl;
  alpha_data.tr=tr;
  alpha_data.delta=delta;
  alpha_data.pi=pi;
    alpha_data.xbeta=xbeta;
  int err=arms_simple(ninit, &xl, &xr,reg_logdalpha,(void*)&alpha_data,dometrop,&xprev, alpha);
  if(err!=0){
  *alpha=xprev;
  }
}

double reg_logdbeta(double beta, void* beta_data){
	struct reg_beta_parm *d = (struct reg_beta_parm*)beta_data;
	double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
	if((testreal(temp))&&(abs(temp)>=0.0)){
			if((d->delta[i]==1)&&(d->pi[i]==0)){
			 temp+=logdWeibloglambda(d->tl[i],d->alpha,d->loglambda[i]+d->x[i]*beta);
			}else if((d->delta[i]==0)&&(d->pi[i]==1)){
			 temp+=logsWeibloglambda(d->tl[i],d->alpha,d->loglambda[i]+d->x[i]*beta);
			}else{
			 double cushion=sWeibloglambda(d->tl[i],d->alpha,d->loglambda[i]+d->x[i]*beta)
			 -sWeibloglambda(d->tr[i],d->alpha,d->loglambda[i]+d->x[i]*beta);
			 if((cushion>0.0)&&(testreal(cushion))){
				temp+=log(cushion);
			 }else{
				temp+=-pow(100.0,100.0);
			}
			}
		 }
	}
	return -log(1.0+pow(beta,2.0)/pow(d->betasl,2.0))+temp;
}

void reg_samplebeta(double* beta, double betasl, double alpha,
 double* tl, double* tr,
 int* delta, int* pi, int size, double*x, double *loglambda) 
{
	double xl = -10.0;
	double xr = 10.0;
	int  ninit = 4;
	int dometrop = 1;
	double xprev = *beta;
	struct reg_beta_parm beta_data;
	beta_data.betasl=betasl;
	beta_data.loglambda = loglambda;
	beta_data.alpha = alpha;
	beta_data.size = size;
	beta_data.tl=tl;
	beta_data.tr=tr;
	beta_data.delta=delta;
	beta_data.pi=pi;
	beta_data.x = x;
	int err= arms_simple(ninit, &xl, &xr,reg_logdbeta, (void*)&beta_data, dometrop, &xprev, beta);
   if(err!=0){
		*beta=xprev;
	}
}


double reg_logdbeta2(double beta, void* beta_data){
	struct reg_beta2_parm *d = (struct reg_beta2_parm*)beta_data;
	double temp = 0.0;
	for (int i = 0; i<(d->size); i++){
	if((testreal(temp))&&(abs(temp)>=0.0)){
			if((d->delta[i]==1)&&(d->pi[i]==0)){
			 temp+=logdWeibloglambda(d->tl[i],d->alpha,d->loglambda[i]+beta);
			}else if((d->delta[i]==0)&&(d->pi[i]==1)){
			 temp+=logsWeibloglambda(d->tl[i],d->alpha,d->loglambda[i]+beta);
			}else{
			 double cushion=sWeibloglambda(d->tl[i],d->alpha,d->loglambda[i]+beta)
			 -sWeibloglambda(d->tr[i],d->alpha,d->loglambda[i]+beta);
			 if((cushion>0.0)&&(testreal(cushion))){
				temp+=log(cushion);
			 }
			}
		 }
	}

	return -pow(beta,2.0)/pow(d->betasl,2.0)/2.0+temp;
}

void reg_samplebeta2(double* beta, double betasl, double alpha,
 double* tl, double* tr,
 int* delta, int* pi, int size, double *loglambda) 
{
	double xl = -10.0;
	double xr = 10.0;
	int  ninit = 4;
	int dometrop = 1;
	double xprev = *beta;
	struct reg_beta2_parm beta_data;
	beta_data.betasl=betasl;
	beta_data.loglambda = loglambda;
	beta_data.alpha = alpha;
	beta_data.size = size;
	beta_data.tl=tl;
	beta_data.tr=tr;
	beta_data.delta=delta;
	beta_data.pi=pi;
	int err= arms_simple(ninit, &xl, &xr,reg_logdbeta2, (void*)&beta_data, dometrop, &xprev, beta);
  if(err!=0){
		*beta=xprev;
	}
}
