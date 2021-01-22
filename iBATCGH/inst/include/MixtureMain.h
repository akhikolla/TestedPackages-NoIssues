#ifndef __MixtureMain_h__
#define __MixtureMain_h__

#include "iBATCGH_RcppExports.h"

double MHR(double phi,int* mR,double* mX,double* mY,double c,double d,double delta,double* mGam,double* mOm1,double* mOm2,double e,double f,arma::icolvec *Rout,double add0,double add1,int* mXI,arma::mat Hs,arma::vec countnotwo,int nR,int choiceg,int selectioncgh,int m,int mg, int g, int nOss,int indep);
double MHxi(int* mXI,arma::mat A,arma::mat Picum,double* mDist,double disfix,double alpha,double* mGam,double* mOm1,double* mOm2,double* mSum,double e,double f,int* mR,double* mY,double c,double d,double delta,arma::colvec mu,arma::colvec sigma,double* mX,arma::mat tran,arma::colvec Pi,arma::imat *xiout,arma::colvec *gammaout,arma::colvec *omega1out,arma::colvec *omega2out, arma::colvec *sumout,double den,double add0,double add1,arma::mat Hs,int choice,int sel,int m,int mg,int g,int nOss, int indep);
void omegagamma(arma::mat distance,arma::imat xi,double disfix,double alpha,arma::colvec *gammaex,arma::colvec *omega1ex,arma::colvec *omega2ex,arma::colvec *sumex);
arma::colvec initial(arma::mat A) ;
arma::rowvec dir(arma::rowvec a);
double truncnorm(double mean,double sd,double lower,double upper);
double truncgamma(double a,double b,double lower);

#endif // __MixtureMain_h__

