#ifndef __MHXImixture_h__
#define __MHXImixture_h__

#include "iBATCGH_RcppExports.h"

double pRXI (int choice,int chg,int m,double add0,double add1,double* mGam,double* mOm1,double* mOm2,double* mGamNew,double* mOm1New,double* mOm2New,int* mR, int mg);
double mml(int* mR,int* mX,double* mY, double c, double d, double delta, int elem,arma::mat Hs,int nOss,int p,int g,int mg);
double ldnwpi(double x, double mean, double sigma);

#define DEN (exp(1.0)-1.0)

#endif // __MHXImixture_h__

