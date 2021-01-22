#ifndef __MHXIProbit_h__
#define __MHXIProbit_h__

#include "iBATCGH_RcppExports.h"

double pRXI2 (int choice,int chg,int m,double alpha0,double alpha1,double* mSum,double* mSumNew,int* mR,int mg);
double mml(int* mR,int* mX,double* mY, double c, double d, double delta, int elem,arma::mat Hs,int nOss,int p,int g,int mg);
double ldnwpi(double x, double mean, double sigma);

#define DEN (exp(1.0)-1.0)

#endif // __MHXIProbit_h__

