#ifndef __MHRmixture_h__
#define __MHRmixture_h__

#include "iBATCGH_RcppExports.h"

double pRAD (int choice,int chg,int m,double add0,double add1,double* mGam,double* mOm1,double* mOm2,int* mRnew,int* mR,int mg);
double pRS (int choice,int chg,int m,double add0,double add1,double* mGam,double* mOm1,double* mOm2,int* mRnew,int* mR,int* mFlag,arma::ivec *flagp,int mg);
double mml(int* mR,int* mX,double* mY, double c, double d, double delta, int elem,arma::mat Hs,int nOss,int p,int g,int mg);

#endif // __MHRmixture_h__

