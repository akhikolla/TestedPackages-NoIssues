#ifndef __MHRProbit_h__
#define __MHRProbit_h__

#include "iBATCGH_RcppExports.h"

double pRAD2 (int choice,int chg,int m,double alpha0,double alpha1,double* mSum,int* mRnew,int* mR,int mg);
double pRS2 (int choice,int chg,int m,double alpha0,double alpha1,double* mSum,int* mRnew,int* mR,int* mFlag,arma::ivec *flagp,int mg);
double mml(int* mR,int* mX,double* mY, double c, double d, double delta, int elem,arma::mat Hs,int nOss,int p,int g,int mg);

#endif // __MHRProbit_h__

