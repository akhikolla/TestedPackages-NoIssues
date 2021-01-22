// [[Rcpp::depends(iBATCGH]]
#include <pRMHXIProbit.h>

// [[Rcpp::interfaces(cpp)]]
//Probability for MH xi
double pRXI2 (int choice,int chg,int m,double alpha0,double alpha1,double* mSum,double* mSumNew,int* mR,int mg){
double out=0;
int choice1=0;
//rgm
	out=out+prProbit(choice,m,chg,alpha0,alpha1,mSumNew,mR,mg)-prProbit(choice,m,chg,alpha0,alpha1,mSum,mR,mg);
//rg(m-1)
if(choice!=0){
	choice1=choice-1;
	out=out+prProbit(choice1,m,chg,alpha0,alpha1,mSumNew,mR,mg)-prProbit(choice1,m,chg,alpha0,alpha1,mSum,mR,mg);
}
//rg(m+1)
if(choice!=(m-1)){
	choice1=choice+1;
	out=out+prProbit(choice1,m,chg,alpha0,alpha1,mSumNew,mR,mg)-prProbit(choice1,m,chg,alpha0,alpha1,mSum,mR,mg);
}
return out;
}

