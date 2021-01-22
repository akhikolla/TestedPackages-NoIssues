// [[Rcpp::depends(iBATCGH]]
#include <pRMHXImixture.h>

// [[Rcpp::interfaces(cpp)]]
//Probability for MH xi
double pRXI (int choice,int chg,int m,double add0,double add1,double* mGam,double* mOm1,double* mOm2,double* mGamNew,double* mOm1New,double* mOm2New,int* mR, int mg){
arma::icolvec R(mR,mg,false);
double add=((R[chg+choice]==0)?add0:add1);
double out=0;
int choice1=0;
//rgm
	out=out+pr2(choice,m,chg,add,mGam,mOm1,mOm2,mGamNew,mOm1New,mOm2New,mR,mg);
//rg(m-1)
if(choice!=0){
	choice1=choice-1;
	out=out+pr2(choice1,m,chg,add,mGam,mOm1,mOm2,mGamNew,mOm1New,mOm2New,mR,mg);
}
//rg(m+1)
if(choice!=(m-2)){
	choice1=choice+1;
	out=out+pr2(choice1,m,chg,add,mGam,mOm1,mOm2,mGamNew,mOm1New,mOm2New,mR,mg);
}
return out;
}

