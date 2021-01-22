// [[Rcpp::depends(iBATCGH]]
#include <pRMHmixture.h>

// [[Rcpp::interfaces(cpp)]]
//Add or Delete step
double pRAD (int choice,int chg,int m,double add0,double add1,double* mGam,double* mOm1,double* mOm2,int* mRnew,int* mR,int mg){
double out=0;
int choice1=0;
//rgm
	out=out+pr(choice,m,chg,add0,add1,mGam,mOm1,mOm2,mRnew,mR,mg);
//rg(m-1)
if(choice!=0){
	choice1=choice-1;
	out=out+pr(choice1,m,chg,add0,add1,mGam,mOm1,mOm2,mRnew,mR,mg);
}
//rg(m+1)
if(choice!=(m-2)){
	choice1=choice+1;
	out=out+pr(choice1,m,chg,add0,add1,mGam,mOm1,mOm2,mRnew,mR,mg);
}
return out;
}

// [[Rcpp::interfaces(cpp)]]
//Swap step
double pRS (int choice,int chg,int m,double add0,double add1,double* mGam,double* mOm1,double* mOm2,int* mRnew,int* mR,int* mFlag,arma::ivec *flagp,int mg){
arma::ivec flag(mFlag,mg,false);
double out=0;
int choice1=0;
//rgm
if(flag[chg+choice]!=1){
	flag[chg+choice]=1;
	out=out+pr(choice,m,chg,add0,add1,mGam,mOm1,mOm2,mRnew,mR,mg);
}
//rg(m-1)
if(choice!=0){
	choice1=choice-1;
	if(flag[chg+choice1]!=1){
		flag[chg+choice1]=1;
		out=out+pr(choice1,m,chg,add0,add1,mGam,mOm1,mOm2,mRnew,mR,mg);
	}
}
//rg(m+1)
if(choice!=(m-2)){
	choice1=choice+1;
	if(flag[chg+choice1]!=1){
		flag[chg+choice1]=1;
		out=out+pr(choice1,m,chg,add0,add1,mGam,mOm1,mOm2,mRnew,mR,mg);
	}
}
*flagp=flag;
return out;
}

