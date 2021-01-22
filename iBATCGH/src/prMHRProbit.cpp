// [[Rcpp::depends(iBATCGH]]
#include <pRMHProbit.h>

// [[Rcpp::interfaces(cpp)]]
//Add or Delete step
double pRAD2 (int choice,int chg,int m,double alpha0,double alpha1,double* mSum,int* mRnew,int* mR,int mg){
double out=0;
int choice1=0;
//rgm
	out=out+prProbit(choice,m,chg,alpha0,alpha1,mSum,mRnew,mg)-prProbit(choice,m,chg,alpha0,alpha1,mSum,mR,mg);
//rg(m-1)
if(choice!=0){
	choice1=choice-1;
	out=out+prProbit(choice1,m,chg,alpha0,alpha1,mSum,mRnew,mg)-prProbit(choice1,m,chg,alpha0,alpha1,mSum,mR,mg);
}
//rg(m+1)
if(choice!=(m-1)){
	choice1=choice+1;
	out=out+prProbit(choice1,m,chg,alpha0,alpha1,mSum,mRnew,mg)-prProbit(choice1,m,chg,alpha0,alpha1,mSum,mR,mg);
}
return out;
}

// [[Rcpp::interfaces(cpp)]]
//Swap step
double pRS2 (int choice,int chg,int m,double alpha0,double alpha1,double* mSum,int* mRnew,int* mR,int* mFlag,arma::ivec *flagp,int mg){
arma::ivec flag(mFlag,mg,false);
double out=0;
int choice1=0;
//rgm
if(flag[chg+choice]!=1){
	flag[chg+choice]=1;
	out=out+prProbit(choice,m,chg,alpha0,alpha1,mSum,mRnew,mg)-prProbit(choice,m,chg,alpha0,alpha1,mSum,mR,mg);
}
//rg(m-1)
if(choice!=0){
	choice1=choice-1;
	if(flag[chg+choice1]!=1){
		flag[chg+choice1]=1;
		out=out+prProbit(choice1,m,chg,alpha0,alpha1,mSum,mRnew,mg)-prProbit(choice1,m,chg,alpha0,alpha1,mSum,mR,mg);
	}
}
//rg(m+1)
if(choice!=(m-1)){
	choice1=choice+1;
	if(flag[chg+choice1]!=1){
		flag[chg+choice1]=1;
		out=out+prProbit(choice1,m,chg,alpha0,alpha1,mSum,mRnew,mg)-prProbit(choice1,m,chg,alpha0,alpha1,mSum,mR,mg);
	}
}
*flagp=flag;
return out;
}

