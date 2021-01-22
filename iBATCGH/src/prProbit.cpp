#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(cpp)]]
double prProbit(int choice,int m,int chg,double alpha0,double alpha1,double* mSum,int* mR,int mg){
arma::icolvec R(mR,mg,false);
arma::colvec sum(mSum,(m+1),false);
double out=0;
double q=0;

	if(choice==0){
		q=alpha0+alpha1*sum[choice+1]*((R[chg+choice+1]==1)?(-1):1);
	}
	if(choice==(m-1)){
		q=alpha0+alpha1*sum[choice]*((R[chg+choice-1]==1)?(-1):1);
	}
	if(choice>0 && choice<(m-1)){
		q=alpha0+alpha1*sum[choice]*((R[chg+choice-1]==1)?(-1):1)+alpha1*sum[choice+1]*((R[chg+choice+1]==1)?(-1):1);
	}
if(R[chg+choice]==0){out=Rf_pnorm5(q,0,1,1,1);}
if(R[chg+choice]==1){out=Rf_pnorm5(q,0,1,0,1);}
return out;
}

