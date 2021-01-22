// [[Rcpp::depends(iBATCGH]]
#include <prMixture.h>

// [[Rcpp::interfaces(cpp)]]
double pr(int choice,int m,int chg,double add0,double add1,double* mGam,double* mOm1,double* mOm2,int* mRnew,int* mR,int mg){
arma::icolvec R(mR,mg,false);
arma::icolvec Rnew(mRnew,mg,false);
arma::colvec gamma(mGam,m,false);
arma::colvec omega1(mOm1,m,false);
arma::colvec omega2(mOm2,m,false);

double out=0;
double add=((R[chg+choice]==0)?add0:add1);
double addnew=((Rnew[chg+choice]==0)?add0:add1);
	if(choice==0){
		out=out+pRedge(gamma[choice],omega2[choice],Rnew[chg+choice],Rnew[chg+choice+1],addnew);
		out=out-pRedge(gamma[choice],omega2[choice],R[chg+choice],R[chg+choice+1],add);
	}
	if(choice==(m-1)){
		out=out+pRedge(gamma[choice],omega1[choice],Rnew[chg+choice],Rnew[chg+choice-1],addnew);
		out=out-pRedge(gamma[choice],omega1[choice],R[chg+choice],R[chg+choice-1],add);
	}
	if(choice>0 && choice<(m-1)){
		out=out+pRmiddle(gamma[choice],omega1[choice],omega2[choice],Rnew[chg+choice],Rnew[chg+choice-1],Rnew[chg+choice+1],addnew);
		out=out-pRmiddle(gamma[choice],omega1[choice],omega2[choice],R[chg+choice],R[chg+choice-1],R[chg+choice+1],add);
	}
return out;
}

