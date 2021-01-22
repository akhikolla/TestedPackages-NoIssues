// [[Rcpp::depends(iBATCGH]]
#include <prMixture.h>

// [[Rcpp::interfaces(cpp)]]
double pr2(int choice,int m,int chg,double add,double* mGam,double* mOm1,double* mOm2,double* mGamNew,double* mOm1New,double* mOm2New,int* mR, int mg){
arma::icolvec R(mR,mg,false);
arma::colvec gamma(mGam,m,false);
arma::colvec omega1(mOm1,m,false);
arma::colvec omega2(mOm2,m,false);
arma::colvec gammanew(mGamNew,m,false);
arma::colvec omega1new(mOm1New,m,false);
arma::colvec omega2new(mOm2New,m,false);

double out=0;
	if(choice==0){
		out=out+pRedge(gammanew[choice],omega2new[choice],R[chg+choice],R[chg+choice+1],add);
		out=out-pRedge(gamma[choice],omega2[choice],R[chg+choice],R[chg+choice+1],add);
	}
	if(choice==(m-1)){
		out=out+pRedge(gammanew[choice],omega1new[choice],R[chg+choice],R[chg+choice-1],add);
		out=out-pRedge(gamma[choice],omega1[choice],R[chg+choice],R[chg+choice-1],add);
	}
	if(choice>0 && choice<(m-1)){
		out=out+pRmiddle(gammanew[choice],omega1new[choice],omega2new[choice],R[chg+choice],R[chg+choice-1],R[chg+choice+1],add);
		out=out-pRmiddle(gamma[choice],omega1[choice],omega2[choice],R[chg+choice],R[chg+choice-1],R[chg+choice+1],add);
	}
return out;
}

