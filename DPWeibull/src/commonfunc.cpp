#include <vector>
#include <Rmath.h>
#include <R.h>
#include <algorithm>
#include <math.h>
#include <Rcpp.h>

#include "arms.c"

using namespace Rcpp ;

struct truncauchy_parm{
	double sl;
};
//find the percentile of each column in a matrix
NumericVector colpercentileRcpp(NumericMatrix x, double percentile) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    int position = nrow*percentile; // Euclidian division
    NumericVector out(ncol);
    for (int j = 0; j < ncol; j++) { 
        NumericVector y = x(_,j); // Copy column -- original will not be mod
        std::nth_element(y.begin(), y.begin() + position, y.end()); 
        out[j] = y[position];  
    }
    return out;
}

NumericVector matrixtimesvector(NumericMatrix x, NumericVector beta){
  NumericVector result(x.nrow());
	for(int i=0;i<x.nrow();i++){
		result(i)=sum(x(i,_)*beta);
	}
	return result;
}

//make a submatrix by making the unwanted column be 0 
NumericMatrix submat_rcpp(NumericMatrix X, LogicalVector condition) { 
    int n=X.nrow(), k=X.ncol();
    NumericMatrix out(sum(condition),k);
    for (int i = 0, j = 0; i < n; i++) {
        if(condition[i]) {
            out(j,_) = X(i,_);
            j = j+1;
        }
    }
    return(out);
}
//sample from a multinomial distribution
IntegerVector oneMultinomCalt(NumericVector probs) {
    int k = probs.size();
    IntegerVector ans(k);
    rmultinom(1, probs.begin(), k, ans.begin());
    return(ans);
}

//concentration parameter generation, use the method from Escobar and West
double nugen(const double nu, const int npts, const int ngrp,const double a, const double b)
{
	double eta=rbeta(1,nu + 1.0, double(npts))[0];
	double pi_eta = (a + double(ngrp) - 1.0) / (double(npts)*(b - log(eta)) + (a + double(ngrp) - 1.0));
	if(rbinom(1,1.0, pi_eta)[0] == 1.0){
		return R::rgamma(a + double(ngrp),1.0/(b - log(eta)));
	}else{
		return R::rgamma(a + double(ngrp)-1.0,1.0/(b - log(eta)));
	}
}
// for fixed lambda, find the lower boundary for alpha
double findbase(double lambda){
	if(lambda>0){
		return std::min(std::max(0.0,log(-log(0.05)/lambda)/log(25.0)),80.0);
	}else{
		return 80.0;
	}
}
// for fixed alpha, find the lower boundary for log lambda
double inversebase(double alpha){
	return -log(0.05)/pow(25.0,alpha);
}
//if the result is a real value
bool testreal(double d) { return d == d; }


double truncauchy(double truncauchy, void* truncauchy_data){
	struct truncauchy_parm *d = (struct truncauchy_parm*)truncauchy_data;
	return -log(1.0+pow(truncauchy,2.0)/pow(d->sl,2.0));
	}
void samptruncauchy(double *xsamp,double sl){
	double xl = -10.0;
	double xr = 10.0;
	int  ninit = 4;
	int dometrop = 1;
	double xprev = *xsamp;
	struct truncauchy_parm truncauchy_data;
	truncauchy_data.sl=sl;
   int err = arms_simple(ninit, &xl, &xr, truncauchy, (void*)&truncauchy_data, dometrop, &xprev, xsamp);
   if(err!=0){
		*xsamp=xprev;
	}
}

