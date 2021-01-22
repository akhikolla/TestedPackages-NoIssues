//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::List splitPrc_cpp(arma::vec X){

	// new C++ splitPrc function
	arma::vec S = arma::sort(X);
	int s = S.size()/2;

	double sumDst = arma::sum(arma::abs(S-S[s]));
	double sumSqr = arma::sum(arma::pow(S-S[s],2));

	return Rcpp::List::create(
		Named("sVal")=S[s],
		Named("kPrc")=sumSqr/std::pow(sumDst,2));
}

// [[Rcpp::export]]
Rcpp::LogicalVector bound2R_cpp(Rcpp::NumericMatrix X, Rcpp::NumericVector Rk){

	int n = X.nrow();
	int m = X.ncol();

	Rcpp::LogicalVector Xb(n);

	for (int i=0; i<n; ++i){
		Xb[i] = true;
		for (int j=0; j<m; ++j){
			Xb[i] &= (X(i,j)>=Rk[j] && X(i,j)<=Rk[m+j]);
		}
	}

	return Xb;
}

// [[Rcpp::export]]
Rcpp::NumericVector getClusters_cpp(arma::mat X, arma::mat U, arma::mat W, arma::mat R){

	int n = X.n_rows;
	int m = X.n_cols;
	int k = std::pow(2, m);

	Rcpp::NumericVector A = Rcpp::wrap(arma::ones(n)*k);

	for(int i=0; i<n; ++i){
		if (all(U.row(i))) {
			arma::uvec w = arma::find(W.row(i)==arma::max(W.row(i)));
			if (w.n_elem==1){
				A[i] = w[0];
			}
			else {
				arma::vec u = arma::randu<vec>(w.n_elem);
				for (int k=0; k<w.n_elem; ++k){
					bool Xb = true;
					for (int j=0; j<m; ++j){
						Xb &= (X(i,j)>=R(w[k],j) && X(i,j)<=R(w[k],m+j));
					}
					if (!Xb) u[k] = 0;
				}
				arma::uvec k = arma::find(u==arma::max(u));
				A[i] = w[k[0]];
			}
		}
	}
	return A+1;
}

arma::Row<double> xWeight(arma::Row<double> x){
	arma::Row<double> w = arma::exp(x - arma::max(x));
	w /= arma::sum(w);
	return w ;
}

// [[Rcpp::export]]
arma::mat dens2wght_cpp(arma::mat L){
	for(int i=0; i<L.n_rows; ++i){
		arma::Row<double> w = L.row(i);
		w = arma::exp(w - arma::max(w));
		L.row(i) = w / arma::sum(w);
	}
	return L ;
}

// [[Rcpp::export]]
double getLkh_cpp(arma::mat L){
	double lkhSum = 0.0;
	for(int i=0; i<L.n_rows; ++i){
		double mxl = arma::max(L.row(i));
		lkhSum += mxl + std::log(arma::sum(arma::exp(L.row(i)-mxl)));
	}
	return lkhSum/L.n_rows;
}
