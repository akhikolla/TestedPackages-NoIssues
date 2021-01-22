#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
NumericVector sample_rcpp(int N, int nsamp) {
	// Inputs: 
	// N: Largest integer to sample from.
	// nsamp: number of samples from [1,2,...,N] with replacement to obtain.
	// Output: 
	// samps: nsamp-length vector of samples from [1,2,...,N] with replacement to obtain.
	NumericVector samps = ceiling(runif(nsamp)*N);
	return samps;
}

//[[Rcpp::export]]
NumericVector Cquantile(NumericVector xx, NumericVector p) {
	NumericVector x = clone(xx);
	x = x.sort();
	int n = x.size();
	int n2 = p.size();
	NumericVector m = 1 - p; 
	NumericVector j = floor(n*p + m);
	NumericVector g = n*p + m - j;
	NumericVector qtls(n2);
	for(int i=0; i < n2; i++){
		qtls[i] = (1-g[i])*x[j[i]-1] + g[i]*x[j[i]];
		}
	return qtls;
}

//[[Rcpp::export]]
List Cdboot_multi(NumericMatrix xxyy, NumericVector lgridlo, NumericVector lgridhi, int B, int B2, int G){
	// Inputs: 
	// xxyy: (n by p+1) matrix for X (design matrix) and response vector y.
	// lgridlo: Lower quantile values of double bootstrap distribution to obtain
	// lgridhi: Upper quantile values of double bootstrap distribution to obtain
	// B: number of first-level bootstrap samples
	// B2: number of double bootstrap samples
	// G: calculate quantile-based empirical coverage at this many grid points
	
	// Outputs: 
	// theta_hat_boot: first-level bootstrap estimates of all slope coefficients
	// theta_qtl_lgrid_lo: (p+1 by B by G by 1) matrix 
	// theta_qtl_lgrid_hi: 
		
	// Cast NumericMatrix etc. to arma types
	mat xy = as<arma::mat>(xxyy);
	
	// Define stuff
	int n = xy.n_rows;
	int nc = xy.n_cols;
	
	xy = join_rows(ones(n,1),xy);

	// Initialize stuff to be used within for loops below
	arma::vec inds(n);
	inds.zeros();
	arma::mat xy_boot(n,nc);
	xy_boot.zeros();
	arma::mat xy_dboot(n,nc);
	xy_dboot.zeros();
	arma::mat theta_hat_boot(B,nc);	
	theta_hat_boot.zeros();
	arma::mat xy_boot2(n,nc);
	xy_boot2.zeros();
	arma::mat Resid_boot(n,nc);	
	Resid_boot.zeros();
	arma::mat theta_hat_dboot(B2,nc);
	theta_hat_dboot.zeros();
	arma::cube theta_qtl_lgrid_lo_element(B, G, 1);
	theta_qtl_lgrid_lo_element.zeros();
	arma::cube theta_qtl_lgrid_hi_element(B, G, 1);
	theta_qtl_lgrid_hi_element.zeros();
	
	arma::field<cube> theta_qtl_lgrid_lo(nc,1);
	for (int idx=0; idx<nc; idx++){
		theta_qtl_lgrid_lo(idx,0) = theta_qtl_lgrid_lo_element;	
		}
	arma::field<cube> theta_qtl_lgrid_hi(nc,1);
	for (int idx=0; idx<nc; idx++){
		theta_qtl_lgrid_hi(idx,0) = theta_qtl_lgrid_hi_element;	
		}

	arma::uvec y(n);
	arma::uvec y2(n);
	IntegerVector frame = seq_len(n);
	arma::mat X_boot(n,nc);	
	X_boot.zeros();
	arma::mat y_boot(n,1);
	y_boot.zeros();
	arma::mat X_dboot(n,nc);
	X_dboot.zeros();
	arma::mat y_dboot(n,1);
	y_dboot.zeros();
	arma::colvec betahat_boot(nc);	
	betahat_boot.zeros();
	arma::colvec betahat_dboot(nc);	
	betahat_dboot.zeros();
	arma::vec qtl_info(G);	
	qtl_info.zeros();
	arma::vec qtl_info2(G);
	qtl_info2.zeros();
	arma::vec theta_dboot_avg(nc);	
	theta_dboot_avg.zeros();
	arma::mat theta_hat_dboot_onecol(B2,1);
	theta_hat_dboot_onecol.zeros();

	
	for(int j = 0; j < B; j++){
		y = as<arma::uvec>(sample_rcpp(n,n))-1;
		xy_boot = xy.rows(y);
		X_boot = xy_boot.cols(0,nc-1);
		y_boot = xy_boot.col(nc);
		betahat_boot = arma::solve(X_boot, y_boot);
		for(int k=0; k < nc; k++){
			theta_hat_boot(j,k) = betahat_boot(k,0);
			}
		
		theta_dboot_avg.zeros();
		for(int k = 0; k < B2; k++){
			y2 = as<arma::uvec>(sample_rcpp(n,n))-1;
			xy_dboot = xy_boot.rows(y2);
			X_dboot = xy_dboot.cols(0,nc-1);
			y_dboot = xy_dboot.col(nc);
			betahat_dboot = arma::solve(X_dboot, y_dboot);
			for(int l=0; l < nc; l++){
				theta_hat_dboot(k,l) = betahat_dboot(l,0);
				theta_dboot_avg(l) += theta_hat_dboot(k,l);
				}
			}
		
		for(int l=0; l < nc; l++){
			theta_dboot_avg(l) = theta_dboot_avg(l) / B2;
			}
					
		for(int m=0; m<nc; m++){
			theta_hat_dboot_onecol = theta_hat_dboot.col(m);			
			qtl_info = as<arma::vec>(Cquantile(as<NumericVector>(wrap(theta_hat_dboot_onecol)),lgridlo));
			qtl_info2 = as<arma::vec>(Cquantile(as<NumericVector>(wrap(theta_hat_dboot_onecol)),lgridhi));				
			for(int l=0; l < G; l++){
				theta_qtl_lgrid_lo(m,0)(j,l,0) = qtl_info(l);
				theta_qtl_lgrid_hi(m,0)(j,l,0) = qtl_info2(l);
				}
			}
		}
		
	// Create output:
	List L = List::create(Named("theta_hat_boot")=theta_hat_boot, 
					Named("theta_qtl_lgrid_lo")=theta_qtl_lgrid_lo, Named("theta_qtl_lgrid_hi")=theta_qtl_lgrid_hi);
	return L;
}
