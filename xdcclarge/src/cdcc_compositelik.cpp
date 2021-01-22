//' @useDynLib xdcclarge, .registration = TRUE
//' @importFrom Rcpp evalCpp

#include <stdio.h>
#include <stdlib.h>
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' This function calculates log-likelihood of cDCC-garch model with composite likelihood method.
//' @param alpha cDCC-GARCH parameter
//' @param beta cDCC-GARCH parameter
//' @param ht matrix of conditional variance vectors (T by N)
//' @param residuals matrix of residual(de-mean) returns  (T by N)
//' @param stdresids matrix of standrdized(De-GARCH) residual returns (T by N)
//' @param uncR unconditional correlation matrix of stdresids (N by N)
//' @param nobs the length of time-series (T)
//' @param ndim the dimension of time-series (N)
//'
//' @return log-likelihood of cDCC-GARCH model(scaler)

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double cdcc_compositelik(double alpha, double beta, arma::mat ht, arma::mat residuals, arma::mat stdresids, arma::mat uncR, int nobs, int ndim){


	// cDCC calc params
	double *rz_lag;
	double tmp = 0.0;
	double rdab;

	rz_lag    = new double[ndim];

	mat hh = zeros(ndim,ndim);
	mat St = zeros(ndim,ndim);
	mat vdvar = zeros(ndim,ndim);

	// new params
	double *nrQbar, *nrQlag, *nrdiagQ;
	double ns11;
	double ns12;
	double ns22;
	double nx11;
	double nx12;
	double nx22;
	double nsdel;
	double nlf_sub;
	double nlfs = 0.0;

	nrQbar     = new double[(ndim-1)*3];
	nrQlag     = new double[(ndim-1)*3];
	nrdiagQ    = new double[(ndim-1)*2];

	double nrinvdiagQ1 = 0.0;
	double nrinvdiagQ2 = 0.0;
	double nrQ1        = 0.0;
	double nrQ2        = 0.0;
	double nrQ3        = 0.0;

	////////////////////////////////////
	// Initialize params
	////////////////////////////////////

	/* a coefficient on rQbar */
	rdab = 1.0 - alpha - beta;

	/* initial values = the average of eps^2 */
	for(int j=0; j<ndim; j++){
		rz_lag[j] = 0.0;
		for(int i=0; i<nobs; i++){
			rz_lag[j] += stdresids(i,j);
		}
		rz_lag[j] = rz_lag[j]/nobs ;
	}

	for(int i=0;i<ndim-1;i++){
		nrQbar[i*3]       = rdab*uncR(i,i);
		nrQbar[i*3+1]     = rdab*uncR(i,i+1);
		nrQbar[i*3+2]     = rdab*uncR(i+1,i+1);
		nrQlag[i*3]       = uncR(i,i);
		nrQlag[i*3+1]     = uncR(i,i+1);
		nrQlag[i*3+2]     = uncR(i+1,i+1);
		nrdiagQ[i*2]      = 1.0;
		nrdiagQ[i*2+1]    = 1.0;
	}


	////////////////////////////////////
	// Main routine
	////////////////////////////////////

	for(int T=0; T<nobs; T++){ // nobs

		// composite likelihood method
		for(int i=0;i<ndim;i++){
			for(int j=0;j<ndim;j++){
				hh(i,j) = sqrt(ht(T,i))*sqrt(ht(T,j));
				vdvar(i,j) = residuals(T,i)*residuals(T,j);

			}
		}

		nlf_sub = 0.0;
		for(int p=0;p<ndim-1;p++){
			tmp            = alpha*rz_lag[p];
			nrQlag[p*3]    = (nrdiagQ[p*2]*(tmp*rz_lag[p]))*nrdiagQ[p*2] + beta*nrQlag[p*3];
			nrQlag[p*3+1]  = (nrdiagQ[p*2]*(tmp*rz_lag[p+1]))*nrdiagQ[p*2+1] + beta*nrQlag[p*3+1];
			nrQ1           = nrQbar[p*3] + nrQlag[p*3];
			nrQ2           = nrQbar[p*3+1] + nrQlag[p*3+1];
			tmp            = alpha*rz_lag[p+1];
			nrQlag[p*3+2]  = (nrdiagQ[p*2+1]*(tmp*rz_lag[p+1]))*nrdiagQ[p*2+1] + beta*nrQlag[p*3+2];
			nrQ3           = nrQbar[p*3+2] + nrQlag[p*3+2];

			nrQlag[p*3]    = nrQ1;
			nrQlag[p*3+1]  = nrQ2;
			nrQlag[p*3+2]  = nrQ3;

			nrdiagQ[p*2]   = sqrt(nrQ1);
			nrdiagQ[p*2+1] = sqrt(nrQ3);
			nrinvdiagQ1    = 1.0/sqrt(nrQ1);
			nrinvdiagQ2    = 1.0/sqrt(nrQ3);

			// Composite Likelihood

			ns11 = hh(p,p)*(nrinvdiagQ1*nrQ1)*nrinvdiagQ1;
			ns12 = hh(p,p+1)*(nrinvdiagQ1*nrQ2)*nrinvdiagQ2;
			ns22 = hh(p+1,p+1)*(nrinvdiagQ2*nrQ3)*nrinvdiagQ2;

			nsdel = ns11*ns22 - ns12*ns12;

			nx11 = vdvar(p,p);
			nx12 = vdvar(p,p+1);
			nx22 = vdvar(p+1,p+1);

			nlf_sub += 0.5*(log(nsdel) + ((ns22*nx11 - 2*ns12*nx12 + ns11*nx22)/nsdel))/nobs;

		}
		nlfs += nlf_sub;


		for(int j=0; j<ndim; j++){
			rz_lag[j] = stdresids(T,j); /* for the next round */
		}

	} // nobs

	delete[] rz_lag;

    return(nlfs);
}
