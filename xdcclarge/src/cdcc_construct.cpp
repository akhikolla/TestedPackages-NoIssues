//' @useDynLib xdcclarge, .registration = TRUE
//' @importFrom Rcpp evalCpp

#include <stdio.h>
#include <stdlib.h>
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' This function constructs DCC-garch model's correlation term(Rt).
//' @param alpha DCC-GARCH parameter
//' @param beta DCC-GARCH parameter
//' @param stdresids matrix of standrdized(De-GARCH) residual returns (T by N)
//' @param uncR unconditional correlation matrix of stdresids (N by N)
//' @param nobs the length of time-series (T)
//' @param ndim the dimension of time-series (N)
//' @param ts how many time series are you taking
//'
//' @return DCC-garch model's correlation term(Rt)

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat cdcc_construct(double alpha, double beta, arma::mat stdresids, arma::mat uncR, int nobs, int ndim, int ts){

	int idx;
	double *rz_lag, *rtmpzz, *rdiagQ, *rzz, *rQlag, *rQbar, *rQ, *rinvdiagQ, *tmprQ, *tmprR;
	double tmp = 0.0;
	double rdab;

	rz_lag    = new double[ndim];
	rtmpzz    = new double[ndim*ndim];
	rdiagQ    = new double[ndim*ndim];
	rzz       = new double[ndim*ndim];
	rQlag     = new double[ndim*ndim];
	rQbar     = new double[ndim*ndim];
	rQ        = new double[ndim*ndim];
	rinvdiagQ = new double[ndim*ndim];
	tmprQ     = new double[ndim*ndim];
	tmprR     = new double[ndim*ndim];

	mat rDCC  = zeros(nobs, ndim*ndim);



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

	/* assigning parameters to Qbar matrix and initial values */
	for(int i=0;i<ndim;i++){
		for(int j=0;j<ndim;j++){
			rdiagQ[j + i*ndim] = 0;
			rinvdiagQ[j + i*ndim] = 0.0;
			rzz[j + i*ndim] = 0.0;
			rtmpzz[j + i*ndim] = 0.0;
			rQlag[j + i*ndim] = uncR(i,j);
			rQbar[j + i*ndim] = rdab*rQlag[j + i*ndim];
			tmprQ[j + i*ndim] = 0.0;
			tmprR[j + i*ndim] = 0.0;

		}
	}

	for(int j=0; j<ndim; j++){
			rdiagQ[j*(1+ndim)] = 1.0;       /* This is an identity matrix */
	}
	
	for(int i=0; i<nobs; i++){ // nobs

		/* outer product of z_(t-1) to create a*zz */
		for(int p=0;p<ndim;p++){
			for(int q=0;q<ndim;q++){
				tmp = alpha*rz_lag[p];
				rtmpzz[q + p*ndim] += tmp*rz_lag[q];
			} // p
		} // q

 		/* normalizing a*zz by diag(Q_t-1)^1/2 */
		/* The first dgemm computes diag(Q_t-1)^1/2*(a*zz) to be saved in rzz. */
		/* The second dgemm carries out rzz*diag(Q_t-1)^1/2 + b*Q_t-1. */
		for(int p=0;p<ndim;p++){
			for(int q=0;q<ndim;q++){
				rzz[q + p*ndim] = 0.0;
				for(int t=0;t<ndim;t++){
					rzz[q + p*ndim] += rdiagQ[t + p*ndim]*rtmpzz[q + t*ndim];
				} // t
			} // p
		} // q
		
		for(int p=0;p<ndim;p++){
			for(int q=0;q<ndim;q++){
				tmp = 0.0;
				for(int t=0;t<ndim;t++){
					tmp += rzz[t + p*ndim]*rdiagQ[q + t*ndim];
				} // t
				rQlag[q + p*ndim] = tmp + beta*rQlag[q + p*ndim];
			} // p
		} // q
		
		/* creating Q matrix for observation t */
		for(int p=0;p<ndim;p++){
			for(int q=0; q<ndim;q++){
				rQ[q + p*ndim]     = rQbar[q + p*ndim] + rQlag[q + p*ndim]; /* Qbar+(a*ee'+b*Q) */
				rQlag[q + p*ndim]  = rQ[q + p*ndim]; /* for the next round */
				rtmpzz[q + p*ndim] = 0.0;
			}
		}

		/* creating a dcc matrix for observation t */
		for(int j=0; j<ndim; j++){
			rdiagQ[(ndim+1)*j] = sqrt(rQ[(ndim+1)*j]);
			rinvdiagQ[(ndim+1)*j] = 1.0/sqrt(rQ[(ndim+1)*j]);
		}
		
		for(int p=0;p<ndim;p++){
			for(int q=0;q<ndim;q++){
				tmprQ[q + p*ndim] = 0.0;
				for(int t=0;t<ndim;t++){
					tmprQ[q + p*ndim] += rinvdiagQ[t + p*ndim]*rQ[q + t*ndim];
				} // t
			} // p
		} // q
		
		for(int p=0;p<ndim;p++){
			for(int q=0;q<ndim;q++){
				tmprR[q + p*ndim] = 0.0;
				for(int t=0;t<ndim;t++){
					tmprR[q + p*ndim] += tmprQ[t + p*ndim]*rinvdiagQ[q + t*ndim];
				} // t
			} // p
		} // q


		if(i >= (nobs-ts)){
			idx = ts-(nobs-i);
			/* tmprR must be saved */
			for(int p=0;p<ndim;p++){
				for(int q=0;q<ndim;q++){
					rDCC(idx, q + p*ndim) = tmprR[q + p*ndim];
				} // p
			} // q
		}
		for(int j=0; j<ndim; j++){
			rz_lag[j] = stdresids(i,j); /* for the next round */
		}


	} // nobs

	delete[] rtmpzz;
	delete[] rz_lag;

  return(rDCC);
}
