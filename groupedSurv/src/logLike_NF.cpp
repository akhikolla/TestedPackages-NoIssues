/*****************************/
/* Jiaxing Lin               */
/* 12-06-2017                */
/* log likelihood function   */
/*****************************/

#include <Rcpp.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppEigen)]]

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;
using namespace Rcpp;

#include "helper_NF.h"
//[[Rcpp::export]]
double logLike_NF(Rcpp::NumericVector Params, Rcpp::NumericMatrix Xmatrix,
                  Rcpp::NumericVector Kivec, Rcpp::NumericVector Deltavec,
                  int ntps) {
  // cast Rcpp Vector and Matrix to Eigen Vector and Matrix
  Eigen::Map<Eigen::VectorXd> params = as<Eigen::Map<Eigen::VectorXd>>(Params);
  Eigen::Map<Eigen::MatrixXd> xmatrix =
      as<Eigen::Map<Eigen::MatrixXd>>(Xmatrix);
  Eigen::Map<Eigen::VectorXd> kivec = as<Eigen::Map<Eigen::VectorXd>>(Kivec);
  Eigen::Map<Eigen::VectorXd> deltavec =
      as<Eigen::Map<Eigen::VectorXd>>(Deltavec);

  // coustomized vectors for computing log likelihood
  Eigen::VectorXd gammapars = params.head(ntps);
  Eigen::VectorXd gammavec(gammapars.size() + 1);
  gammavec << gammapars;
  gammavec(gammapars.size()) = INFINITY;

  Eigen::VectorXd thetavec(params.size() - ntps);
  thetavec << params.tail(params.size() - ntps);

  Eigen::VectorXd xitheta = xmatrix * thetavec;
  Eigen::VectorXd gammai(kivec.size());
  
	//for (int i = 0; i < kivec.size(); i++)
	//	gammai(i) = gammavec(kivec(i) - 1); // adjust for c++ index from zero
  

//  for (int i = 0; i < kivec.size(); i++)
//	{
// 		if(kivec(i) == INFINITY && deltavec(i) == 0) // adjust for c++ index from zero
//				kivec(i) = ntps + 1;
//  }
    //cout << "like out" <<endl;
	Eigen::VectorXd firstterm(deltavec.size());
  for (int i = 0; i < deltavec.size(); i++) {
    if (deltavec(i) == 1)
		{
	 		gammai(i) = gammavec(kivec(i) - 1);			
      firstterm(i) = log(1 - exp(-exp(gammai(i) + xitheta(i))));
	  }else{  
			firstterm(i) = 0;
		}
    //cout << "like out" <<endl;
	}
  // cout << "gammai: \n" << gammai << endl;
  // cout << "xitheta: \n " << xitheta  << endl;
  // cout << "firstterm:\n" << firstterm <<endl;

  Eigen::VectorXd rowSums(deltavec.size());
  Eigen::MatrixXd sumTerms = sumKim1mat(gammapars, xitheta, kivec, deltavec);
  // cout << "sumterms:\n" << sumTerms << endl;
  Eigen::VectorXd logLike = firstterm - sumTerms.rowwise().sum();
  // cout << "logLike for each patient:\n" <<logLike <<endl;
  // cout << logLike.sum() << endl;

  return logLike.sum();
}
