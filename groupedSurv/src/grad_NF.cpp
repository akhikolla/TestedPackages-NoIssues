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
Rcpp::NumericVector grad_NF(Rcpp::NumericVector Params,
                            Rcpp::NumericMatrix Xmatrix,
                            Rcpp::NumericVector Kivec,
                            Rcpp::NumericVector Deltavec, int ntps) {
  // map Rcpp Vector and Matrix to Eigen Vector and Matrix
  // without using extra memory
  Eigen::Map<Eigen::VectorXd> params = as<Eigen::Map<Eigen::VectorXd>>(Params);
  Eigen::Map<Eigen::MatrixXd> xmatrix =
      as<Eigen::Map<Eigen::MatrixXd>>(Xmatrix);
  Eigen::Map<Eigen::VectorXd> kivec = as<Eigen::Map<Eigen::VectorXd>>(Kivec);
  Eigen::Map<Eigen::VectorXd> deltavec =
      as<Eigen::Map<Eigen::VectorXd>>(Deltavec);

	// gammavec
  Eigen::VectorXd gammavec = params.head(ntps);

  // thetavec
  Eigen::VectorXd thetavec = params.tail(params.size() - ntps);

  // xitheta
  Eigen::VectorXd xitheta = xmatrix * thetavec;

	//cout << "grad out" <<endl;
  // summat
  Eigen::MatrixXd summat = sumKim1mat(gammavec, xitheta, kivec, deltavec);

  // dfrac
  Eigen::VectorXd dfrac = deltaifrac(gammavec, xitheta, kivec, deltavec);

  // wrt \thea
  Eigen::MatrixXd grad_thetaTemp(xmatrix.rows(), xmatrix.cols());
  Eigen::VectorXd temp = dfrac - summat.rowwise().sum();
  for (int i = 0; i < xmatrix.rows(); i++)
    for (int j = 0; j < xmatrix.cols(); j++)
      grad_thetaTemp(i, j) = xmatrix(i, j) * temp(i);
  Eigen::VectorXd grad_theta = grad_thetaTemp.colwise().sum();

  // wrt \gamma
  Eigen::VectorXd grad_alpha = -summat.colwise().sum();
  for (int i = 0; i < kivec.size(); i++)
	{  
		if(deltavec(i) != 0){
 			grad_alpha(kivec(i) - 1) =
      	  grad_alpha(kivec(i) - 1) + dfrac(i); // taken care of ki in C++ and R
    }
	}
	//cout << "grad out" <<endl;
  // concatenate
  Eigen::VectorXd grad_Xd(grad_theta.size() + grad_alpha.size());
  grad_Xd << grad_alpha, grad_theta;

  Rcpp::NumericVector grad(wrap(grad_Xd));

  return grad;
}
