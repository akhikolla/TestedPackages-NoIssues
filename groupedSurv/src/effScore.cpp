/*****************************/
/* Jiaxing Lin               */
/* 12-06-2017                */
/* log likelihood function   */
/*****************************/

#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppEigen)]]

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;
using namespace Rcpp;

#include "effScore.h"

//[[Rcpp::export]]
List effScore(double beta, Rcpp::NumericVector Params, Rcpp::NumericMatrix G,
                  Rcpp::NumericMatrix Xmatrix, Rcpp::NumericVector Kivec, 
                  Rcpp::NumericVector Deltavec, int ntps, int nCores, bool reScore) {
  // cast Rcpp Vector and Matrix to Eigen Vector and Matrix
  Eigen::Map<Eigen::VectorXd> params = 
      as<Eigen::Map<Eigen::VectorXd>>(Params);
  Eigen::Map<Eigen::MatrixXd> xmatrix =
      as<Eigen::Map<Eigen::MatrixXd>>(Xmatrix);
  Eigen::Map<Eigen::VectorXd> kivec = 
      as<Eigen::Map<Eigen::VectorXd>>(Kivec);
  Eigen::Map<Eigen::VectorXd> deltavec =
      as<Eigen::Map<Eigen::VectorXd>>(Deltavec);
  Eigen::Map<Eigen::MatrixXd> g =
      as<Eigen::Map<Eigen::MatrixXd>>(G);

  Eigen::VectorXd scoreStatic(g.cols());
  scoreStatic.fill(0);
  Eigen::MatrixXd UiFull(g.rows(), g.cols());
   
  #ifdef _OPENMP
  omp_set_num_threads(nCores); // define number of threads.	
  #endif
  #ifdef _OPENMP
  #pragma omp parallel for 
  #endif
   
  for(int i = 0; i < g.cols(); i++)
  {  
     Eigen::VectorXd g_tmp = g.col(i);
     Eigen::VectorXd USNPs(g_tmp.size());
     
     // default to return statistics only		 
     scoreStatic(i) = effScore_NF_S(beta, params, g_tmp, xmatrix, kivec, deltavec, ntps, reScore, USNPs);
     
		 // options for return score matrix
		 if(reScore)
        UiFull.col(i)  = USNPs; 
  }

  Rcpp::NumericVector sStatics(wrap(scoreStatic)); 
  
	if(reScore) 
  {	Rcpp::NumericMatrix uScore(wrap(UiFull));
    return List::create(
                        Named("sStatics") = sStatics,
                        Named("uScore")  = uScore
                     );
  }else{
  	return List::create(Named("sStatics") = sStatics);
                    
  }
}
