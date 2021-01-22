/*****************************/
/* Jiaxing Lin               */
/* 05-10-2016                */
/* Main function to find the */
/* sigm for the maximum like-*/
/* lihood function for multi-*/
/* families.                 */
/*****************************/
//[[Rcpp::depends(BH)]]

#define _USE_MATH_DEFINES
#define BOOST_BIND_NO_PLACEHOLDERS

#include <iostream>
#include <math.h>
#include "global.h"
#include <boost/math/tools/minima.hpp>
#include "fam_LLVar.h"
#include "Rcpp.h"
#include "famSize.h"
#include <functional> // functional

using namespace std::placeholders;
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
double varEst( std::vector<std::string > fam_group,
             Rcpp::NumericVector alpha,
             Rcpp::NumericVector dtime,
             Rcpp::NumericVector delta,
             Rcpp::NumericVector g,       
             double beta,  // The value of beta
             double lower, // Lower bound of opt regime
             double upper,  // Upper bound of opt regime
			 std::vector<std::string > f_ind,
			 int m) 
{
  double* av = new double[alpha.size()];
  int*    dt = new int[dtime.size()];
  int* Delta = new int[delta.size()];
  double*     G = new double[g.size()];
  double* logat = new double[alpha.size()];
  int*  famsize = new int[m];
 
  for(int i=0; i<alpha.size(); i++)
    av[i]=alpha[i];
  for(int i=0; i<alpha.size(); i++)
    logat[i]=log(av[i]);
  for(int i=0; i<dtime.size(); i++)
    dt[i]=dtime[i];
  for(int i=0; i<g.size(); i++)
    G[i]=g[i];
  for(int i=0; i<delta.size(); i++)
    Delta[i]=delta[i];  	 

  // generate number of subject for each family and put in famsize array.
  famSize(famsize, fam_group, fam_group.size());
 
  int    dttmp[4] = {0};
  int Deltatmp[4] = {0};
  double  Gtmp[4] = {0};
   
  global_alpha_v_       = av;       
  global_Dtime_         = dttmp;   
  global_Delta_         = Deltatmp;
  global_G_             = Gtmp;
  global_beta_          = &beta; 
  global_log_alpha_v_   = logat;

  typedef std::pair<double, double> Result;
  Result res = boost::math::tools::brent_find_minima(
              std::bind(fam_LLVar, _1, famsize, dt, Delta, G, m, f_ind),
              lower, upper, 10);
  double min = res.first;
	
  return min;
  /*
   return the max likelihood and 
   corresponding values for beta
   */
}




