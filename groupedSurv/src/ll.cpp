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
#include <iostream>
#include <math.h>
#include "global.h"
#include <boost/bind.hpp>
#include "fminbr.h"
#include <ctime>
#include "fam_LL.h"
#include "Rcpp.h"
#include "famSize.h"

using boost::bind;
typedef boost::function<double(double x)> bindtype;

using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
double  ll(  std::vector<std::string > fam_group,
             Rcpp::NumericVector alpha,
             Rcpp::NumericVector dtime,
             Rcpp::NumericVector delta,
             Rcpp::NumericVector g,       
             double beta,  // The value of beta
             double var, // Lower bound of opt regime
			 std::vector<std::string > f_ind, int m) 
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
   
  global_alpha_v_       = av;       
  global_Dtime_         = dt;   
  global_Delta_         = Delta;
  global_G_             = G;
  global_beta_          = &beta; 
  global_log_alpha_v_   = logat;
  
  return fam_LL(var, famsize, dt, Delta, G, m, f_ind);
}




