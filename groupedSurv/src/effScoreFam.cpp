#include <stdio.h>
#include <math.h>
#include <vector>
#include "cuba.h"
#include "Rcpp.h"
#include "global.h"
#include "integrand.h"
#include "famSize.h"
#include <RcppEigen.h>
//#include <gsl/gsl_integration.h>
#include "fsigma.h"
#include "falpha.h"
#include "fbeta.h"

//[[Rcpp::depends(RcppEigen)]]

using Eigen::VectorXd;
using Eigen::MatrixXd;

using namespace Rcpp;
using namespace std;

#define NCOMP 13
#define COM 3

//[[Rcpp::export]]
List effScoreFam(double beta, double sigma2,
                   std::vector<std::string > fam_group,
                   Rcpp::NumericVector i_at,
                   Rcpp::NumericVector i_dt,
                   Rcpp::NumericVector i_Delta,
                   Rcpp::NumericVector i_G,
				   std::vector<std::string > f_ind, int m)
{
  int alpha_size = i_at.size();
  
  //Create two set of variables that will be used in the 
  //integratation part for both likelihood and its partial derivatives
  //with respect to beta.
  //int comp,
  int nregions, neval, fail;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];

  // Convert R input objects into c++ objects.
  double* at	= new double[i_at.size()];
  double* logat = new double[i_at.size()];
  int*    dt 	= new int[i_dt.size()];
  int* Delta 	= new int[i_Delta.size()];
  double*  G  	= new double[i_G.size()];
  int*  famsize = new int[m];
 
  for(int i=0; i<i_at.size(); i++)
				at[i]=i_at[i];	
  for(int i=0; i<i_at.size(); i++)
        logat[i]=log(i_at[i]);
  for(int i=0; i<i_dt.size(); i++)
        dt[i]=i_dt[i];
  for(int i=0; i<i_G.size(); i++)
        G[i]=i_G[i];
  for(int i=0; i<i_Delta.size(); i++)
        Delta[i]=i_Delta[i];

  // automatic generate the famSize array
  famSize(famsize, fam_group, fam_group.size());
		
  // handle global variables
  global_alpha_v_		= at;
  global_log_alpha_v_	= logat;
  global_Dtime_			= dt;
  global_G_				= G;
  global_Delta_         = Delta;
  global_beta_          = &beta;
  global_sigma2_        = &sigma2;
  global_alpha_size_	= &alpha_size;

  // Create some verctors to store the results.
  double l_i = 0;
  //double pdevm[NCOMP]={0}; 
  

  // U_Eta vector
  VectorXd U_Eta(alpha_size+1);
  // C_Beta_Eta vector 
  VectorXd C_Beta_Eta(alpha_size+1);
  C_Beta_Eta.fill(0);
  // V_Eta_Eta matrix
  MatrixXd V_Eta_Eta(alpha_size+1, alpha_size+1);
  V_Eta_Eta.fill(0);
  
  NumericVector betaScore(m);
  NumericMatrix etaScore(alpha_size+1,m);

	
  int alpha_scoreindex = 0;
  int flag = 0;
  global_2off_flag_ = &flag;
  global_alpha_index_ = &alpha_scoreindex;
  int curPid = 0;
  // Do the integration.
  for(int k = 0; k < m; k++){
	if(famsize[k] == 2 && f_ind[curPid].compare(f_ind[curPid+1]) == 0)
        flag = 1;
	for(int j = 0; j < famsize[k]; j++){
      global_Dtime_[j] = dt[j+curPid];
      global_G_[j]     = G[j+curPid];
      global_Delta_[j] = Delta[j+curPid];
  	}
    // shift the pid position using famsize.
	curPid += famsize[k];
    // compute the effective score for each family.
	if(famsize[k] >= 1)
	{        
		Cuhre(famsize[k], NCOMP,
          Integrand,
          1e-6, 1e-12,
          0, 0, 50000,
          0,
          &nregions, &neval, &fail,
          integral, error, prob);

		flag = 0;
	}
  // change computing of integral for one member family to use cuhre to 
	// avoid problem with windows
	/*else
	{
		double ll = 0, dfbeta = 0, dfalpha = 0, dfsigma = 0;
		gsl_integration_workspace * w
    	= gsl_integration_workspace_alloc (10000);					
   		double result, error=0.0001;
	    double alpha = 1.0;
	
		// compute the likelihood for individual case
    	gsl_function F;
    	F.function = &fll;
    	F.params = &alpha;
	
		gsl_integration_qag (&F,-30*sqrt(sigma2),
                    30*sqrt(sigma2), 1e-6, 1e-4, 10000,1,
           			w, &ll, &error);	
    	//gsl_integration_qagi (&F, 0, 1e-7, 10000,
        //	w, &ll, &error);
		// compute the partial derivative of likelihood regarding to \sigma^2.
		F.function = &fsigma;
  	    gsl_integration_qag (&F,-30*sqrt(sigma2),
                    30*sqrt(sigma2), 1e-6, 1e-4, 10000,1,
           			w, &dfsigma, &error);


		//gsl_integration_qagi (&F, 0, 1e-7, 10000,
        //    w, &dfsigma, &error);
		// compute the partial derivative of likelihood regarding to \beta.
		F.function = &fbeta;
		//gsl_integration_qagi (&F, 0, 1e-7, 10000,
        //    w, &dfbeta, &error);
        gsl_integration_qag (&F,-30*sqrt(sigma2),
                    30*sqrt(sigma2), 1e-6, 1e-4, 10000,1,
           			w, &dfbeta, &error);

		// compute the partial derivative of likelihood regarding to \alpha
		//F.function = &falpha;
        //gsl_integration_qagi (&F, 0, 1e-7, 10000,
        //    w, &dfalpha, &error);
		

		integral[0] = ll;
		integral[1] = dfsigma;
		integral[2] = dfbeta; 
		
		for(int alpha_i = 0; alpha_i < alpha_size; alpha_i++)
		{
			if(alpha_i > global_Dtime_[0]-1)
				integral[3+alpha_i] = 0;	
		 	else{
				F.function = &falpha;
       	 		alpha_scoreindex = alpha_i;
				gsl_integration_qag (&F,-30*sqrt(sigma2),
                    30*sqrt(sigma2), 1e-6, 1e-4, 10000,1,
           			w, &dfalpha, &error);
				integral[3+alpha_i] = dfalpha;
			}
		}
		gsl_integration_workspace_free (w);
	}	
  */
	// compute the partial derivtive of likelihood for each family	
	l_i +=log(integral[0]);	

	// compute the score for individual family
    betaScore(k) = integral[2]/integral[0];  //beta
	  etaScore(0, k) = integral[1]/integral[0]; //sigma^2
  	for(int i = 0; i < alpha_size; i++) 	
		etaScore(i+1,k) = integral[3+i]/integral[0]; //alpha 
  	
	// compute the terms for efficient scores
	// U_beta terms
	double U_Beta = integral[2]/integral[0];
	// U_Eta terms
	U_Eta[0] = integral[1]/integral[0]; // sigma^2
	for(int i = 0; i < alpha_size; i++)
		U_Eta[1+i] = integral[3+i]/integral[0]; //alpha
	
	// CbetaEta and VEtaEta are obtained from means for all the families
	// (sum up here and take mean at the end).    
  	C_Beta_Eta += U_Beta*U_Eta;
  	V_Eta_Eta += U_Eta*U_Eta.transpose();
  }	
 
  NumericVector C_Beta_Eta_out(C_Beta_Eta.size()-1);
  for(int i = 0; i < C_Beta_Eta.size()-1; i++)
     C_Beta_Eta_out[i] = C_Beta_Eta[i]/m;
  
	NumericMatrix V_Eta_Eta_out(C_Beta_Eta.size()-1,C_Beta_Eta.size()-1);
  for(int i = 0; i < C_Beta_Eta.size()-1; i++)
	  for(int j = 0; j < C_Beta_Eta.size()-1; j++)
     V_Eta_Eta_out(i,j) = V_Eta_Eta(i,j)/m;
  
  NumericMatrix etaScore_out(alpha_size, m);
  for(int i = 0; i < alpha_size; i++)
	  for(int j = 0; j < m; j++)
     etaScore_out(i,j) = etaScore(i,j);

  //return results as list 
  return List::create(
  			    Named("LogLik") = l_i,
						Named("beta_score")  = betaScore,
						Named("eta_score") = etaScore_out,
						Named("C_Beta_Eta")  = C_Beta_Eta_out,
						Named("V_Eta_Eta")  = V_Eta_Eta_out   
					 );
  
}


