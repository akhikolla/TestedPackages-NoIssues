/*
 This is the routine calling the 
 Cuba Cuhre doing the integration
 for the likelihood function for
 mutiple trio families.
 Jiaxing 20170308
 */

#include <math.h>
#include "cuba.h"
#include "integrand_est.h"
#include "global.h"
#include <iostream>
#include <vector>

#define NCOMP 1
#define PI  M_PI

double fam_LLVar(double sigma2, int fam_size[], int dt[], 
                 int Delta[], double G[], int m, vector<string> &f_ind)
{
  
  int flag = 0;
  global_2off_flag_ = &flag;
  global_sigma2_ = &sigma2;
  
  int curPid = 0;
  int  nregions, neval, fail;
  double integral[NCOMP]={0}, error[NCOMP]={0}, prob[NCOMP]={0}; 
  double sum_LL = 0.00;
  
  for(int k = 0; k < m; k++){
    if(fam_size[k] == 2 && f_ind[curPid].compare(f_ind[curPid+1]) == 0)
      flag = 1;
    for(int j = 0; j < fam_size[k]; j++){
      global_Dtime_[j] = dt[j+curPid];
      global_G_[j]     = G[j+curPid];
      global_Delta_[j] = Delta[j+curPid];
    }
    
    curPid +=fam_size[k]; // update the current PId.
    if(fam_size[k] != 1)
     Cuhre(fam_size[k], NCOMP,
            Integrand_est,
            1e-10, 1e-15,
            0, 0, 100000,
            0,
            &nregions, &neval, &fail,
            integral, error, prob
            );
    
    switch(fam_size[k])
    {      
      case 4: 
        sum_LL += log((double)integral[0]
                      /(*global_sigma2_)/(*global_sigma2_)
                      /(2*PI)/(2*PI)/.5
                    );
        break;
      case 3:
        sum_LL += log((double)integral[0]
                      /(*global_sigma2_)/sqrt((*global_sigma2_))
                      /(2*PI)/sqrt(2*PI)/sqrt(0.5)
                      );
        break;
      case 2:
		if(flag == 1)
        sum_LL += log((double)integral[0]/(*global_sigma2_)/(2*PI));
        else
        sum_LL += log((double)integral[0]/(*global_sigma2_)/(2*PI)/sqrt(.75));
        break;
		  case 1:
				flag = 0;
				sum_LL += log((double)integral[0]/sqrt(*global_sigma2_)/sqrt(2*PI));
		}
	
	//if(fam_size[k] == 1)
	//	sum_LL += log(IntInd()/sqrt(*global_sigma2_)/sqrt(2*PI));
   // flag = 0;
  }
  return -sum_LL; 		// return the negative log likelihood 
}


