/*
 Jiaxing Lin
 Dec 7 2017
 */

/*
 Functions for repeated terms in partial derivatives
 Sum i=1 to n of Sum j=1 to k_i-1
 x_i and g_i terms can be pulled out
 */
#ifndef __helper_ef_h__
#define __helper_ef_h__

#include <math.h>
#include <limits>
/*static Eigen::MatrixXd sumKim1mat(const VectorXd &gammavec,
                                  const VectorXd &xitheta,
                                  const VectorXd &kivec, const VectorXd &deltavec) {
  Eigen::MatrixXd sumterms(xitheta.size(), gammavec.size());
  Eigen::VectorXd kim1(kivec.size());
  kim1.fill(-1);
  kim1 = kim1 + kivec;

	for (int i = 0; i < kivec.size(); i++)
  {
    if(kivec(i) == INFINITY && deltavec(i) == 0) // adjust for c++ index from zero
       kim1(i) = gammavec.size() ;
  }

  for (int i = 0; i < xitheta.size(); i++) {
    for (int j = 0; j < gammavec.size(); j++) {
      if (j >= kim1(i) || gammavec(j) == -INFINITY || gammavec(j) == INFINITY)
        sumterms(i, j) = 0;
      else
        sumterms(i, j) = exp(xitheta(i) + gammavec(j));
    }
  }
  return sumterms;
}

// delta_i * fraction of exp()s
static VectorXd deltaifrac(const VectorXd &gammavecNI, const VectorXd &xitheta,
                           const VectorXd &kivec, const VectorXd &deltavec) {
  // add INFINITY as the last element of gammavec
  Eigen::VectorXd gammavec(gammavecNI.size() + 1);
  gammavec << gammavecNI;
  gammavec(gammavecNI.size()) = INFINITY;

  // get the last gamma value for each sample.
  Eigen::VectorXd gammai(kivec.size());
  for (int i = 0; i < kivec.size(); i++)
    gammai(i) = gammavec(kivec(i) - 1);

  // get the exp() terms
  Eigen::VectorXd ex(kivec.size());
  Eigen::VectorXd exex(kivec.size());
  Eigen::VectorXd out(kivec.size());
  for (int i = 0; i < kivec.size(); i++) {
		if(deltavec(i) == 0)
		{
			out(i) = 0;
			continue;
		}
    if(gammai(i) == INFINITY||gammai(i) == -INFINITY)
    {
      out(i) = 0;
      continue;
    }
    ex(i) = exp(gammai(i) + xitheta(i));
    exex(i) = exp(-ex(i));
    //if (deltavec(i) == 1)
      out(i) = exex(i) * ex(i) / (1 - exex(i));
    //else
      //out(i) = 0;
  }
  return out;
}
*/
// delta_i * fraction of exp()s for Eff score
inline VectorXd deltaifrac_EF(const VectorXd &gammavecNI, const VectorXd &gib, const VectorXd &xitheta,
                           const VectorXd &kivec, const VectorXd &deltavec) {
  long long inf = std::numeric_limits<long>::infinity();
    // add INFINITY as the last element of gammavec
  Eigen::VectorXd gammavec(gammavecNI.size() + 1);
  gammavec << gammavecNI;
  gammavec(gammavecNI.size()) = inf;

  // get the last gamma value for each sample.
  Eigen::VectorXd gammai(kivec.size());
  for (int i = 0; i < kivec.size(); i++)
    gammai(i) = gammavec(kivec(i) - 1);

  // get the exp() terms
  Eigen::VectorXd ex(kivec.size());
  ex.fill(0);
  Eigen::VectorXd exex(kivec.size());
  exex.fill(0);
  Eigen::VectorXd out(kivec.size());
  out.fill(0);
  for (int i = 0; i < kivec.size(); i++) {
		if(deltavec(i) == 0)
		{
			out(i) = 0;
			continue;
		}
    if(gammai(i) == inf||gammai(i) == -inf)
    {
      out(i) = 0;
      continue;
    }

    ex(i) = exp(gammai(i) + gib(i)+ xitheta(i));
    if (deltavec(i) == 1)
      out(i) =  ex(i) / ( exp(ex(i)) - 1);
    else
      out(i) = 0;
    
  }
  return out;
}



inline Eigen::MatrixXd sumKim1mat_EF(const VectorXd &gammavec,const VectorXd &gib,
                                  const VectorXd &xitheta,
                                  const VectorXd &kivec, const VectorXd &deltavec) {
  Eigen::MatrixXd sumterms(xitheta.size(), gammavec.size());
  //  cout << "gammavec.size: \n" << gammavec.size() <<endl;
  Eigen::VectorXd kim1(kivec.size());
  kim1.fill(-1);
  kim1 = kim1 + kivec;
  long long inf = std::numeric_limits<long>::infinity();
  for (int i = 0; i < kivec.size(); i++)
  {
    if(kivec(i) == inf && deltavec(i) == 0) // adjust for c++ index from zero
       kim1(i) = gammavec.size();
  }
  
  for (int i = 0; i < xitheta.size(); i++) {
    for (int j = 0; j < gammavec.size(); j++) {
      if (j >= kim1(i) || gammavec(j) == -inf || gammavec(j) == inf)
        sumterms(i, j) = 0;
      else
        sumterms(i, j) = exp(xitheta(i) + gib(i) + gammavec(j));
    }
  }
  return sumterms;
}





#endif
