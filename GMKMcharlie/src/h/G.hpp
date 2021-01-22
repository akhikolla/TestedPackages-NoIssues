# pragma once
# include <Rcpp.h>
# include "macros.hpp"


// gaussian kernel class
template<typename indtype, typename valtype>
struct G
{
  bool updateAlpha, updateMean, updateSigma;
  // valtype piConstMul;
  valtype alpha; // Weight on this kernel.
  valtype sqrtOfDet; // Square root of the covariance matrix's determinant.
  vec<valtype> mu; // Kernel's center.
  vec<valtype> cholU; // Upper triangle matrix from Cholesky decomposition of the covariance matrix.
  vec<valtype> ptr; // Densities this kernel projects at a set of points in space. Could be an empty vector.


  // G(){}
  // G(indtype d) { piConstMul = std::pow(2.0 * M_PI, d * (-0.5)); alpha = 1; }


  inline bool computeCholUandSqrtOfDet(vec<valtype> &sigma);
  inline bool computeCholUandSqrtOfDet(valtype *sigma);


  inline void shallowCopy(G &x) // Copy parameters of the kernel but not the densities it projects.
  {
    updateAlpha = x.updateAlpha;
    updateMean = x.updateMean;
    updateSigma = x.updateSigma;
    alpha = x.alpha;
    sqrtOfDet = x.sqrtOfDet;
    mu = x.mu;
    cholU = x.cholU;
  }


  // For debugging
  inline void print(bool convertToSigma = 1)
  {
    Rcpp::Rcout << alpha << "    ";
    for(indtype i = 0, iend = mu.size(); i < iend;++i)
    {
      Rcpp::Rcout << mu[i] <<" ";
    }

    Rcpp::Rcout << "    ";
    if(convertToSigma)
    {
      vec<valtype> tmp(mu.size() * mu.size());
      vec<valtype*> tmp2(mu.size());
      triCholToFullSigma(&cholU[0], &tmp[0], mu.size(), &tmp2[0]);
      for(indtype i = 0, iend = tmp.size(); i < iend; ++i)
      {
        Rcpp::Rcout << tmp[i] << " ";
      }
    }
    else
    {
      for(indtype i = 0, iend = cholU.size(); i < iend; ++i)
      {
        Rcpp::Rcout << cholU[i] <<" ";
      }
    }
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "density container: ";
    for(indtype i = 0, iend = ptr.size(); i < iend; ++i)
    {
      Rcpp::Rcout << ptr[i] << " ";
    }
    Rcpp::Rcout << "\n\n";
  }


  // M is a tempoarary vector of size d, the dimensionality of a data point
  // inline valtype densityEval(valtype *x, indtype d, valtype *M, indtype withAlphaWithPi = 3)
  inline valtype densityEval(valtype *x, indtype d, valtype *M, valtype piConstMul)
  {
    if(sqrtOfDet <= 0) return 0;


    valtype *mean = &*mu.begin();
    valtype *cl = &*cholU.begin();


    // mahalanobis distance computation
    valtype *Mst = M;
    valtype *Mbegin = Mst;
    valtype *Mend = M + d;


    *Mst = (*x - *mean) / *cl;
    valtype mahanaD = *Mst * *Mst;
    indtype rowLen = 1;
    for(;;)
    {
      ++Mst;
      if(Mst >= Mend) break;
      ++x;
      ++mean;
      cl += rowLen;
      ++rowLen;
      valtype numerator = *x - *mean - std::inner_product(Mbegin, Mst, cl, 0.0);
      valtype denominator = *(cl + (Mst - Mbegin));
      if(denominator == 0.0) *Mst = std::numeric_limits<valtype>::max();
      // if(numerator + denominator == 0.0) *Mst = 1.0;
      else *Mst = numerator / denominator;
      // *Mst = (*x - *mean - std::inner_product(Mbegin, Mst, cl, 0.0)) / *(cl + (Mst - Mbegin));
      // if(!std::isfinite(*Mst))
      // {
      //   Rcpp::Rcout << "*Mst not finite" << ", ";
      //   Rcpp::Rcout << "fenzi = " << *x - *mean - std::inner_product(Mbegin, Mst, cl, 0.0) << ", ";
      //   Rcpp::Rcout << "fenmu = " << *(cl + (Mst - Mbegin)) << "\n";
      //   Rcpp::Rcout << "sqrtOfDet = " << sqrtOfDet << "\n";
      // }
      mahanaD += *Mst * *Mst;
    }


    valtype rst = std::exp(-mahanaD / 2) / sqrtOfDet;
    // if(!std::isfinite(mahanaD)) Rcpp::Rcout << "rst = " << rst << "\n";
    // if(!std::isfinite(rst))
    // {
    //   Rcpp::Rcout << "-mahanaD / 2 = " << -mahanaD / 2 << ", ";
    //   Rcpp::Rcout << "std::exp(-mahanaD / 2) = " << std::exp(-mahanaD / 2) << ", ";
    //   Rcpp::Rcout << "sqrtOfDet = " << sqrtOfDet << "\n";
    // };
    return rst * alpha * piConstMul;
  }
};











