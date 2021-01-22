// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
# include <RcppArmadillo.h>
# include <RcppParallel.h>
# include "h/dnyTasking.hpp"
# include "h/macros.hpp"
# include "h/G.hpp"
# include "h/funs.hpp"
using namespace RcppParallel;
using namespace Rcpp;


// [[Rcpp::export]]
double testGdensity(
    NumericVector x, NumericVector mu, NumericVector sigma, double alpha = 1)
{
  G<int, double> gaussian;
  gaussian.alpha = alpha;
  gaussian.mu.assign(mu.begin(), mu.end());
  int d = mu.size();
  // Rcout << "d = " << d << "\n";
  vec<double> triSigma(std::size_t(d) * (d + 1) / 2);
  fullSigmaToTriSigma(&*sigma.begin(), &triSigma.front(), d);
  if(false)
  {
    for(int i = 0, iend = triSigma.size(); i < iend; ++i)
      Rcout << triSigma[i] << ", ";
    Rcout << "\n\n";
  }
  bool b = gaussian.computeCholUandSqrtOfDet(triSigma);
  if(!b)
  {
    Rcout << "problematic\n";
    return 0;
  }
  vec<double> tmp(mu.size());
  double pi_ = std::pow(2.0 * M_PI, d * (-0.5));
  return gaussian.densityEval(&x[0], mu.size(), &tmp[0], pi_);
}


/*
template<typename indtype, typename valtype>
void update1G(G<indtype, valtype> &gaussian, valtype &diffS, // valtype eps,
              indtype d, indtype Xsize, valtype *X,
              const valtype *pointWeight, valtype *rowSum, valtype *errorS,
              valtype *buffer, bool updateAlpha = true)
  // buffer is a temporary container of size Xsize + d * (d + 1) / 2 + d
{
  indtype triMatSize = std::size_t(d) * (d + 1) / 2;
  valtype *&W = buffer; // The first Xsize elements of buffer stores weights
  valtype sumW = 0;
  for(indtype i = 0; i < Xsize; ++i)
  {
    W[i] = 0;
    // if(rowSum[i] >= invInf) W[i] = gaussian.ptr[i] / rowSum[i];
    if(rowSum[i] > 0) W[i] = gaussian.ptr[i] / rowSum[i];
    W[i] *= pointWeight[i]; // branch prediction will optimize this away
    sumW += W[i];
  }


  valtype alphaNew = gaussian.alpha;
  if(updateAlpha) alphaNew = sumW / Xsize;


  diffS = 0;
  diffS += relaErr<valtype> (alphaNew, gaussian.alpha);
  gaussian.alpha = alphaNew;


  for(indtype i = 0; i < Xsize; ++i) W[i] /= sumW;


  valtype *oldmu = buffer + Xsize; // buffer[Xsize] -> buffer[Xsize + d] stores the old mean
  std::copy(gaussian.mu.begin(), gaussian.mu.end(), oldmu);


  // Update mean
  valtype *mean = &gaussian.mu.front();
  std::fill(mean, mean + d, 0);
  for(indtype i = 0; i < Xsize; ++i)
  {
    valtype *val= X + i * d;
    for(indtype j = 0; j < d; ++j)
    {
      mean[j] += W[i] * val[j];
    }
  }


  // relative difference
  diffS += relaErrSum(mean, oldmu, d);


  valtype *Sigma = buffer + Xsize; // buffer[Xsize] -> buffer[Xsize + triMatSize]
  // stores the new covariance (upper tri)
  std::fill(Sigma, Sigma + triMatSize, 0);


  // update sigma
  valtype *x_mu = Sigma + triMatSize; // buffer[Xsize + triMatSize] stores x - mu, temporary
  for(indtype i = 0; i < Xsize; ++i)
  {
    valtype *val = X + i * d;
    for(indtype t = 0; t < d; ++t)
    {
      x_mu[t] = val[t] - mean[t];
    }


    valtype *s = Sigma;
    for(indtype p = 0; p < d; ++p)
    {
      for(indtype q = 0; q <= p; ++q) // upper triangle
      {
        *s += x_mu[p] * x_mu[q] * W[i];
        ++s;
      }
    }
  }


  // bool b = gaussian.computeCholUandSqrtOfDet(Sigma);
  bool b = gaussian.computeCholUandSqrtOfDet(Sigma);
  // After this function, Sigma contains the old covariance matrix
  if(!b)
  {
    gaussian.alpha = 0;
    diffS += triMatSize;
    return;
  }


  // relative difference
  valtype *&cholPrevous = Sigma; // Now Sigma stores the old covariance matrix (upper tri)
  diffS += relaErrSum(&gaussian.cholU[0], cholPrevous, triMatSize);
}
*/


template<typename indtype, typename valtype>
void update1G(G<indtype, valtype> &gaussian, valtype &diffS, // valtype eps,
              indtype d, indtype Xsize,
              valtype *X, const valtype *pointWeight,
              valtype *rowSum, valtype *errorS, valtype *buffer)
  // buffer is a temporary container of size Xsize + d * (d + 1) / 2 + d
{
  indtype triMatSize = std::size_t(d) * (d + 1) / 2;
  valtype *&W = buffer; // The first Xsize elements of buffer stores weights
  valtype sumW = 0;
  for(indtype i = 0; i < Xsize; ++i)
  {
    W[i] = 0;
    if(rowSum[i] > 0) W[i] = gaussian.ptr[i] / rowSum[i];
    W[i] *= pointWeight[i]; // Branch prediction will optimize this away
    sumW += W[i];
  }


  diffS = 0;
  // std::cout << "gaussian.updateAlpha = " + std::to_string(gaussian.updateAlpha) + ", ";
  if(gaussian.updateAlpha)
  {
    valtype alphaNew = sumW / Xsize;
    // Rcout << "alphaNew = " << alphaNew << ", ";
    diffS += relaErr<valtype> (alphaNew, gaussian.alpha);
    gaussian.alpha = alphaNew;
  }


  if(gaussian.updateMean or gaussian.updateSigma)
  {
    for(indtype i = 0; i < Xsize; ++i) W[i] /= sumW;
  }


  // Update mean
  valtype *mean = &gaussian.mu.front();
  if(gaussian.updateMean)
  {
    valtype *oldmu = buffer + Xsize; // buffer[Xsize] -> buffer[Xsize + d] stores the old mean
    std::copy(gaussian.mu.begin(), gaussian.mu.end(), oldmu);
    std::fill(mean, mean + d, 0);
    for(indtype i = 0; i < Xsize; ++i)
    {
      valtype *val= X + i * d;
      for(indtype j = 0; j < d; ++j)
      {
        mean[j] += W[i] * val[j];
      }
    }
    // Relative difference
    diffS += relaErrSum(mean, oldmu, d);
  }


  if(gaussian.updateSigma)
  {
    valtype *Sigma = buffer + Xsize; // buffer[Xsize] -> buffer[Xsize + triMatSize]
    // Stores the new covariance (upper tri)
    std::fill(Sigma, Sigma + triMatSize, 0);
    // Update sigma
    valtype *x_mu = Sigma + triMatSize; // buffer[Xsize + triMatSize] stores x - mu, temporary
    for(indtype i = 0; i < Xsize; ++i)
    {
      valtype *val = X + i * d;
      for(indtype t = 0; t < d; ++t) x_mu[t] = val[t] - mean[t];
      valtype *s = Sigma;
      for(indtype p = 0; p < d; ++p)
      {
        for(indtype q = 0; q <= p; ++q) // upper triangle
        {
          *s += x_mu[p] * x_mu[q] * W[i];
          ++s;
        }
      }
    }
    bool b = gaussian.computeCholUandSqrtOfDet(Sigma);
    // After this function, Sigma contains the old covariance matrix
    if(!b)
    {
      gaussian.alpha = 0;
      diffS += triMatSize;
      return;
    }
    // Relative difference
    valtype *&cholPrevous = Sigma; // Now Sigma stores the old covariance matrix (upper tri)
    diffS += relaErrSum(&gaussian.cholU[0], cholPrevous, triMatSize);
  }


}






template<typename indtype, typename valtype>
struct updateParaConventional: public Worker
{
  indtype d, Xsize, gmodelSize;
  // valtype eps;
  valtype *X, *pointWeight, *rowSum;
  G<indtype, valtype> *gmodel;
  valtype *errorS;
  vec<valtype> *tmpCntr;
  bool updateAlpha;
  dynamicTasking *dT;


  void operator()(std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI)) break;
      update1G(gmodel[objI], errorS[objI], // eps,
               d, Xsize, X, pointWeight,
               rowSum, errorS, &tmpCntr[st][0]);
    }
  }


  updateParaConventional(
    indtype d,
    indtype Xsize,
    indtype gmodelSize,
    valtype *X,
    valtype *pointWeight,
    valtype *rowSum,
    G<indtype, valtype> *gmodel, valtype *errorS, indtype NofCPU, bool updateAlpha):
    d(d), Xsize(Xsize), gmodelSize(gmodelSize), // eps(eps),
    X(X), pointWeight(pointWeight), rowSum(rowSum), gmodel(gmodel), errorS(errorS),
    updateAlpha(updateAlpha)
  {
    dynamicTasking dt(NofCPU, gmodelSize); dT = &dt;
    vec<vec<valtype> > tmp(NofCPU, vec<valtype> (
        std::size_t(d) * (d + 1) / 2 + Xsize + d));
    tmpCntr = &tmp[0];
    parallelFor(0, NofCPU, *this);
  }
};




// X[, i] is the i_th data point
// [[Rcpp::export]]
NumericMatrix findSpreadedMeanWrapper(NumericMatrix X, int K, int maxCore = 7)
{
  return findSpreadedMean(X, K, maxCore);
}




// [[Rcpp::export]]
NumericMatrix makeCovariancesWrapper(NumericMatrix X, int K)
{
  return makeCovariances01(X, K);
}




// dat[, i] is the i_th data point
// mu[, i] is the center of the i_th Gaussian componenet
// sigma[, i] is the full covariance matrix of the i_th Gaussian component
List paraGmmFullinit(
    NumericMatrix dat,
    NumericVector alpha,
    NumericMatrix mu,
    NumericMatrix sigma,
    NumericVector xweight,
    int Nthreads,
    int maxit,
    double maxEigenRatio,
    double eps,
    double annihilationEPS,
    double tlimit,
    int verbose,
    LogicalVector updateAlpha,
    LogicalVector updateMean,
    LogicalVector updateSigma,
    bool paraConvergeMaxErr,
    bool loglikehoodConverge,
    int loglikehoodConvergeBlock)
{
  int d = dat.nrow();
  int Xsize = dat.ncol();


  double *X = &*dat.begin();


  // gaussian kernel initialization
  vec<G<int, double> > Gvec(alpha.size());
  {
    double *muptr = &*mu.begin();
    double *sigmaptr = &*sigma.begin();
    for(int i = 0, iend = Gvec.size(); i < iend; ++i, muptr += d, sigmaptr += d * d)
    {
      Gvec[i].updateAlpha = updateAlpha[i];
      // Rcout << "initialization: Gvec[i].updateAlpha = " << Gvec[i].updateAlpha << "\n";
      Gvec[i].updateMean = updateMean[i];
      Gvec[i].updateSigma = updateSigma[i];


      Gvec[i].alpha = alpha[i];
      Gvec[i].mu.assign(muptr, muptr + d);
      // Extract tri-sigma, perform cholesky decomposition
      {
        Gvec[i].cholU.resize(std::size_t(d) * (d + 1) / 2);
        fullSigmaToTriSigma(sigmaptr, &Gvec[i].cholU[0], d);
        Gvec[i].computeCholUandSqrtOfDet(Gvec[i].cholU);
      }
    }


    // Annihilate if alpha < annihilationEPS
    annihilateGinVec(Gvec, annihilationEPS);


    // Density container initialization
    for(int i = 0, iend = Gvec.size(); i < iend; ++i)
      Gvec[i].ptr.resize(Xsize);
  }


  // EM
  NumericVector rowSum(Xsize);
  {
    double *xw = &xweight[0];
    int gaussianParaN = std::size_t(d + 1) * d / 2 + d + 1;
    double endTime = std::clock() + Nthreads * tlimit * CLOCKS_PER_SEC;
    vec<double> auxCntr;
    vec<double> errorS;
    double previousLL = -1e300;
    int NllunderEps = 0; // The last `NllunderEps` iterations resulted in loglikelihood less than `eps`.
    for(int iter = 0, iterEnd = maxit; iter < iterEnd; ++iter)
    {
      int GvecSize = Gvec.size();
      cmptDensity<int, double> (d, Xsize, GvecSize, X, &Gvec[0], Nthreads);
      cmptRowSum<int, double> (Xsize, GvecSize, &Gvec.front(), &rowSum[0], auxCntr, Nthreads);
      errorS.assign(GvecSize, 0);


      int totalParameters = GvecSize * gaussianParaN;
      updateParaConventional<int, double> (
          d, Xsize, GvecSize, // eps,
          X, xw, &rowSum[0], &Gvec.front(), &errorS.front(), Nthreads, updateAlpha);


      annihilateGinVec(Gvec, annihilationEPS); // Erase components that have small weights.
      cleanGaussianKernelNotMeetingEigenRatio(Gvec, d, maxEigenRatio, Nthreads);
      // The function does nothing if maxEigenRatio <= 0.
      // Erase components that have outrageous eigen ratios.


      int totalParametersAfterCleansing = Gvec.size() * gaussianParaN;


      double avgRelativeErr = 0;
      if(totalParametersAfterCleansing == totalParameters)
      {
        if(!loglikehoodConverge)
        {
          if(paraConvergeMaxErr) avgRelativeErr = *std::max_element(errorS.begin(), errorS.end());
          else avgRelativeErr = std::accumulate(errorS.begin(), errorS.end(), 0.0) / totalParameters;
        }
        else
        {
          double ll = 0;
          for(int i = 0, iend = rowSum.size(); i < iend; ++i) ll += std::log(rowSum[i]);
          if(ll < (((previousLL > 0) * 2 - 1) * eps + 1) * previousLL) ++NllunderEps;
          else NllunderEps = 0;
          previousLL = ll;
        }
      }
      else
      {
        if(verbose)
          Rcout << "Eigenvalue ratios exceed or component weights fall below thresholds, " <<
            "Gaussian components annihilated.\n";
        avgRelativeErr = eps * 2;
        NllunderEps = 0;
      }


      if(!loglikehoodConverge)
      {
        if(avgRelativeErr < eps or double(std::clock()) > endTime) break;
      }
      else
      {
        if(NllunderEps >= loglikehoodConvergeBlock or double(std::clock()) > endTime) break;
      }


      // if(verbose != 0 and (iter + 1) % verbose == 0)
      if(verbose != 0)
      {
        if(!loglikehoodConverge)
        {
          std::string meanMax = "";
          if(paraConvergeMaxErr) meanMax += ", max absolute relative change in parameters = ";
          else meanMax += ", mean absolute relative change in parameters = ";
          Rcout << "iteration = " << iter + 1 << meanMax << avgRelativeErr <<
              ", number of remaining kernels = " << Gvec.size() << "\n";
        }
        else
        {
          Rcout << "iteration = " << iter + 1 << ", log-likelihood = " << previousLL <<
            ", each of the last " << NllunderEps << " iterations yielded less than " <<
              eps * 100 << "% increase in log-likelihood\n";
        }
      }
    }
  }


  // out
  if(true)
  {
    NumericVector w(Gvec.size());
    NumericMatrix mu(d, Gvec.size()), sigma(d * d, Gvec.size());
    int imu = 0, isigma = 0;
    vec<double*> tmp;
    for(int i = 0, iend = Gvec.size(); i < iend; ++i, imu += d, isigma += d * d)
    {
      w[i] = Gvec[i].alpha;
      std::copy(Gvec[i].mu.begin(), Gvec[i].mu.end(), &mu[imu]);
      tmp.resize(d);
      triCholToFullSigma(&Gvec[i].cholU.front(), &sigma[isigma], d, &tmp[0]);
    }


    IntegerVector clust(Xsize);
    for(int i = 0; i < Xsize; ++i)
    {
      int whichMax = 0;
      double denMax = Gvec[0].ptr[i];
      for(int j = 1, jend = Gvec.size(); j < jend; ++j)
      {
        if(denMax >= Gvec[j].ptr[i]) continue;
        whichMax = j;
        denMax = Gvec[j].ptr[i];
      }
      clust[i] = whichMax;
    }


    // if(!updateAlpha and mu.size() < alpha.size())
    // {
    //   std::fill(w.begin(), w.end(), 1.0 / mu.size());
    //   double r = (mu.size() + 0.0) / alpha.size();
    //   for(int i = 0, iend = rowSum.size(); i < iend; ++i) rowSum[i] *= r;
    // }
    // Compute rowSum and w.
    if(true)
    {
      for(int i = 0, iend = w.size(); i < iend; ++i) w[i] = Gvec[i].alpha;
      // for(int i = 0, iend = w.size(); i < iend; ++i) Rcout << Gvec[i].alpha << ", ";
      double r = 1.0 / std::accumulate(w.begin(), w.end(), 0.0);
      for(int i = 0, iend = w.size(); i < iend; ++i) w[i] *= r;
      for(int i = 0, iend = rowSum.size(); i < iend; ++i) rowSum[i] *= r;
    }


    return List::create(Named("alpha") = w, Named("mu") = mu,
                        Named("sigma") = sigma, Named("fitted") = rowSum,
                        Named("clusterMember") = clust);
  }
}




// xweight does not need to sum up to 1
// [[Rcpp::export]]
List paraGmm(
    NumericMatrix X,
    NumericVector Xw,
    int G,
    NumericVector alpha,
    NumericMatrix mu,
    NumericMatrix sigma,
    double eigenRatioLim,
    double convergenceEPS,
    double alphaEPS,
    int maxIter,
    double tlimit,
    int verbose,
    int maxCore,
    LogicalVector updateAlpha,
    LogicalVector updateMean,
    LogicalVector updateSigma,
    bool paraConvergeMaxErr,
    bool loglikehoodConverge,
    int loglikehoodConvergeBlock
  )
{
  if(alpha.size() == 0)
  {
    alpha = NumericVector(G, 1.0 / G);
  }
  if(mu.size() == 0)
  {
    mu = findSpreadedMean(X, G, maxCore);
  }
  if(sigma.size() == 0)
  {
    sigma = makeCovariances01(X, G);
  }


  // Validate initial covariance matrices and maxEigenRatio.
  if(true)
  {
    int d = mu.nrow();
    arma::mat tmp(d, d);
    arma::colvec tmpv(d);
    int dd = d * d;
    for(int i = 0, iend = mu.ncol(); i < iend; ++i)
    {
      std::copy(&sigma[0] + i * dd, &sigma[0] + i * dd + dd, &tmp[0]);
      arma::eig_sym(tmpv, tmp);
      if(tmpv[0] <= 0)
      {
        Rcout << "The " << i << "th covariance matrix is not positive-definite. Quit.";
        return List::create();
      }
      double tmpRatio = tmpv[d - 1] / tmpv[0];
      if(eigenRatioLim > 0 and tmpRatio > eigenRatioLim)
      {
        Rcout << "The " << i << "th covariance matrix's max:min eigen ratio exceeds threshold. Quit.";
        return List::create();
      }
    }
  }




  LogicalVector alltrue(alpha.size()), allfalse(alpha.size());
  for(int i = 0, iend = alltrue.size(); i < iend; ++i)
  {
    alltrue[i] = true;
    allfalse[i] = false;
  }


  vec<unsigned char> updateAlpha_, updateMean_, updateSigma_;
  if(updateAlpha.size() == 1)
  {
    if(updateAlpha[0]) updateAlpha = alltrue;
    else updateAlpha = allfalse;
  }


  if(updateMean.size() == 1)
  {
    if(updateMean[0]) updateMean = alltrue;
    else updateMean = allfalse;
  }


  if(updateSigma.size() == 1)
  {
    if(updateSigma[0]) updateSigma = alltrue;
    else updateSigma = allfalse;
  }


  return paraGmmFullinit(
    X, alpha, mu, sigma, Xw, maxCore, maxIter,
    eigenRatioLim, convergenceEPS, alphaEPS, tlimit, verbose,
    updateAlpha, updateMean, updateSigma,
    paraConvergeMaxErr, loglikehoodConverge, loglikehoodConvergeBlock);
}
































