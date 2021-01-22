// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

#include <cmath>

using namespace Rcpp;
using namespace RcppParallel;

struct OneColTask : public Worker
{
  const RVector<double> mySeq1;
  const double myMean1;
  const double mySd1;
  const double lengthX1;
  const int n_sim1;
  
  RMatrix<double> myRetMat1;
  
  OneColTask(const NumericVector mySeq, const double myMean, const double mySd, const double lengthX, const int n_sim, NumericMatrix myRetMat) 
  : mySeq1(mySeq), myMean1(myMean), mySd1(mySd), lengthX1(lengthX), n_sim1(n_sim), myRetMat1(myRetMat) {}
  
 
  void operator()(std::size_t begin, std::size_t end) {
            
    std::vector<double> mySims (lengthX1);
    
    for (std::size_t j = begin; j < end; j++) {
           
      for(std::size_t i = 0; i < lengthX1; i++){
        mySims[i] = R::rnorm(myMean1, mySd1);
      }

      std::sort(mySims.begin(), mySims.end());

       for(std::size_t k = 0; k < mySeq1.length(); k++){ 
          double weights = (mySeq1[k] * lengthX1) - floor(mySeq1[k] * lengthX1);

          myRetMat1.column(j)[k] = weights * mySims[(int)std::max(floor(mySeq1[k] * lengthX1),0.0)] + (1-weights) * mySims[(int)std::min(ceil(mySeq1[k] * lengthX1),lengthX1-1)];

       }
     }
    
  }
};


// [[Rcpp::export]]
NumericMatrix myQQNormIntern(NumericVector mySeq, double myMean, double mySd, double lengthX, int n_sim) {
  
  NumericMatrix myReturnMatrix(mySeq.length(), n_sim);
    
  OneColTask oneColTask(mySeq,myMean,mySd,lengthX,n_sim,myReturnMatrix);
  
  parallelFor(0, myReturnMatrix.ncol(), oneColTask);
  
  return myReturnMatrix;
}
