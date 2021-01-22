
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "IncDTW.h"

#include <cmath>
#include <list>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace Rcpp;
using namespace RcppParallel;
using namespace std;



// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


struct wdv_dtw_par : public Worker {

   // Input matrix to read from must be of the RMatrix<T> form
   // if using Rcpp objects
   // const RMatrix<double> mat;
   const arma::vec Q;
   const std::vector<arma::vec> vov;


   // Output vector to write to must be of the RVector<T> form
   // if using Rcpp objects
   RVector<double> output;

   bool normalize;
   
   // const Rcpp::List params;
   std::string step_pattern;
   int ws;
   double threshold;
   


   // initialize from Rcpp input and output matrixes (the RMatrix class
   // can be automatically converted to from the Rcpp matrix type)
   wdv_dtw_par(arma::vec Q, 
              std::vector<arma::vec> vov, 
              Rcpp::NumericVector output,
              bool normalize,
              std::string step_pattern,
              int ws, 
              double threshold)
      : Q(Q), vov(vov), output(output), normalize(normalize), 
        step_pattern(step_pattern), ws(ws), threshold(threshold){}

   // Function call operator to iterate over a specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {

      double norm = 1;
      double n = Q.n_rows;
      double m;

      for (std::size_t j = begin; j < end; ++j) {
         if(normalize){
            m = vov.at(j).n_rows;
            norm = 1/(n+m);
         }
         // output[j] = fdtw(Q , vov.at(j), params ) * norm;
         output[j] = multp_dtw2vec_ws_ea(Q , vov.at(j), step_pattern, ws , threshold) * norm;
      }
   }
};


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


struct wdv_dtw_par_mv : public Worker {
   
   // Input matrix to read from must be of the RMatrix<T> form
   // if using Rcpp objects
   // const RMatrix<double> mat;
   const arma::mat Q;
   const std::vector<arma::mat> vom;
   
   
   // Output vector to write to must be of the RVector<T> form
   // if using Rcpp objects
   RVector<double> output;
   
   bool normalize;
   std::string step_pattern;
   std::string dist_method;
   int ws;
   double threshold;
   
   
   
   // initialize from Rcpp input and output matrixes (the RMatrix class
   // can be automatically converted to from the Rcpp matrix type)
   wdv_dtw_par_mv(arma::mat Q, 
                  std::vector<arma::mat> vom, 
                  Rcpp::NumericVector output,
                  bool normalize, 
                  std::string step_pattern, 
                  std::string dist_method, 
                  int ws, double threshold)
      : Q(Q), vom(vom), output(output), normalize(normalize), step_pattern(step_pattern),
        dist_method(dist_method), ws(ws), threshold(threshold){}
   
   // Function call operator to iterate over a specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {
      
      double norm = 1;
      double n = Q.n_rows;
      double m;
      
      for (std::size_t j = begin; j < end; ++j) {
         if(normalize){
            m = vom.at(j).n_rows;
            norm = 1/(n+m);
         }
         // output[j] = fdtw(Q , vom.at(j), params ) * norm;
         output[j] = multp_dtw2vec_mv_ws_ea(Q , vom.at(j), step_pattern, dist_method, ws, threshold ) * norm;
      }
   }
};


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
Rcpp::NumericVector parallel_dv_dtw(arma::vec Q, Rcpp::List lot, bool normalize, 
                                    std::string step_pattern, 
                                    int ws, double threshold){

   // allocate the vector/matrix we will return
   Rcpp::NumericVector output(lot.length());


   // Convert lot -> vom/vov, list of time series -> vector of matrices/ vector of vectors
    std::vector<arma::vec> vov;
      for (List::iterator it = lot.begin(); it != lot.end(); ++it) {
         vov.push_back(as<arma::vec >(*it));
      }

      // create the worker
      wdv_dtw_par wdv(Q, vov, output, normalize, step_pattern, ws, threshold);
      // call it with parallelFor
      parallelFor(0, lot.length(), wdv);

   return output;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
Rcpp::NumericVector parallel_dv_dtw_mv(arma::mat Q, Rcpp::List lot, bool normalize, 
                                       std::string step_pattern, std::string dist_method, 
                                       int ws, double threshold)
{
   // allocate the vector/matrix we will return
   Rcpp::NumericVector output(lot.length());
   
   
   // Convert lot -> vom/vov, list of time series -> vector of matrices/ vector of vectors
   std::vector<arma::mat> vom;
   for (List::iterator it = lot.begin(); it != lot.end(); ++it) {
      vom.push_back(as<arma::mat >(*it));
   }
   
   // create the worker
   wdv_dtw_par_mv wdv(Q, vom, output, normalize, step_pattern, dist_method, ws, threshold);
   // call it with parallelFor
   parallelFor(0, lot.length(), wdv);
   
   return output;
}



