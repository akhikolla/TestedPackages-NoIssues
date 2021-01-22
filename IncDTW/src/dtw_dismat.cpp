
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

struct wdm_ws_ea : public Worker {

   const std::vector<arma::vec> vov;
   std::vector<int> ii;
   std::vector<int> jj;
   RVector<double> output;
   bool normalize;
   std::string step_pattern;
   int ws;
   double threshold;

   // initialize from Rcpp input and output matrixes (the RMatrix class
   // can be automatically converted to from the Rcpp matrix type)
   wdm_ws_ea(std::vector<arma::vec> vov, std::vector<int> ii, std::vector<int> jj,
             Rcpp::NumericVector output, 
             bool normalize,  
             std::string step_pattern, 
             int ws, double threshold)
      : vov(vov), ii(ii), jj(jj), output(output), normalize(normalize),
        step_pattern(step_pattern), ws(ws), threshold(threshold){}

   // Function call operator to iterate over a specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {

      double norm = 1;
      double n, m;
      int i,j;

      for (std::size_t k = begin; k < end; ++k) {
         i = ii[k];
         j = jj[k];
         if(normalize){
            n = vov.at(i).n_rows;
            m = vov.at(j).n_rows;
            norm = 1/(n+m);
         }
         output[k] = multp_dtw2vec_ws_ea(vov.at(i) , vov.at(j), step_pattern, ws, threshold) * norm;
      }
   }
};


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


struct wdm_mv_ws_ea : public Worker {
   
   const std::vector<arma::mat> vom;
   std::vector<int> ii;
   std::vector<int> jj;
   RVector<double> output;
   bool normalize;
   std::string step_pattern;
   std::string dist_method;
 
   int ws;
   double threshold;
   
   // initialize from Rcpp input and output matrixes (the RMatrix class
   // can be automatically converted to from the Rcpp matrix type)
   wdm_mv_ws_ea(std::vector<arma::mat> vom, std::vector<int> ii, std::vector<int> jj,
                      Rcpp::NumericVector output, 
                      bool normalize, std::string step_pattern, 
                      std::string dist_method, int ws, double threshold)
      : vom(vom), ii(ii), jj(jj), output(output), normalize(normalize), 
        step_pattern(step_pattern), dist_method(dist_method), ws(ws), threshold(threshold){}
   
   // Function call operator to iterate over a specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {
      
      double norm = 1;
      double n, m;
      int i,j;
      
      for (std::size_t k = begin; k < end; ++k) {
         i = ii[k];
         j = jj[k];
         if(normalize){
            n = vom.at(i).n_rows;
            m = vom.at(j).n_rows;
            norm = 1/(n+m);
         }
         output[k] = multp_dtw2vec_mv_ws_ea(vom.at(i) , vom.at(j), step_pattern, 
                                          dist_method, ws, threshold) * norm;
      }
   }
};


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
Rcpp::NumericVector parallel_dm_dtw(Rcpp::List lot,
                                    std::vector<int> ii,
                                    std::vector<int> jj,
                                    bool normalize, std::string step_pattern,
                                    int ws, double threshold ){
 
   // allocate the vector/matrix we will return
   int N = lot.length();
   int NN = N * (N-1)/2;
   Rcpp::NumericVector output(NN);


   // Convert lot -> vom/vov, list of time series -> vector of matrices/ vector of vectors
   std::vector<arma::vec> vov;
   for (List::iterator it = lot.begin(); it != lot.end(); ++it) {
      vov.push_back(as<arma::vec >(*it));
   }

   // create the worker
   wdm_ws_ea wdm(vov, ii, jj, output, normalize, step_pattern, ws, threshold);
   parallelFor(0, NN, wdm);
   
   return output;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
Rcpp::NumericVector parallel_dm_dtw_mv(Rcpp::List lot,
                                       std::vector<int> ii,
                                       std::vector<int> jj,
                                       bool normalize, 
                                       std::string step_pattern,
                                       std::string dist_method,
                                       int ws, double threshold ){
   
   // allocate the vector/matrix we will return
   int N = lot.length();
   int NN = N * (N-1)/2;
   Rcpp::NumericVector output(NN);
   
   
   // Convert lot -> vom/vov, list of time series -> vector of matrices/ vector of vectors
   std::vector<arma::mat> vom;
   for (List::iterator it = lot.begin(); it != lot.end(); ++it) {
      vom.push_back(as<arma::mat >(*it));
   }
   
   // create the worker
   wdm_mv_ws_ea wdm(vom, ii, jj, output, normalize, step_pattern, dist_method, ws, threshold);
   parallelFor(0, NN, wdm);
   

   return output;
}


