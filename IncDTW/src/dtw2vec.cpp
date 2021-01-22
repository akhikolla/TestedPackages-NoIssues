#include <RcppArmadillo.h>
#include "IncDTW.h"
using namespace Rcpp;
using namespace std;
//using namespace stdlib; 
// [[Rcpp::plugins(cpp11)]]


XPtr<funcPtr_dist> select_dist(std::string dist_method) {
   
   if (dist_method == "norm1"){
      return(XPtr<funcPtr_dist>(new funcPtr_dist(&dist_norm1)));
      
   }else if (dist_method == "norm2_square"){
      return(XPtr<funcPtr_dist>(new funcPtr_dist(&dist_norm2_square)));
      
   }else if (dist_method == "norm2"){
      return(XPtr<funcPtr_dist>(new funcPtr_dist(&dist_norm2)));
      
   }else{
      return XPtr<funcPtr_dist>(R_NilValue); // runtime error as NULL no XPtr
   }
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


typedef double (*funcPtr_step_vec)(const double gcm10, 
                                   const double gcm11, 
                                   const double gcm01, 
                                   const double cm00);

XPtr<funcPtr_step_vec> selectVecStep(std::string step_pattern) {
   if (step_pattern == "symmetric1")
      return(XPtr<funcPtr_step_vec>(new funcPtr_step_vec(&mystep_symmetric1)));
   else if (step_pattern == "symmetric2")
      return(XPtr<funcPtr_step_vec>(new funcPtr_step_vec(&mystep_symmetric2)));
   else
      return XPtr<funcPtr_step_vec>(R_NilValue); // runtime error as NULL no XPtr
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


typedef std::function<double( const arma::mat&, const arma::mat&, int, int, int )> distFunction;


distFunction selectDistFunction(std::string dist_method) 
{
   distFunction f;
   if(dist_method == "norm1"){
      f = std::bind(&dist_norm1 , 
                    std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
   }else if(dist_method == "norm2"){
      f = std::bind(&dist_norm2 , 
                    std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
   }else if(dist_method == "norm2_square"){
      f = std::bind(&dist_norm2_square , 
                    std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
   }
   return f;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


typedef std::function<double( double, double, double, double )> vecStepFunction;


vecStepFunction selectVecStepFunction(std::string step_pattern)
{
   vecStepFunction f;
   if(step_pattern == "symmetric1"){
      f = std::bind(&mystep_symmetric1 , std::placeholders::_1, std::placeholders::_2,
                                         std::placeholders::_3, std::placeholders::_4);
   }else if(step_pattern == "symmetric2"){
      f = std::bind(&mystep_symmetric2 , std::placeholders::_1, std::placeholders::_2,
                                         std::placeholders::_3, std::placeholders::_4);
   }
   return f;
}




// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec_cm (Rcpp::NumericMatrix cm, std::string step_pattern)
{
   
   int nx = cm.nrow();
   int ny = cm.ncol();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double ret;
   
   // first column
   *p1 = cm(0,0);
   for(int i=1; i<nx; i++){
      p1[i] = cm(i,0) + p1[i-1];
   }

   // set step function
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;

  for(int j=1; j < ny; j++){
      *p2 = cm(0,j) + *(p1);
      
      for(int i=1; i<nx; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i, j));
      }
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_dtw2vec_cm_inc (Rcpp::NumericVector gcm_lc, Rcpp::NumericMatrix cm, 
                         std::string step_pattern)
{
   // x ... time series that is fixed
   // y ... time series of new observations, ONLY new observations to be appended
   // gcm_lc ... last column of old GCM
   
   int n = cm.nrow();
   int m = cm.ncol();// = nnewObs
   
   double * p1 = new double [n];
   double * p2 = new double [n];
   double * ptmp; 
   
   double mynan;
   NumericVector gcm_lr_new(m);
   NumericVector gcm_lc_new(n);
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   // set step sunction
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   
   if(n != gcm_lc.size()){
      return mynan;
   }
   
   // first column
   for(int i=0; i<n; i++){
      p1[i] = gcm_lc[i];
   }
   
   for(int j=0; j < m; j++){
      *p2 = cm(0,j) + *(p1);
      
      for(int i=1; i<n; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,j));
      }
      gcm_lr_new[j] = *(p2+n-1);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   for(int i=0; i<n; i++){
      gcm_lc_new[i] = *(p1+i);
   }
   
   
   List ret;
   ret["gcm_lr_new"] = gcm_lr_new;
   ret["gcm_lc_new"] = gcm_lc_new;
   ret["distance"] = *(p1+n-1);
   delete[] p1;
   delete[] p2;
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// [[Rcpp::export]]
double cpp_dtw2vec_cm_ws_ea (Rcpp::NumericMatrix cm, 
                          std::string step_pattern, int ws, double threshold)
{
   
   int iBegin = 0;
   int iEnd = 0;
   
   int nx = cm.nrow();
   int ny = cm.ncol();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double mynan;
   double ret;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   // set step sunction
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   
   //initialize a and b with NAN
   for(int i=0; i < nx; i++){
      p1[i] = mynan;
      p2[i] = mynan;
   }
   
   // first column
   p1[0] = cm(0,0);
   if(p1[0] > threshold)  return(mynan); 
   
   iEnd   = std::min(nx, ws+1);
   for(int i=1; i < iEnd; i++){
      p1[i] = cm(i,0) + p1[i-1];
      if(p1[i] > threshold) p1[i] = mynan;
   }
   
   
   for(int j=1; j < ny; j++){
      nanCounter = 0;
      iBegin = j-ws;
      if(iBegin <= 0){
         *p2 = cm(0,j) + *(p1);
         if(*(p2) > threshold){
            *(p2) = mynan;
            nanCounter ++;
         }
         iBegin = 1;   
      }else if (iBegin == 1){
         nanCounter = 1;
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
      }else{
         nanCounter = iBegin;
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
         *(p2+iBegin -2) = mynan;//must not be available for finding the cheapest path
      }
      
      iEnd   = j+ws+1;
      if(iEnd >= nx){
         iEnd = nx;   
      }else{
         *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
      }
      
      for (int i = iBegin; i < iEnd; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,j));
         if((*(p2+i) > threshold) | (*(p2+i) != *(p2+i))){
            *(p2+i) = mynan;
            nanCounter ++;
         }
      }
      
      if(nanCounter == nx) return(mynan);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}



// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_dtw2vec_cm_ws_inc(NumericVector gcm_lc, Rcpp::NumericMatrix cm, 
                        std::string step_pattern, int ws, int ny){
   
   // cm ... cost matrix for new observations of y
   // gcm_lc ... last column of old GCM
   // int ws ... window size
   // int ny ... length of time series y exclusive new observations,
   //             length(c(y,newObs)) = length(newObs) + ny
   
   int n = cm.nrow();
   int m = cm.ncol();// = nnewObs
   
   double * p1 = new double [n];
   double * p2 = new double [n];
   double * ptmp; 
   
   int iBegin = 0;
   int iEnd = 0;
   
   double mynan;
   NumericVector gcm_lr_new(m);
   NumericVector gcm_lc_new(n);
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   // set step sunction
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   
   if(n != gcm_lc.size()){
      return mynan;
   }
   
   // first column
   for(int i=0; i<n; i++){
      p1[i] = gcm_lc[i];
      p2[i] = mynan;//initialize b with NAN
   }
   
   for(int j=0; j < m; j++){
      iBegin = ny+j-ws;
      if(iBegin <= 0){
         *p2 = cm(0,j) + *(p1);
         iBegin = 1;
      }else if (iBegin == 1){
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
      }else{
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
         *(p2+iBegin -2) = mynan;//must not be available for finding the cheapest path
      }
      
      iEnd   = ny+j+ws+1;
      if(iEnd >= n){
         iEnd = n;
      }else{
         *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
      }
      
      for (int i = iBegin; i < iEnd; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,j));
      
      }
      gcm_lr_new[j] = *(p2+n-1);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   for(int i=0; i<n; i++){
      gcm_lc_new[i] = *(p1+i);
   }
   
   
   List ret;
   ret["gcm_lr_new"] = gcm_lr_new;
   ret["gcm_lc_new"] = gcm_lc_new;
   ret["distance"] = *(p1+n-1);
   delete[] p1;
   delete[] p2;
   return ret ;
   
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec (const arma::vec& x, const arma::vec& y, 
                    std::string step_pattern)
{
   
   int nx = x.size();
   int ny = y.size();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double ret;
   
   // first column
   *p1 = abs(x[0]-y[0]);
   for(int i=1; i<nx; i++){
      p1[i] = abs(x[i]-y[0]) + p1[i-1];
   }
   
   // set step function
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   
   for(int j=1; j < ny; j++){
      *p2 = abs(x[0]-y[j]) + *(p1);
      
      for(int i=1; i<nx; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), abs(x[i]-y[j]));
      }
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec_ws (const arma::vec& x, const arma::vec& y, 
                       std::string step_pattern, int ws)
{
   
   int iBegin = 0;
   int iEnd = 0;
   
   int nx = x.size();
   int ny = y.size();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double mynan;
   double ret;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   // set step sunction
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   
   //initialize a and b with NAN
   for(int i=0; i < nx; i++){
      p1[i] = mynan;
      p2[i] = mynan;
   }
   
   // first column
   p1[0] = abs(x[0]-y[0]);
   iEnd   = std::min(nx, ws+1);
   for(int i=1; i < iEnd; i++){
      p1[i] = abs(x[i]-y[0]) + p1[i-1];
   }
   
   for(int j=1; j < ny; j++){
      iBegin = j-ws;
      if(iBegin <= 0){
         *p2 = abs(x[0]-y[j]) + *(p1);
         iBegin = 1;   
      }else if (iBegin == 1){
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
      }else{
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
         *(p2+iBegin -2) = mynan;//must not be available for finding the cheapest path
      }
      
      
      iEnd   = j+ws+1;
      if(iEnd >= nx){
         iEnd = nx;   
      }else{
         *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
      }
      
      for (int i = iBegin; i < iEnd; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), abs(x[i]-y[j]));
      }
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec_ea (const arma::vec& x, const arma::vec& y, 
                       std::string step_pattern, double threshold)
{
   int nx = x.size();
   int ny = y.size();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double mynan;
   double ret;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   // set step sunction
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   
   // first column
   p1[0] = abs(x[0]-y[0]);
   if(p1[0] > threshold)  return(mynan); 
   
   for(int i=1; i<nx; i++){
      p1[i] = abs(x[i]-y[0]) + p1[i-1];
      if(p1[i] > threshold) p1[i] = mynan;
   }

   for(int j=1; j < ny; j++){
      
      nanCounter = 0;
      *p2 = abs(x[0]-y[j]) + *(p1);
      if(*(p2) > threshold){
         *(p2) = mynan;
         nanCounter ++;
      }
      
      for(int i=1; i<nx; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), abs(x[i]-y[j]));
         
         if((*(p2+i) > threshold) | (*(p2+i) != *(p2+i))){
            *(p2+i) = mynan;
            nanCounter ++;
         } 
      }
      if(nanCounter == nx) return(mynan);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec_ws_ea (const arma::vec& x, const arma::vec& y, 
                          std::string step_pattern, int ws, double threshold)
{
   
   int iBegin = 0;
   int iEnd = 0;
   
   int nx = x.size();
   int ny = y.size();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double mynan;
   double ret;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   // set step sunction
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   
   //initialize a and b with NAN
   for(int i=0; i < nx; i++){
      p1[i] = mynan;
      p2[i] = mynan;
   }
   
   // first column
   p1[0] = abs(x[0]-y[0]);
   if(p1[0] > threshold)  return(mynan); 
   
   iEnd   = std::min(nx, ws+1);
   for(int i=1; i < iEnd; i++){
      p1[i] = abs(x[i]-y[0]) + p1[i-1];
      if(p1[i] > threshold) p1[i] = mynan;
   }
   
   
   for(int j=1; j < ny; j++){
      nanCounter = 0;
      iBegin = j-ws;
      if(iBegin <= 0){
         *p2 = abs(x[0]-y[j]) + *(p1);
         if(*(p2) > threshold){
            *(p2) = mynan;
            nanCounter ++;
         }
         iBegin = 1;   
      }else if (iBegin == 1){
         nanCounter = 1;
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
      }else{
         nanCounter = iBegin;
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
         *(p2+iBegin -2) = mynan;//must not be available for finding the cheapest path
      }
      
      iEnd   = j+ws+1;
      if(iEnd >= nx){
         iEnd = nx;   
      }else{
         *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
      }
      
      for (int i = iBegin; i < iEnd; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), abs(x[i]-y[j]));
         if((*(p2+i) > threshold) | (*(p2+i) != *(p2+i))){
            *(p2+i) = mynan;
            nanCounter ++;
         }
      }
      
      if(nanCounter == nx) return(mynan);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec_mv (const arma::mat& x, const arma::mat& y, 
                       std::string step_pattern, std::string dist_method)
{
   int ncol  = x.n_cols;
   int nx = x.n_rows;
   int ny = y.n_rows;
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp;
   double ret;

   SEXP dist_SEXP = select_dist(dist_method);
   XPtr<funcPtr_dist> xpfun(dist_SEXP);
   funcPtr_dist dist_fun = *xpfun;

   // set step sunction
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   
   // set step function with std::function
   // vecStepFunction step_fun = selectVecStepFunction(step_pattern);
   // distFunction dist_fun = selectDistFunction(dist_method);
   
   // first column
   *p1 = dist_fun(x, y, 0, 0, ncol);
   for(int i=1; i<nx; i++){
      p1[i] = dist_fun(x, y, i, 0, ncol) + p1[i-1];
   }


   for(int j=1; j < ny; j++){
      *p2 = dist_fun(x, y, 0, j, ncol) + *(p1);

      for(int i=1; i<nx; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), dist_fun(x, y, i, j, ncol));
      }
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }

   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;

   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec_mv_ws_ea (const arma::mat& x, const arma::mat& y, 
                             std::string step_pattern, std::string dist_method, int ws, double threshold)
{

   int iBegin = 0;
   int iEnd = 0;

   int ncol  = x.n_cols;
   int nx = x.n_rows;
   int ny = y.n_rows;
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp;
   double mynan;
   double ret;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();

   SEXP dist_SEXP = select_dist(dist_method);
   XPtr<funcPtr_dist> xpfun(dist_SEXP);
   funcPtr_dist dist_fun = *xpfun;


   // set step sunction
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   


   //initialize a and b with NAN
   for(int i=0; i < nx; i++){
      p1[i] = mynan;
      p2[i] = mynan;
   }

   // first column
   p1[0] = dist_fun(x, y, 0, 0, ncol);
   if(p1[0] > threshold)  return(mynan);

   iEnd   = std::min(nx, ws+1);
   for(int i=1; i < iEnd; i++){
      p1[i] =  dist_fun(x, y, i, 0, ncol) + p1[i-1];
      if(p1[i] > threshold) p1[i] = mynan;
   }


   for(int j=1; j < ny; j++){
      nanCounter = 0;
      iBegin = j-ws;
      if(iBegin <= 0){
         *p2 = dist_fun(x, y, 0, j, ncol) + *(p1);
         if(*(p2) > threshold){
            *(p2) = mynan;
            nanCounter ++;
         }
         iBegin = 1;
      }else if (iBegin == 1){
         nanCounter = 1;
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
      }else{
         nanCounter = iBegin;
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
         *(p2+iBegin -2) = mynan;//must not be available for finding the cheapest path
      }

      iEnd   = j+ws+1;
      if(iEnd >= nx){
         iEnd = nx;
      }else{
         *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
      }

      for (int i = iBegin; i < iEnd; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), dist_fun(x, y, i, j, ncol));
         if((*(p2+i) > threshold) | (*(p2+i) != *(p2+i))){
            *(p2+i) = mynan;
            nanCounter ++;
         }
      }

      if(nanCounter == nx) return(mynan);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;

   return (ret);
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec_v32 (const arma::vec& x, const arma::vec& y)
{
   
   int nx = x.size();
   int ny = y.size();
   int i;
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double ret;
   
   // first column
   p1[0] = abs(x[0]-y[0]);
   for(i=1; i<nx; i++){
      p1[i] = abs(x[i]-y[0]) + p1[i-1];
   }
   
   for(int j=1; j < ny; j++){
      *p2 = abs(x[0]-y[j]) + *(p1);
      
      for( i=1; i < nx-32; i+= 32){
         *(p2+i)   = abs(x[i]-y[j])   + mymin(*(p2+i-1), mymin(*(p1+i)  , *(p1+i-1)));
         *(p2+i+1) = abs(x[i+1]-y[j]) + mymin(*(p2+i)  , mymin(*(p1+i+1), *(p1+i)));
         *(p2+i+2) = abs(x[i+2]-y[j]) + mymin(*(p2+i+1), mymin(*(p1+i+2), *(p1+i+1)));
         *(p2+i+3) = abs(x[i+3]-y[j]) + mymin(*(p2+i+2), mymin(*(p1+i+3), *(p1+i+2)));
         *(p2+i+4) = abs(x[i+4]-y[j]) + mymin(*(p2+i+3), mymin(*(p1+i+4), *(p1+i+3)));
         *(p2+i+5) = abs(x[i+5]-y[j]) + mymin(*(p2+i+4), mymin(*(p1+i+5), *(p1+i+4)));
         *(p2+i+6) = abs(x[i+6]-y[j]) + mymin(*(p2+i+5), mymin(*(p1+i+6), *(p1+i+5)));
         *(p2+i+7) = abs(x[i+7]-y[j]) + mymin(*(p2+i+6), mymin(*(p1+i+7), *(p1+i+6)));
         *(p2+i+8) = abs(x[i+8]-y[j]) + mymin(*(p2+i+7), mymin(*(p1+i+8), *(p1+i+7)));
         *(p2+i+9) = abs(x[i+9]-y[j]) + mymin(*(p2+i+8), mymin(*(p1+i+9), *(p1+i+8)));
         *(p2+i+10) = abs(x[i+10]-y[j]) + mymin(*(p2+i+9), mymin(*(p1+i+10), *(p1+i+9)));
         *(p2+i+11) = abs(x[i+11]-y[j]) + mymin(*(p2+i+10), mymin(*(p1+i+11), *(p1+i+10)));
         *(p2+i+12) = abs(x[i+12]-y[j]) + mymin(*(p2+i+11), mymin(*(p1+i+12), *(p1+i+11)));
         *(p2+i+13) = abs(x[i+13]-y[j]) + mymin(*(p2+i+12), mymin(*(p1+i+13), *(p1+i+12)));
         *(p2+i+14) = abs(x[i+14]-y[j]) + mymin(*(p2+i+13), mymin(*(p1+i+14), *(p1+i+13)));
         *(p2+i+15) = abs(x[i+15]-y[j]) + mymin(*(p2+i+14), mymin(*(p1+i+15), *(p1+i+14)));
         
         *(p2+i+16) = abs(x[i+16]-y[j]) + mymin(*(p2+i+15), mymin(*(p1+i+16), *(p1+i+15)));
         *(p2+i+17) = abs(x[i+17]-y[j]) + mymin(*(p2+i+16), mymin(*(p1+i+17), *(p1+i+16)));
         *(p2+i+18) = abs(x[i+18]-y[j]) + mymin(*(p2+i+17), mymin(*(p1+i+18), *(p1+i+17)));
         *(p2+i+19) = abs(x[i+19]-y[j]) + mymin(*(p2+i+18), mymin(*(p1+i+19), *(p1+i+18)));
         *(p2+i+20) = abs(x[i+20]-y[j]) + mymin(*(p2+i+19), mymin(*(p1+i+20), *(p1+i+19)));
         *(p2+i+21) = abs(x[i+21]-y[j]) + mymin(*(p2+i+20), mymin(*(p1+i+21), *(p1+i+20)));
         *(p2+i+22) = abs(x[i+22]-y[j]) + mymin(*(p2+i+21), mymin(*(p1+i+22), *(p1+i+21)));
         *(p2+i+23) = abs(x[i+23]-y[j]) + mymin(*(p2+i+22), mymin(*(p1+i+23), *(p1+i+22)));
         
         *(p2+i+24) = abs(x[i+24]-y[j]) + mymin(*(p2+i+23), mymin(*(p1+i+24), *(p1+i+23)));
         *(p2+i+25) = abs(x[i+25]-y[j]) + mymin(*(p2+i+24), mymin(*(p1+i+25), *(p1+i+24)));
         *(p2+i+26) = abs(x[i+26]-y[j]) + mymin(*(p2+i+25), mymin(*(p1+i+26), *(p1+i+25)));
         *(p2+i+27) = abs(x[i+27]-y[j]) + mymin(*(p2+i+26), mymin(*(p1+i+27), *(p1+i+26)));
         *(p2+i+28) = abs(x[i+28]-y[j]) + mymin(*(p2+i+27), mymin(*(p1+i+28), *(p1+i+27)));
         *(p2+i+29) = abs(x[i+29]-y[j]) + mymin(*(p2+i+28), mymin(*(p1+i+29), *(p1+i+28)));
         *(p2+i+30) = abs(x[i+30]-y[j]) + mymin(*(p2+i+29), mymin(*(p1+i+30), *(p1+i+29)));
         *(p2+i+31) = abs(x[i+31]-y[j]) + mymin(*(p2+i+30), mymin(*(p1+i+31), *(p1+i+30)));
         
      }
      for (; i < nx; i++) {
         *(p2+i)   = abs(x[i]-y[j])   + mymin(*(p2+i-1), mymin(*(p1+i)  , *(p1+i-1)));
      }
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_dtw2vec_inc (NumericVector x, NumericVector newObs, NumericVector gcm_lc,
                      std::string step_pattern)
{
   // x ... time series that is fixed
   // y ... time series of new observations, ONLY new observations to be appended
   // gcm_lc ... last column of old GCM
   
   int nx = x.size();
   int nnewObs = newObs.size();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp;
   double mynan;
   NumericVector gcm_lr_new(nnewObs);
   NumericVector gcm_lc_new(nx);
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   // set step sunction
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   
   if(nx != gcm_lc.size()){
      return mynan;
   }
   
   // first column
   for(int i=0; i<nx; i++){
      p1[i] = gcm_lc[i];
   }
   
   for(int j=0; j < nnewObs; j++){
      *p2 = abs(x[0]-newObs[j]) + *(p1);
      
      for(int i=1; i<nx; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), abs(x[i]-newObs[j]));
      }
      gcm_lr_new[j] = *(p2+nx-1);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   for(int i=0; i<nx; i++){
      gcm_lc_new[i] = *(p1+i);
   }
   
   
   List ret;
   ret["gcm_lr_new"] = gcm_lr_new;
   ret["gcm_lc_new"] = gcm_lc_new;
   ret["distance"] = *(p1+nx-1);
   delete[] p1;
   delete[] p2;
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_dtw2vec_inc_ws (NumericVector x, NumericVector newObs, NumericVector gcm_lc,
                         int ws, int ny, std::string step_pattern)
{
   // x ... time series that is fixed
   // newObs ... time series of new observations, ONLY new observations to be appended
   // gcm_lc ... last column of old GCM
   // int ws ... window size
   // int ny ... length of time series y exclusive new observations,
   //             length(c(y,newObs)) = length(newObs) + ny
   
   int nx = x.size();
   int nnewObs = newObs.size();
   int iBegin = 0;
   int iEnd = 0;
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp;
   double mynan;
   NumericVector gcm_lr_new(nnewObs);
   NumericVector gcm_lc_new(nx);
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   // set step sunction
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   
   if(nx != gcm_lc.size()){
      return mynan;
   }
   
   // first column
   for(int i=0; i<nx; i++){
      p1[i] = gcm_lc[i];
      p2[i] = mynan;//initialize b with NAN
   }
   
   for(int j=0; j < nnewObs; j++){
      iBegin = ny+j-ws;
      if(iBegin <= 0){
         *p2 = abs(x[0]-newObs[j]) + *(p1);
         iBegin = 1;
      }else if (iBegin == 1){
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
      }else{
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
         *(p2+iBegin -2) = mynan;//must not be available for finding the cheapest path
      }
      
      iEnd   = ny+j+ws+1;
      if(iEnd >= nx){
         iEnd = nx;
      }else{
         *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
      }
      
      for (int i = iBegin; i < iEnd; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), abs(x[i]-newObs[j]));
      }
      gcm_lr_new[j] = *(p2+nx-1);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   for(int i=0; i<nx; i++){
      gcm_lc_new[i] = *(p1+i);
   }
   
   List ret;
   ret["gcm_lr_new"] = gcm_lr_new;
   ret["gcm_lc_new"] = gcm_lc_new;
   ret["distance"] = *(p1+nx-1);
   delete[] p1;
   delete[] p2;
   
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_dtw2vec_inc_mv (const arma::mat& x, const arma::mat& newObs,
                         const arma::vec& gcm_lc, std::string dist_method,
                         std::string step_pattern)
{
   
   int ncol  = x.n_cols;
   int nx = x.n_rows;
   int nnewObs = newObs.n_rows;
   int lcs = gcm_lc.size();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp;
   double mynan;
   
   NumericVector gcm_lr_new(nnewObs);
   NumericVector gcm_lc_new(nx);
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   SEXP dist_SEXP = select_dist(dist_method);
   XPtr<funcPtr_dist> xpfun(dist_SEXP);
   funcPtr_dist dist_fun = *xpfun;
   
   // set step sunction
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   
   if(nx != lcs ){
      return mynan;
   }
   
   // first column
   for(int i=0; i<nx; i++){
      p1[i] = gcm_lc[i];
   }
   
   for(int j=0; j < nnewObs; j++){
      *p2 = dist_fun(x, newObs, 0, j, ncol) + *(p1);
      
      for(int i=1; i<nx; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), dist_fun(x, newObs, i, j, ncol));
      }
      gcm_lr_new[j] = *(p2+nx-1);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   for(int i=0; i<nx; i++){
      gcm_lc_new[i] = *(p1+i);
   }
   
   
   List ret;
   ret["gcm_lr_new"] = gcm_lr_new;
   ret["gcm_lc_new"] = gcm_lc_new;
   ret["distance"] = *(p1+nx-1);
   delete[] p1;
   delete[] p2;
   return ret ;
}



// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_dtw2vec_inc_mv_ws (const arma::mat& x, const arma::mat& newObs,
                            const arma::vec& gcm_lc, std::string dist_method,
                            int ws, int ny, std::string step_pattern)
{
   
   int ncol  = x.n_cols;
   int nx = x.n_rows;
   int nnewObs = newObs.n_rows;
   int lcs = gcm_lc.size();
   int iBegin = 0;
   int iEnd = 0;
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp;
   double mynan;
   
   NumericVector gcm_lr_new(nnewObs);
   NumericVector gcm_lc_new(nx);
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   SEXP dist_SEXP = select_dist(dist_method);
   XPtr<funcPtr_dist> xpfun(dist_SEXP);
   funcPtr_dist dist_fun = *xpfun;
   
   // set step sunction
   SEXP step_SEXP = selectVecStep(step_pattern);
   XPtr<funcPtr_step_vec> xpfun_step(step_SEXP);
   funcPtr_step_vec step_fun = *xpfun_step;
   
   if(nx != lcs){
      return mynan;
   }
   
   // first column
   for(int i=0; i<nx; i++){
      p1[i] = gcm_lc[i];
      p2[i] = mynan;//initialize b with NAN
   }
   
   for(int j=0; j < nnewObs; j++){
      iBegin = ny+j-ws;
      if(iBegin <= 0){
         *p2 = dist_fun(x, newObs, 0, j, ncol) + *(p1);
         iBegin = 1;
      }else if (iBegin == 1){
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
      }else{
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
         *(p2+iBegin -2) = mynan;//must not be available for finding the cheapest path
      }
      
      iEnd   = ny+j+ws+1;
      if(iEnd >= nx){
         iEnd = nx;
      }else{
         *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
      }
      
      for (int i = iBegin; i < iEnd; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), dist_fun(x, newObs, i, j, ncol));
      }
      gcm_lr_new[j] = *(p2+nx-1);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   for(int i=0; i<nx; i++){
      gcm_lc_new[i] = *(p1+i);
   }
   
   
   List ret;
   ret["gcm_lr_new"] = gcm_lr_new;
   ret["gcm_lc_new"] = gcm_lc_new;
   ret["distance"] = *(p1+nx-1);
   delete[] p1;
   delete[] p2;
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


double multp_dtw2vec_ws_ea (const arma::vec& x, const arma::vec& y, 
                          std::string step_pattern, int ws, double threshold)
{
  
   int iBegin = 0;
   int iEnd = 0;
   
   int nx = x.size();
   int ny = y.size();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double mynan;
   double ret;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   // set step function with std::function
   vecStepFunction step_fun = selectVecStepFunction(step_pattern);
   
   //initialize a and b with NAN
   for(int i=0; i < nx; i++){
      p1[i] = mynan;
      p2[i] = mynan;
   }
   
   // first column
   p1[0] = abs(x[0]-y[0]);
   if(p1[0] > threshold)  return(mynan); 
   
   iEnd   = std::min(nx, ws+1);
   for(int i=1; i < iEnd; i++){
      p1[i] = abs(x[i]-y[0]) + p1[i-1];
      if(p1[i] > threshold) p1[i] = mynan;
   }
   
   
   
   for(int j=1; j < ny; j++){
      nanCounter = 0;
      iBegin = j-ws;
      if(iBegin <= 0){
         *p2 = abs(x[0]-y[j]) + *(p1);
         if(*(p2) > threshold){
            *(p2) = mynan;
            nanCounter ++;
         }
         iBegin = 1;   
      }else if (iBegin == 1){
         nanCounter = 1;
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
      }else{
         nanCounter = iBegin;
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
         *(p2+iBegin -2) = mynan;//must not be available for finding the cheapest path
      }
      
      iEnd   = j+ws+1;
      if(iEnd >= nx){
         iEnd = nx;   
      }else{
         *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
      }
      
      for (int i = iBegin; i < iEnd; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), abs(x[i]-y[j]));
         if((*(p2+i) > threshold) | (*(p2+i) != *(p2+i))){
            *(p2+i) = mynan;
            nanCounter ++;
         }
      }
      
      if(nanCounter == nx) return(mynan);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


double multp_dtw2vec_mv_ws_ea (const arma::mat& x, const arma::mat& y, 
                             std::string step_pattern, std::string dist_method, int ws, double threshold)
{
   
   
   int iBegin = 0;
   int iEnd = 0;
   
   int ncol  = x.n_cols;
   int nx = x.n_rows;
   int ny = y.n_rows;
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp;
   double mynan;
   double ret;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   distFunction dist_fun = selectDistFunction(dist_method);
   vecStepFunction step_fun = selectVecStepFunction(step_pattern);
   
   
   //initialize a and b with NAN
   for(int i=0; i < nx; i++){
      p1[i] = mynan;
      p2[i] = mynan;
   }
   
   // first column
   p1[0] = dist_fun(x, y, 0, 0, ncol);
   if(p1[0] > threshold)  return(mynan);
   
   iEnd   = std::min(nx, ws+1);
   for(int i=1; i < iEnd; i++){
      p1[i] =  dist_fun(x, y, i, 0, ncol) + p1[i-1];
      if(p1[i] > threshold) p1[i] = mynan;
   }
   
   
   for(int j=1; j < ny; j++){
      nanCounter = 0;
      iBegin = j-ws;
      if(iBegin <= 0){
         *p2 = dist_fun(x, y, 0, j, ncol) + *(p1);
         if(*(p2) > threshold){
            *(p2) = mynan;
            nanCounter ++;
         }
         iBegin = 1;
      }else if (iBegin == 1){
         nanCounter = 1;
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
      }else{
         nanCounter = iBegin;
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
         *(p2+iBegin -2) = mynan;//must not be available for finding the cheapest path
      }
      
      iEnd   = j+ws+1;
      if(iEnd >= nx){
         iEnd = nx;
      }else{
         *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
      }
      
      for (int i = iBegin; i < iEnd; i++){
         *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), dist_fun(x, y, i, j, ncol));
         if((*(p2+i) > threshold) | (*(p2+i) != *(p2+i))){
            *(p2+i) = mynan;
            nanCounter ++;
         }
      }
      
      if(nanCounter == nx) return(mynan);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}




