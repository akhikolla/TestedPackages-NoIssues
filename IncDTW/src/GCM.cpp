#include <RcppArmadillo.h>
#include "IncDTW.h"
using namespace Rcpp;
using namespace std;

XPtr<funcPtr_dist> select_dist2(std::string dist_method) {

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


struct gcmOneStep {
   double g;
   int d;
} ;


gcmOneStep gcm_step_symm1(double gcm10,  // vertical
                           double gcm11, // diagonal
                           double gcm01, // horizontal
                           double cm00)  // cost
{
   // z > nan for z != nan is required by C the standard
   int nan10 = isnan(gcm10), nan01 = isnan(gcm01);
   gcmOneStep ret;
   
   if(!nan10 && !nan01){
      
      if(gcm11 <= gcm10 && gcm11 <= gcm01){
         ret.g = cm00 + gcm11;
         ret.d = 1;
      } else if(gcm10 <= gcm11 && gcm10 <= gcm01){
         ret.g = cm00 + gcm10;
         ret.d = 3;
      }else{
         ret.g = cm00 + gcm01;
         ret.d = 2;
      }
      
   } else if(nan10 && nan01){
      
      ret.g = cm00 + gcm11;
      ret.d = 1;
      
   } else if (nan10){
      
      if(gcm11 <= gcm01){
         ret.g = cm00 + gcm11;
         ret.d = 1;
      } else{
         ret.g = cm00 + gcm01;
         ret.d = 2;
      }
      
   } else{// if (std::isnan(gcm(i, j-1))){
      
      if(gcm11 <= gcm10){
         ret.g = cm00 + gcm11;
         ret.d = 1;
      } else{
         ret.g = cm00 + gcm10;
         ret.d = 3;
      }
      
   }
   
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


gcmOneStep gcm_step_symm2(double gcm10,  // vertical
                          double gcm11, // diagonal
                          double gcm01, // horizontal
                          double cm00)  // cost
{
   // z > nan for z != nan is required by C the standard
   int nan10 = isnan(gcm10), nan01 = isnan(gcm01);
   gcmOneStep ret;
   
   if(!nan10 && !nan01){
      
      if(gcm11 + cm00 <= gcm10 && gcm11 + cm00 <= gcm01){
         ret.g = 2 * cm00 + gcm11;
         ret.d = 1;
      } else if(gcm10 <= gcm11 + cm00 && gcm10 <= gcm01){
         ret.g = cm00 + gcm10;
         ret.d = 3;
      }else{
         ret.g = cm00 + gcm01;
         ret.d = 2;
      }
      
   } else if(nan10 && nan01){
      
      ret.g = 2 * cm00 + gcm11;
      ret.d = 1;
      
   } else if (nan10){
      
      if(gcm11 + cm00 <= gcm01){
         ret.g = 2 * cm00 + gcm11;
         ret.d = 1;
      } else{
         ret.g = cm00 + gcm01;
         ret.d = 2;
      }
      
   } else{// if (std::isnan(gcm(i, j-1))){
      
      if(gcm11 + cm00 <= gcm10){
         ret.g = 2 * cm00 + gcm11;
         ret.d = 1;
      } else{
         ret.g = cm00 + gcm10;
         ret.d = 3;
      }
      
   }
   
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


typedef gcmOneStep (*funcPtr_step)(const double gcm10, 
                                  const double gcm11, 
                                  const double gcm01, 
                                  const double cm00);


XPtr<funcPtr_step> selectGcmStep(std::string step_pattern) {
   if (step_pattern == "symmetric1")
      return(XPtr<funcPtr_step>(new funcPtr_step(&gcm_step_symm1)));
   else if (step_pattern == "symmetric2")
      return(XPtr<funcPtr_step>(new funcPtr_step(&gcm_step_symm2)));
   else
      return XPtr<funcPtr_step>(R_NilValue); // runtime error as NULL no XPtr
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
NumericVector cpp_znorm(NumericVector x, double sd_threshold,
                        Rcpp::Nullable< Rcpp::NumericVector > mu_in = R_NilValue,
                        Rcpp::Nullable< Rcpp::NumericVector > sd_in = R_NilValue)
 {
   
   int N = x.size();
   NumericVector y(N);
   double mu;
   double sd;
   
   if (mu_in.isNull() && sd_in.isNull()) {
      mu = 0;
      sd = 0;
      
      for(int i = 0; i < N; i++ ){
         mu += x[i];
      }
      mu /= N;
   
      for(int i = 0; i < N; i++ ){
         y[i] = (x[i] - mu);
         sd += y[i] * y[i];
      }
      sd /= (N-1);
      sd = sqrt(sd);
      // Rcout << "\nsd: " << sd;
      if(sd > sd_threshold){
         for(int i = 0; i < N; i++ ){
            y[i] /= sd;
         }
      }
      
   }else if (mu_in.isNull() && !sd_in.isNull()) {
      mu = 0;
      sd = Rcpp::as< Rcpp::NumericVector >(sd_in)[0];
      
      for(int i = 0; i < N; i++ ){
         mu += x[i];
      }
      mu /= N;
      
      if(sd > sd_threshold){
         for(int i = 0; i < N; i++ ){
            y[i] = (x[i] - mu)/sd;
         }
      }else{
         for(int i = 0; i < N; i++ ){
            y[i] = (x[i] - mu);
         }
      }
      
   }else if (!mu_in.isNull() && sd_in.isNull()) {
      sd = 0;
      mu = Rcpp::as< Rcpp::NumericVector >(mu_in)[0];
      
      for(int i = 0; i < N; i++ ){
         y[i] = (x[i] - mu);
         sd += y[i] * y[i];
      }
      sd /= (N-1);
      sd = sqrt(sd);
      
      if(sd > sd_threshold){
         for(int i = 0; i < N; i++ ){
            y[i] /= sd;
         }
      }
      
   }else{
      
      mu = Rcpp::as< Rcpp::NumericVector >(mu_in)[0];
      sd = Rcpp::as< Rcpp::NumericVector >(sd_in)[0];
      
      if(sd > sd_threshold){
         for(int i = 0; i < N; i++ ){
            y[i] = (x[i] - mu)/sd;
         }
      }else{
         for(int i = 0; i < N; i++ ){
            y[i] = (x[i] - mu);
         }
      }
   }
   
   return y;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
NumericVector cpp_norm01(NumericVector x, double sd_threshold,
                         Rcpp::Nullable< Rcpp::NumericVector > min_in = R_NilValue,
                         Rcpp::Nullable< Rcpp::NumericVector > max_in = R_NilValue)
{
   int N = x.size();
   double denom;
   NumericVector y(N);
   
   double xmin;
   double xmax;
   
   if (min_in.isNull() && max_in.isNull()) {
      xmin = x[0];
      xmax = x[0];
      for(int i = 1; i < N; i++ ){
         if( x[i] < xmin){
            xmin = x[i];
         }
         if( x[i] > xmax){
            xmax = x[i];
         };
      }
      
   } else if (min_in.isNull() && !max_in.isNull()) {
      xmin = x[0];
      xmax = Rcpp::as< Rcpp::NumericVector >(max_in)[0];
      
      for(int i = 1; i < N; i++ ){
         if( x[i] < xmin){
            xmin = x[i];
         }
      }
      
   }else if (!min_in.isNull() && max_in.isNull()) {
      xmax = x[0];
      xmin = Rcpp::as< Rcpp::NumericVector >(min_in)[0];
      
      for(int i = 1; i < N; i++ ){
         if( x[i] > xmax){
            xmax = x[i];
         };
      }
      
   }else{
      xmin = Rcpp::as< Rcpp::NumericVector >(min_in)[0];
      xmax = Rcpp::as< Rcpp::NumericVector >(max_in)[0];
   }
   
   denom = xmax - xmin;
   if(denom > sd_threshold){
      for(int i = 0; i < N; i++ ){
         y[i] = (x[i] - xmin)/ denom;
      }
   }else{
      for(int i = 0; i < N; i++ ){
         y[i] = (x[i] - xmin);
      }
   }
   
   
   return y;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
NumericMatrix cpp_cm(const arma::mat& x,
                     const arma::mat& y,
                     std::string dist_method, int ws, int nPrevObs)
{
   int nx = x.n_rows;
   int ncol = x.n_cols;
   int ny = y.n_rows;
   int iBegin;
   int iEnd;
   
   
   SEXP dist_SEXP = select_dist2(dist_method);
   XPtr<funcPtr_dist> xpfun(dist_SEXP);
   funcPtr_dist dist_fun = *xpfun;
   
   NumericMatrix cm(nx, ny);
   if(ws == -1){
      for(int j = 0; j < ny; j++){
         for(int i = 0; i < nx; i++){
            cm(i, j) = dist_fun(x, y, i, j, ncol);
         }
      }
   } else{
      std::fill( cm.begin(), cm.end(), NumericVector::get_na() );
      for (int j = 0; j < ny; j++){
         iBegin = std::max(0 , nPrevObs + j - ws);
         iEnd   = std::min(nx, nPrevObs + j + ws + 1);
         
         for (int i = iBegin; i < iEnd; i++){
            cm(i, j) = dist_fun(x, y, i, j, ncol);
         }
      }
   }

   return cm ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
NumericMatrix cpp_diffm(const NumericVector& x,
                        const NumericVector& y,
                        int ws, int nPrevObs){
   
   int nx = x.size();
   int ny = y.size();
   int iBegin;
   int iEnd;
   
   NumericMatrix diffm(nx, ny);
   if(ws == -1){
      for(int j = 0; j < ny; j++){
         for(int i = 0; i < nx; i++){
            diffm(i, j) = x[i] - y[j];
         }
      }
   } else{
      std::fill( diffm.begin(), diffm.end(), NumericVector::get_na() );
      for (int j = 0; j < ny; j++){
         iBegin = std::max(0 , nPrevObs + j - ws);
         iEnd   = std::min(nx, nPrevObs + j + ws + 1);
         
         for (int i = iBegin; i < iEnd; i++){
            diffm(i, j) = x[i] - y[j];
         }
      }
   }

   return diffm ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List GCM_Sakoe_cpp(Rcpp::NumericMatrix cM, int ws, std::string step_pattern){
   // Global Cost Matrix
   int n = cM.nrow();
   int m = cM.ncol();
   int iBegin = 0;
   int iEnd = 0;
   int wsI = 0;
   int wsJ = 0;
   double cost = 0;
   gcmOneStep tmp;
   
   SEXP step_SEXP = selectGcmStep(step_pattern);
   XPtr<funcPtr_step> xpfun(step_SEXP);
   funcPtr_step step_fun = *xpfun;

   NumericMatrix gcm(n, m);
   IntegerMatrix dm(n, m);
   std::fill( gcm.begin(), gcm.end(), NumericVector::get_na() );
   // std::fill( dm.begin(), dm.end(), NumericVector::get_na() );
   std::fill( dm.begin(), dm.end(), NA_INTEGER );
   
   wsI = std::min(n, ws+1);
   wsJ = std::min(m, ws+1);
   
   gcm(0,0) = cM(0,0);
   for(int i =1; i < wsI; i++){
      dm(i,0)=3;
      gcm(i,0) = cM(i, 0) + gcm(i - 1, 0);
   }
   for(int j =1; j < wsJ; j++){
      dm(0,j)=2;
      gcm(0, j) = cM(0, j) + gcm(0, j - 1);
   }
   
   for (int j = 1; j < m; j++){
      iBegin = std::max(1, j-ws);
      iEnd   = std::min(n, j+ws+1);

      for (int i = iBegin; i < iEnd; i++){
         cost = cM(i,j);
         tmp = step_fun(gcm(i-1, j), 
                        gcm(i-1, j-1),
                        gcm(i  , j-1),
                        cost);
         gcm(i, j) = tmp.g;
         dm(i,j) = tmp.d;
      }
   }

   List ret;
   ret["gcm"] = gcm;
   ret["dm"] = dm;
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List IGCM_Sakoe_cpp(Rcpp::NumericMatrix gcmN, //global costmatrix with new empty columns
                    Rcpp::IntegerMatrix dmN, //direction matrix with new empty columns
                    Rcpp::NumericMatrix cmN,//local cost matrix of new observations and old constant vector
                    int ws,
                    std::string step_pattern){ 
   
   
   // C has one to many new columns that are not covered by the global cost matrix gcm so far
   // indexC gives the first index of the new observations
   // 
   int n = gcmN.nrow();
   int m = gcmN.ncol();
   int j;
   int Nnew = cmN.ncol();//number of new observations
   double cost;
   int iBegin;
   int iEnd;
   int wsJ = 0;
   gcmOneStep tmp;
   
   SEXP step_SEXP = selectGcmStep(step_pattern);
   XPtr<funcPtr_step> xpfun(step_SEXP);
   funcPtr_step step_fun = *xpfun;
   
   wsJ = std::min(m, ws+1);
   
   // first row
   if((m-Nnew) < ws){
      j = m-Nnew;
      do{
         gcmN(0,j) = gcmN(0,j-1) + cmN(0,(j-m+Nnew)); 
         dmN(0,j)  = 2;   
         j = j + 1;
      } while (j < wsJ);
   }
   
   // remaining
   
   for (int j = (m-Nnew); j < m; j++){
      iBegin = std::max(1, j-ws);
      iEnd   = std::min(n, j+ws+1);
      
      for (int i = iBegin; i < iEnd; i++){
         
         cost = cmN(i,(j-m+Nnew));
         tmp = step_fun(gcmN(i-1, j), 
                        gcmN(i-1, j-1),
                        gcmN(i  , j-1),
                        cost);
         gcmN(i, j) = tmp.g;
         dmN(i,j) = tmp.d;
      }
   }
   
   List ret;
   ret["gcm"] = gcmN;
   ret["dm"] = dmN;
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List GCM_cpp(Rcpp::NumericMatrix cM, std::string step_pattern){
   // Global Cost Matrix
   int n = cM.nrow();
   int m = cM.ncol();
   gcmOneStep tmp;
   
   SEXP step_SEXP = selectGcmStep(step_pattern);
   XPtr<funcPtr_step> xpfun(step_SEXP);
   funcPtr_step step_fun = *xpfun;
   
   double cost;
   //std::vector<double> lc;//local costs
   NumericMatrix gcm(n, m);
   IntegerMatrix dm(n, m);
   
   gcm(0,0) = cM(0,0);
   for(int i =1; i < n; i++){
      dm(i,0)=3;
      gcm(i,0) = cM(i, 0) + gcm(i - 1, 0);
   }
   for(int j =1; j < m; j++){
      dm(0,j)=2;
      gcm(0, j) = cM(0, j) + gcm(0, j - 1);
   }
   // dm(0,0) = NAN;
   dm(0,0) = NA_INTEGER;
   
   for (int i = 1; i < n; i++){
      for (int j = 1; j < m; j++){
         cost = cM(i,j);
         tmp = step_fun(gcm(i-1, j), 
                        gcm(i-1, j-1),
                        gcm(i  , j-1),
                        cost);
         gcm(i, j) = tmp.g;
         dm(i,j) = tmp.d;
         
      }
   }
   
   //List z = List::create( DTW ) ;
   List ret;
   ret["gcm"] = gcm;
   ret["dm"] = dm;
   return ret ;
}




// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// [[Rcpp::export]]
List IGCM_cpp(Rcpp::NumericMatrix gcmN, //global costmatrix with new empty columns
             Rcpp::IntegerMatrix dmN, //direction matrix with new empty columns
             Rcpp::NumericMatrix cmN,
             std::string step_pattern){ //local cost matrix of new observations and old constant vector
          
          
   // C has one to many new columns that are not covered by the global cost matrix gcm so far
   // indexC gives the first index of the new observations
   // 
   int n = gcmN.nrow();
   int m = gcmN.ncol();
   int Nnew = cmN.ncol();//number of new observations
   double cost;
   
   gcmOneStep tmp;
   
   SEXP step_SEXP = selectGcmStep(step_pattern);
   XPtr<funcPtr_step> xpfun(step_SEXP);
   funcPtr_step step_fun = *xpfun;
   
   for (int j = (m-Nnew); j < m; j++){
      gcmN(0,j) = gcmN(0,j-1) + cmN(0,(j-m+Nnew)); 
      dmN(0,j)  = 2;
   }
   //cost = n*m;
   
   for (int j = (m-Nnew); j < m; j++){
      for (int i = 1; i < n; i++){
         cost = cmN(i,(j-m+Nnew));
         tmp = step_fun(gcmN(i-1, j), 
                        gcmN(i-1, j-1),
                        gcmN(i  , j-1),
                        cost);
         gcmN(i, j) = tmp.g;
         dmN(i,j) = tmp.d;
      }
   }
   
   List ret;
   ret["gcm"] = gcmN;
   ret["dm"] = dmN;
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List BACKTRACK_cpp(Rcpp::IntegerMatrix dm){//direction matrix with new empty columns
              
   int n = dm.nrow();
   int m = dm.ncol();
   int i = n;
   int j = m;
   int step;
   vector<int> ii;
   vector<int> jj;
   vector<int> wp;
   
  ii.push_back(i);
  jj.push_back(j);
   
   do{
      step = dm(i-1,j-1);
      if(step == 1){
         i = i - 1;
         j = j - 1;
      } else if ( step == 2){
         j = j - 1;
      } else if ( step == 3){
         i = i - 1;
      } else{
         i = 99;
         j = 99;
      }
      ii.push_back(i);
      jj.push_back(j);
      wp.push_back(step);
   } while (i > 1 || j > 1);
   
   // ii.push_back(1);
   // jj.push_back(1);
   
   List ret;
   ret["ii"] = ii;
   ret["wp"] = wp;
   ret["jj"] = jj;
   return ret ;
}



// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List BACKTRACK2IN_cpp(Rcpp::IntegerMatrix dm, Rcpp::NumericMatrix diffM){
   
   int n = dm.nrow();
   int m = dm.ncol();
   int i = n;
   int j = m;
   int step;
   vector<int> ii;
   vector<int> jj;
   vector<int> wp;//warping path
   vector<double> diffp;//path of differences   
   diffp.push_back(diffM(n-1, m-1));
   
   ii.push_back(i);
   jj.push_back(j);
   do{
      step = dm(i-1,j-1);
      if(step == 1){
         i = i - 1;
         j = j - 1;
      } else if ( step == 2){
         j = j - 1;
      } else if ( step == 3){
         i = i - 1;
      } else{
         i = 99;
         j = 99;
      }
      ii.push_back(i);
      jj.push_back(j);
      wp.push_back(step);
      diffp.push_back(diffM(i-1, j-1));
   } while (i > 1 || j > 1);
   
   List ret;
   ret["ii"] = ii;
   ret["jj"] = jj;
   ret["wp"] = wp;
   ret["diffp"] = diffp;
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List BACKTRACK2II_cpp(Rcpp::IntegerMatrix dm, Rcpp::IntegerMatrix diffM){
   
   int n = dm.nrow();
   int m = dm.ncol();
   int i = n;
   int j = m;
   int step;
   vector<int> ii;
   vector<int> jj;
   vector<int> wp;//warping path
   vector<int> diffp;//path of differences   
   diffp.push_back(diffM(n-1, m-1));
   
   ii.push_back(i);
   jj.push_back(j);
   do{
      step = dm(i-1,j-1);
      if(step == 1){
         i = i - 1;
         j = j - 1;
      } else if ( step == 2){
         j = j - 1;
      } else if ( step == 3){
         i = i - 1;
      } else{
         i = 99;
         j = 99;
      }
      ii.push_back(i);
      jj.push_back(j);
      wp.push_back(step);
      diffp.push_back(diffM(i-1, j-1));
   } while (i > 1 || j > 1);
   
   List ret;
   ret["ii"] = ii;
   ret["jj"] = jj;
   ret["wp"] = wp;
   ret["diffp"] = diffp;
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
Rcpp::NumericMatrix normmat(Rcpp::NumericMatrix x){
   
   int n = x.nrow();
   int m = x.ncol();
   int i,j;
   for(i = 0; i < n; i++){
      for(j = 0; j < m; j++){
         x(i, j) = x(i, j)/(2 + i + j);
      }
   }
   return x ;
}


