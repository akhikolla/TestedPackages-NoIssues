#include <RcppArmadillo.h>
#include "IncDTW.h"
using namespace Rcpp;
using namespace std;




// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double get_lb(const NumericMatrix &tube, const NumericVector &x, 
              int j0, int jsup){
   double lb = 0;
   int k;
   
   for(int j = j0; j < jsup; j++){
      k = j - j0;
      if(x[j] > tube(k, 1)){
         
         lb += x[j] - tube(k, 1);
         
      }else if(x[j] < tube(k, 0)){
         
         lb += tube(k, 0)- x[j];
         
      }
      
   }
   return lb;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double get_lb_znorm(const NumericMatrix &tube, const NumericVector &x, 
                    double mu, double sigma, double threshold,
                    int j0, int jsup){
   double lb = 0;
   int k;
   double z;
   int j = j0;
   
   while((j < jsup) && (lb < threshold)){
      
      k = j - j0;
      z = (x[j] - mu)/sigma;
      
      if(z > tube(k, 1)){
         
         lb += z - tube(k, 1);
         
      }else if(z < tube(k, 0)){
         
         lb += tube(k, 0)- z;
         
      }
      j += 1;
   }
   
   return lb;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double get_lb_mv1(const NumericMatrix &tube, const NumericMatrix &x, 
                  int j0, int jsup, int nc){
   double lb = 0;
   int j, c, k;
   for(j = j0; j < jsup; j++){
      k = j - j0;
      for(c = 0; c < nc; c++){
         if(x(j,c) > tube(k, 2*c+1)){
            
            lb += x(j,c) - tube(k, 2*c+1);
            
         }else if(x(j,c) < tube(k, 2*c)){
            
            lb += tube(k, 2*c)- x(j,c);
            
         }
      }
   }
   return lb;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double get_lb_mv2(const NumericMatrix &tube, const NumericMatrix &x, 
                  int j0, int jsup, int nc){
   double lb = 0;
   double tmp, dimsum;
   int j, c, k;
   for(j = j0; j < jsup; j++){
      dimsum = 0;
      k = j-j0;
      
      for(c = 0; c < nc; c++){
         if(x(j,c) > tube(k, 2*c+1)){
            tmp = x(j,c) - tube(k, 2*c+1);
            dimsum += tmp * tmp;
            
         }else if(x(j,c) < tube(k, 2*c)){
            tmp = tube(k, 2*c)- x(j,c);
            dimsum += tmp * tmp;
            
         }
      }
      lb += sqrt(dimsum);
   }
   return lb;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double get_lb_mv22(const NumericMatrix &tube, const NumericMatrix &x, 
                   int j0, int jsup, int nc){
   double lb = 0;
   double tmp;
   int j, c, k;
   
   for(j = j0; j < jsup; j++){
      k = j - j0;
      for(c = 0; c < nc; c++){
         if(x(j,c) > tube(k, 2*c+1)){
            tmp = x(j,c) - tube(k, 2*c+1);
            lb += tmp * tmp;
            
         }else if(x(j,c) < tube(k, 2*c)){
            tmp = tube(k, 2*c)- x(j,c);
            lb += tmp * tmp;
            
         }
      }
   }
   return lb;
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


double get_lb_mv1_z(const NumericMatrix &tube, 
                    const NumericMatrix &x, 
                    const NumericVector &mu, 
                    const NumericVector &sigma, 
                    double threshold,
                    int j0, int jsup, int nc){
   
   double lb = 0;
   int c, k;
   double z;
   int j = j0;
   
   while((j < jsup) && (lb < threshold)){
      k = j - j0;
      
      for(c = 0; c < nc; c++){
         z = (x(j, c) - mu[c])/sigma[c];
         
         if(z > tube(k, 2*c+1)){
            
            lb += z - tube(k, 2*c+1);
            
         }else if(z < tube(k, 2*c)){
            
            lb += tube(k, 2*c) - z;
            
         }
      }
      j += 1;
   }
   
   return lb;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


double get_lb_mv2_z(const NumericMatrix &tube, 
                    const NumericMatrix &x, 
                    const NumericVector &mu, 
                    const NumericVector &sigma, 
                    double threshold,
                    int j0, int jsup, int nc){
   double lb = 0;
   double tmp, dimsum;
   int  c, k;
   double z;
   int j = j0;
   
   while((j < jsup) && (lb < threshold)){
      dimsum = 0;
      k = j - j0;
      
      for(c = 0; c < nc; c++){
         z = (x(j, c) - mu[c])/sigma[c];
         
         if(z > tube(k, 2*c+1)){
            
            tmp = z - tube(k, 2*c+1);
            dimsum += tmp * tmp;
            
         }else if(z < tube(k, 2*c)){
            
            tmp = tube(k, 2*c)- z;
            dimsum += tmp * tmp;
            
         }
      }
      lb += sqrt(dimsum);
      j += 1;
   }
   return lb;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


double get_lb_mv22_z(const NumericMatrix &tube, 
                     const NumericMatrix &x, 
                     const NumericVector &mu, 
                     const NumericVector &sigma, 
                     double threshold,
                     int j0, int jsup, int nc){
   double lb = 0;
   double tmp;
   int c, k;
   double z;
   int j = j0;
   
   while((j < jsup) && (lb < threshold)){
      
      k = j - j0;
      
      for(c = 0; c < nc; c++){
         z = (x(j, c) - mu[c])/sigma[c];
         
         if(z > tube(k, 2*c+1)){
            
            tmp = z - tube(k, 2*c+1);
            lb += tmp * tmp;
            
         }else if(z < tube(k, 2*c)){
            
            tmp = tube(k, 2*c)- z;
            lb += tmp * tmp;
            
         }
      }
      
      j += 1;
   }
   return lb;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

void cpp_cm(NumericMatrix &cm,
            const NumericVector &y, const NumericVector &h,
            int i0, int i1, int nh){
   // cpp_norm01_noreturn
   
   for(int i = 0; i < nh ; i++){
      for(int j = i0; j < i1; j++){
         cm(i, j) = abs(y[j] - h[i]);
      }
   }
   
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


void cpp_cm1_mv(NumericMatrix &cm,
                const NumericMatrix &y, const NumericMatrix &h,
                int j0, int j1, int nh, int nc){
   
   int c, i, j;
   double tmp;
   
   for(i = 0; i < nh ; i++){
      for(j = j0; j < j1; j++){
         tmp = 0;
         for(c = 0; c < nc; c++){
            tmp += abs(y(j,c) - h(i,c));
         }
         cm(i, j) = tmp;
      }
   }
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


void cpp_cm2_mv(NumericMatrix &cm,
                const NumericMatrix &y, const NumericMatrix &h,
                int j0, int j1, int nh, int nc){
   
   int c, i, j;
   double tmp_sum, tmp2;
   
   for(i = 0; i < nh ; i++){
      for(j = j0; j < j1; j++){
         tmp_sum = 0;
         for(c = 0; c < nc; c++){
            tmp2 = (y(j,c) - h(i,c));
            tmp_sum += tmp2 * tmp2;
         }
         cm(i, j) = sqrt(tmp_sum);
      }
   }
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


void cpp_cm2square_mv(NumericMatrix &cm,
                      const NumericMatrix &y, const NumericMatrix &h,
                      int j0, int j1, int nh, int nc){
   
   int c, i, j;
   double tmp_sum, tmp2;
   
   for(i = 0; i < nh ; i++){
      for(j = j0; j < j1; j++){
         tmp_sum = 0;
         for(c = 0; c < nc; c++){
            tmp2 = (y(j,c) - h(i,c));
            tmp_sum += tmp2 * tmp2;
         }
         cm(i, j) = tmp_sum;
      }
   }
}



// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double dist1_mv_z(const NumericMatrix &h, 
                  const NumericMatrix &x, 
                  const NumericVector &mu, 
                  const NumericVector &sigma, 
                  int ih, int ix, int nc){
   double ret = 0;
   for(int k = 0; k < nc; k++){
      ret += abs( (x(ix, k) - mu[k])/sigma(k) - h(ih, k) );
   }
   return (ret);
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double dist2_mv_z(const NumericMatrix &h, 
                  const NumericMatrix &x, 
                  const NumericVector &mu, 
                  const NumericVector &sigma, 
                  int ih, int ix, int nc){
   double ret = 0;
   double tmp = 0;
   for(int k = 0; k < nc; k++){
      tmp = (x(ix, k) - mu[k])/sigma(k) - h(ih, k) ;
      ret += tmp * tmp;
   }
   ret = sqrt(ret);
   return (ret);
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double dist22_mv_z(const NumericMatrix &h, 
                   const NumericMatrix &x, 
                   const NumericVector &mu, 
                   const NumericVector &sigma, 
                   int ih, int ix, int nc){
   double ret = 0;
   double tmp = 0;
   for(int k =0; k < nc; k++){
      tmp = (x(ix, k) - mu[k])/sigma(k) - h(ih, k) ;
      ret += tmp * tmp;
   }
   ret = ret * 1;
   return (ret);
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


typedef double (*funcPtr_step_run)(const double gcm10, 
                const double gcm11, 
                const double gcm01, 
                const double cm00);

XPtr<funcPtr_step_run> selectVecStepRun(std::string step_pattern) {
   if (step_pattern == "symmetric1")
      return(XPtr<funcPtr_step_run>(new funcPtr_step_run(&mystep_symmetric1)));
   else if (step_pattern == "symmetric2")
      return(XPtr<funcPtr_step_run>(new funcPtr_step_run(&mystep_symmetric2)));
   else
      return XPtr<funcPtr_step_run>(R_NilValue); // runtime error as NULL no XPtr
}



// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


typedef void (*funcPtr_cm)(NumericMatrix &cm,
              const NumericMatrix &y, 
              const NumericMatrix &h,
              int j0, int j1, int nh, int nc);

XPtr<funcPtr_cm> select_cm(std::string dist_method) {
   if (dist_method == "norm1")
      return(XPtr<funcPtr_cm>(new funcPtr_cm(&cpp_cm1_mv)));
   else if (dist_method == "norm2")
      return(XPtr<funcPtr_cm>(new funcPtr_cm(&cpp_cm2_mv)));
   else if (dist_method == "norm2_square")
      return(XPtr<funcPtr_cm>(new funcPtr_cm(&cpp_cm2square_mv)));
   else
      return XPtr<funcPtr_cm>(R_NilValue); // runtime error as NULL no XPtr
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


typedef double (*funcPtr_lb)(const NumericMatrix &tube, 
                const NumericMatrix &x, 
                int j0, int jsup, int nc);

XPtr<funcPtr_lb> select_lb(std::string dist_method) {
   if (dist_method == "norm1")
      return(XPtr<funcPtr_lb>(new funcPtr_lb(&get_lb_mv1)));
   else if (dist_method == "norm2")
      return(XPtr<funcPtr_lb>(new funcPtr_lb(&get_lb_mv2)));
   else if (dist_method == "norm2_square")
      return(XPtr<funcPtr_lb>(new funcPtr_lb(&get_lb_mv22)));
   else
      return XPtr<funcPtr_lb>(R_NilValue); // runtime error as NULL no XPtr
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


typedef double (*funcPtr_dist_mv_z)(const NumericMatrix &h, 
                const NumericMatrix &x, 
                const NumericVector &mu, 
                const NumericVector &sigma, 
                int ih, int ix, int nc);

XPtr<funcPtr_dist_mv_z> select_dist_mv_z(std::string dist_method) {
   if (dist_method == "norm1")
      return(XPtr<funcPtr_dist_mv_z>(new funcPtr_dist_mv_z(&dist1_mv_z)));
   else if (dist_method == "norm2")
      return(XPtr<funcPtr_dist_mv_z>(new funcPtr_dist_mv_z(&dist2_mv_z)));
   else if (dist_method == "norm2_square")
      return(XPtr<funcPtr_dist_mv_z>(new funcPtr_dist_mv_z(&dist22_mv_z)));
   else
      return XPtr<funcPtr_dist_mv_z>(R_NilValue); // runtime error as NULL no XPtr
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


typedef double (*funcPtr_lb_z)(const NumericMatrix &tube, 
                const NumericMatrix &x, 
                const NumericVector &mu, 
                const NumericVector &sigma, 
                double threshold,
                int j0, int jsup, int nc);

XPtr<funcPtr_lb_z> select_lb_z(std::string dist_method) {
   if (dist_method == "norm1")
      return(XPtr<funcPtr_lb_z>(new funcPtr_lb_z(&get_lb_mv1_z)));
   else if (dist_method == "norm2")
      return(XPtr<funcPtr_lb_z>(new funcPtr_lb_z(&get_lb_mv2_z)));
   else if (dist_method == "norm2_square")
      return(XPtr<funcPtr_lb_z>(new funcPtr_lb_z(&get_lb_mv22_z)));
   else
      return XPtr<funcPtr_lb_z>(R_NilValue); // runtime error as NULL no XPtr
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

NumericVector cpp_range(NumericVector x, int i0, int i1){
   NumericVector range(2);
   range[0] = x[i0];
   range[1] = x[i0];
   
   for(int i = i0+1; i < i1; i++){
      // min
      if(x[i] < range[0]){
         range[0] = x[i];
      }
      // max
      if(x[i] > range[1]){
         range[1] = x[i];
      }
   }
   return range;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


void cpp_set_range_mv(NumericMatrix & range,
                      const NumericMatrix &x, int c, int i0, int i1){
   
   double xmin = x(i0, c);
   double xmax = x(i0, c);
   
   for(int i = i0+1; i < i1; i++){
      // min
      if(x(i,c) < xmin){
         xmin = x(i,c);
      }
      // max
      if(x(i,c) > xmax){
         xmax = x(i,c);
      }
   }
   range(0,c) = xmin;
   range(1,c) = xmax;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
NumericMatrix cpp_get_tube( const NumericVector &h, int ws){
   // LU ... matrix of the dimension nx h 2
   
   int i, j, j0, j1;
   int nh = h.size();
   double tmp_min;
   double tmp_max;
   NumericMatrix tube(nh, 2);
   
   for(i = 0; i < nh; i++){
      j0 = max(0, i - ws);
      j1 = min(nh, i + ws);
      // Rcout<< "i: " << i <<" j0: "<< j0<< " j1: "<<j1 <<"\n";
      tmp_min = h[j0];
      tmp_max = h[j0];
      for(j = j0+1; j < j1; j++){
         if(h[j] < tmp_min ){
            tmp_min = h[j];
         }
         if(h[j] > tmp_max ){
            tmp_max = h[j];
         }
      }
      tube(i,0) = tmp_min;
      tube(i,1) = tmp_max;
   }
   
   return tube;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
void cpp_set_tube(NumericMatrix &tube, const NumericVector &h, int ws){
   // LU ... matrix of the dimension nx h 2
   int i, j, j0, j1;
   int nh = h.size();
   double tmp_min;
   double tmp_max;
   
   
   for(i = 0; i < nh; i++){
      j0 = max(0, i - ws);
      j1 = min(nh, i + ws);
      tmp_min = h[j0];
      tmp_max = h[j0];
      for(j = j0+1; j < j1; j++){
         if(h[j] < tmp_min ){
            tmp_min = h[j];
         }
         if(h[j] > tmp_max ){
            tmp_max = h[j];
         }
      }
      tube(i,0) = tmp_min;
      tube(i,1) = tmp_max;
   }
   
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


void cpp_set_tube_mv(NumericMatrix &tube, const NumericMatrix &h, int ws){
   // LU ... matrix of the dimension nx h 2*nc
   int c,i, j, j0, j1;
   int nh = h.nrow();
   int nc = h.ncol();
   double tmp_min;
   double tmp_max;
   // NumericMatrix range(nh, 2*nc);
   
   for(i = 0; i < nh; i++){
      j0 = max(0, i - ws);
      j1 = min(nh, i + ws);
      
      for(c = 0; c < nc; c++){
         tmp_min = h(j0,c);
         tmp_max = h(j0,c);
         for(j = j0+1; j < j1; j++){
            if(h(j, c) < tmp_min ){
               tmp_min = h(j, c);
            }
            if(h(j, c) > tmp_max ){
               tmp_max = h(j, c);
            }
         }
         tube(i,2*c)   = tmp_min;
         tube(i,2*c+1) = tmp_max;
      }
   }
   
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



// [[Rcpp::export]]
NumericMatrix cpp_get_tube_mv(const NumericMatrix &h, int ws){
   // LU ... matrix of the dimension nx h 2*nc
   int c,i, j, j0, j1;
   int nh = h.nrow();
   int nc = h.ncol();
   double tmp_min;
   double tmp_max;
   // NumericMatrix range(nh, 2*nc);
   NumericMatrix tube(nh, 2*nc); 
   
   for(i = 0; i < nh; i++){
      j0 = max(0, i - ws);
      j1 = min(nh, i + ws);
      
      for(c = 0; c < nc; c++){
         tmp_min = h(j0,c);
         tmp_max = h(j0,c);
         for(j = j0+1; j < j1; j++){
            if(h(j, c) < tmp_min ){
               tmp_min = h(j, c);
            }
            if(h(j, c) > tmp_max ){
               tmp_max = h(j, c);
            }
         }
         tube(i,2*c)   = tmp_min;
         tube(i,2*c+1) = tmp_max;
      }
   }
   return tube;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double get_mean(const NumericVector &x, int i0, int imax){
   double mu = 0;
   double nh = (imax - i0 + 1);
   for(int i = i0; i <= imax; i++){
      mu += x[i];
   }
   mu /= nh;
   return mu;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


double get_sigma2(const NumericVector &x, int i0, int imax, double mu2){
   double sigma2 = 0;
   double nh = (imax - i0 + 1);
   for(int i = i0; i <= imax; i++){
      sigma2 += x[i] * x[i];
   }
   sigma2 = sigma2/(nh-1) - mu2 * nh/(nh-1);
   return sigma2;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


void set_mean_sigma_mv(NumericVector &mu,
                       NumericVector &sigma,
                       NumericVector &sigma2,
                       const NumericMatrix &x, 
                       int i0, int imax, int nc){
   
   double nh = (imax - i0 + 1);
   double tmp = 0;
   double mu2 = 0;
   for(int c = 0; c < nc; c++){
      tmp = 0;
      // set mean
      for(int i = i0; i <= imax; i++){
         tmp += x(i, c);
      }
      mu[c] = tmp/nh;
      mu2 = mu[c] * mu[c];
      
      // set sigma
      tmp = 0;
      for(int i = i0; i <= imax; i++){
         tmp += x(i,c) * x(i,c);
      }
      sigma2[c] = tmp/(nh-1) - mu2 * nh/(nh-1);
      
      // if(sigma2[c] < var_precision){
      if(sigma2[c] < 0.0000000001){
         sigma[c] = 1; 
      } else{
         sigma[c] = sqrt(sigma2[c]);
      }
      
   }
}




//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

void cpp_norm01(const NumericVector &x, NumericVector &y,
                int i0, int i1, double xmin, double xmax){
   // cpp_norm01_noreturn
   
   double nominator = xmax - xmin;
   
   // double sd_precision = 0.000000001; // 1e-9
   // if(nominator < sd_precision) nominator = 1;
   if(nominator < 0.000000001) nominator = 1;
   
   for(int i = i0; i < i1; i++){
      y[i] = (x[i] - xmin)/nominator ;
   }
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


void cpp_norm01_mv(const NumericMatrix &x, NumericMatrix &y,
                   int c, int i0, int i1, double xmin, double xmax){
   
   
   double nominator = xmax - xmin;
   
   // double sd_precision = 0.000000001; // 1e-9
   // if(nominator < sd_precision) nominator = 1;
   if(nominator < 0.000000001) nominator = 1;
   
   for(int i = i0; i < i1; i++){
      y(i,c) = (x(i,c) - xmin)/nominator ;
   }
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

struct kNN_Info {
   double vmax; // max(kNN_val)
   int which_vmax; // which(kNN_val == max(kNN_val))
   int imax;// max(kNN_ix)
   int which_imax; //index of kNN_ix which is == imax, so: which(kNN_imax == imax)
   int nr_detected; //number of detected fits
   int nr_looking4; //number of fits we are looking for
} ;

void debug_print_kNN(int debug, std::string position, int j, kNN_Info &ki, 
                     NumericVector &kNN_val, IntegerVector &kNN_ix){
   if(debug == 1){
      Rcout << "j = "<< j << ", " << position << ": ";
      Rcout << " imax: "<< ki.imax<< " "<< 
         " which_imax: "<<ki.which_imax << " "<< 
            " vmax: "<<ki.vmax << " "<< 
               " which_vmax: "<<ki.which_vmax << " "<< 
                  " nr_detected: "<<ki.nr_detected <<" "<< 
                     " nr_looking4: "<<ki.nr_looking4 <<"\n "; 
      Rcout << kNN_val << "\n";
      Rcout << kNN_ix << "\n";
   }
}


kNN_Info fill_kNN_Info(const List &kNN_inf_list){
   
   // this function gets only called if (kNN_inf_list["nr_looking4"] > 0)
   
   kNN_Info kNN_inf;
   // as<int>(kNN_inf_list["nr_looking4"]);
   kNN_inf.vmax        = as<double>(kNN_inf_list["vmax"]);      // default = R_PosInf
   kNN_inf.which_vmax  = as<int>(kNN_inf_list["which_vmax"]);   // default = 0
   kNN_inf.which_imax  = as<int>(kNN_inf_list["which_imax"]);   // default = 0
   kNN_inf.imax        = as<int>(kNN_inf_list["imax"]);         // default = 0
   kNN_inf.nr_detected = as<int>(kNN_inf_list["nr_detected"]);  // default = 1
   kNN_inf.nr_looking4 = as<int>(kNN_inf_list["nr_looking4"]);  // default = kNNk
   
   if(kNN_inf.nr_detected < kNN_inf.nr_looking4){
      kNN_inf.which_vmax = -99;   
   }
   
   return kNN_inf;
}


void kick_vmax_kNN(NumericVector &kNN_val, IntegerVector &kNN_ix, 
                   kNN_Info &ki, double new_val, int new_ix){
   // replace the k-th nearest neighbor, so the one with the maximum distance
   // of all k nearest neighbors
   
   int j=0;
   if(ki.nr_detected < ki.nr_looking4){
      for(int i = 0; i < kNN_val.size(); i++){
         if (kNN_ix[i] == -99) {
            j = i;  
            break;
         }
      }
      
      kNN_val[j] = new_val;
      kNN_ix[j] = new_ix;
      ki.imax = new_ix;
      ki.which_imax = j;
      
      ki.nr_detected += 1;
      
   }else{
      kNN_val[ki.which_vmax] = new_val;
      kNN_ix[ki.which_vmax] = new_ix;
      ki.imax = new_ix;
      ki.which_imax = ki.which_vmax;
   }
   
   // update vmax info
   if(ki.nr_detected == ki.nr_looking4){
      double bsf = kNN_val[0];
      int new_i = 0;
      for(int i = 1; i< kNN_val.size(); i++){
         if(kNN_val[i] > bsf){
            bsf = kNN_val[i];
            new_i = i;
         }
      }
      ki.vmax = bsf;
      ki.which_vmax = new_i;
   }else{
      
      ki.vmax = R_PosInf;
      ki.which_vmax = -99;
   }
   
}


void kick_imax_kNN(NumericVector &kNN_val, IntegerVector &kNN_ix, 
                   kNN_Info &ki, double new_val, int new_ix){
   // replace the most recent nearest neighbor, so the one with 
   // the biggest index
   
   kNN_val[ki.which_imax] = new_val;
   kNN_ix[ki.which_imax] = new_ix;
   ki.imax = new_ix;
   // ki.which_imax remains the same
   
   // update vmax info
   if(ki.nr_detected == ki.nr_looking4){
      if(new_val > ki.vmax){
         ki.vmax = new_val;
         ki.which_vmax = ki.which_imax ;
         
      }else{
         double bsf = kNN_val[0];
         int new_i = 0;
         for(int i = 1; i< kNN_val.size(); i++){
            if(kNN_val[i] > bsf){
               bsf = kNN_val[i];
               new_i = i;
            }
         }
         ki.vmax = bsf;
         ki.which_vmax = new_i;
      }
      
      // }else if(ki.which_imax == ki.which_vmax){
      //    ki.vmax = new_val;
      // }
   }
}


void kick_vmax_kNN_lot(NumericVector &kNN_val, IntegerVector &kNN_ix, 
                       IntegerVector &kNN_lot_ix, kNN_Info &ki, double new_val, 
                       int new_ix, int lot_ix){
   // replace the k-th nearest neighbor, so the one with the maximum distance
   // of all k nearest neighbors
   
   int j=0;
   if(ki.nr_detected < ki.nr_looking4){
      for(int i = 0; i < kNN_val.size(); i++){
         if (kNN_ix[i] == -99) {
            j = i;  
            break;
         }
      }
      
      kNN_val[j] = new_val;
      kNN_ix[j] = new_ix;
      kNN_lot_ix[j] = lot_ix;
      ki.imax = new_ix;
      ki.which_imax = j;
      
      ki.nr_detected += 1;
      
   }else{
      kNN_val[ki.which_vmax] = new_val;
      kNN_ix[ki.which_vmax] = new_ix;
      kNN_lot_ix[ki.which_vmax] = lot_ix;
      ki.imax = new_ix;
      ki.which_imax = ki.which_vmax;
   }
   
   // update vmax info
   if(ki.nr_detected == ki.nr_looking4){
      double bsf = kNN_val[0];
      int new_i = 0;
      for(int i = 1; i< kNN_val.size(); i++){
         if(kNN_val[i] > bsf){
            bsf = kNN_val[i];
            new_i = i;
         }
      }
      ki.vmax = bsf;
      ki.which_vmax = new_i;
      // Rcout << "\n new ki.vmax set: " << ki.vmax << " _ " << new_i;
   }else{
      
      ki.vmax = R_PosInf;
      ki.which_vmax = -99;
   }
   
}


void kick_imax_kNN_lot(NumericVector &kNN_val, IntegerVector &kNN_ix,
                       IntegerVector &kNN_lot_ix, kNN_Info &ki, double new_val, 
                       int new_ix, int lot_ix){
   // replace the most recent nearest neighbor, so the one with 
   // the biggest index
   
   kNN_val[ki.which_imax] = new_val;
   kNN_ix[ki.which_imax] = new_ix;
   kNN_lot_ix[ki.which_imax] = lot_ix;
   ki.imax = new_ix;
   // ki.which_imax remains the same
   
   // update vmax info
   if(ki.nr_detected == ki.nr_looking4){
      if(new_val > ki.vmax){
         ki.vmax = new_val;
         ki.which_vmax = ki.which_imax ;
         
      }else{
         double bsf = kNN_val[0];
         int new_i = 0;
         for(int i = 1; i< kNN_val.size(); i++){
            if(kNN_val[i] > bsf){
               bsf = kNN_val[i];
               new_i = i;
            }
         }
         ki.vmax = bsf;
         ki.which_vmax = new_i;
      }
      
      // }else if(ki.which_imax == ki.which_vmax){
      //    ki.vmax = new_val;
      // }
   }
}



void initialize_kNN(kNN_Info &ki, NumericVector &kNN_val, 
                 IntegerVector &kNN_ix, IntegerVector &kNN_lot_ix,
                 int lot_ix, int kNNk, double initial_bsfiw, int overlap_size){
     
   if(lot_ix == 1){
      // for the first time series of lot set kNN info this way, 
      
      kNN_val[0] = initial_bsfiw;
      kNN_ix[0] = 0;
      kNN_lot_ix[0] = lot_ix;
      
      for(int k = 1; k < kNNk; k++){
         kNN_val[k] = R_PosInf;
         kNN_ix[k] = -99;
      }
   }else{
      
      // artificially set imax
      ki.imax = -overlap_size - 1;
      
      if(ki.nr_detected < ki.nr_looking4){
         kick_vmax_kNN_lot(kNN_val, kNN_ix, kNN_lot_ix, 
                           ki, initial_bsfiw, 0, lot_ix);
         
      }else{
         if(initial_bsfiw < ki.vmax){
            kick_vmax_kNN_lot(kNN_val, kNN_ix, kNN_lot_ix, 
                              ki, initial_bsfiw, 0, lot_ix);
         }
      }
   }
  
}


void update_kNN(NumericVector &kNN_val, IntegerVector &kNN_ix, 
                double new_val, int new_ix){
   
   int n = kNN_val.size();
   int k = (int) n/2;
   int k_prev=0;
   int i0 = 0;
   int i1 = n;
   
   if(new_val <= kNN_val[0]){
      k = 0;
   }else{
      while (k != k_prev){
         k_prev = k;
         if(new_val < kNN_val[k]){
            i1 = k;
            k = (int) i0 + (i1-i0)/2;
         }else{
            i0 = k;
            k = (int) i0 + (i1-i0)/2;
         }
      }
      k += 1;
   }
   
   
   
   for(int i = n-1; i>=k ; i--){
      kNN_val[i] = kNN_val[i-1];
      kNN_ix[i] = kNN_ix[i-1];
   }
   kNN_val[k] = new_val;
   kNN_ix[k] = new_ix;
   
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
NumericVector cpp_kNN_rev(const NumericVector &disvec, int w, int debug = 0){
   
   int nd = disvec.size();
   double bsfir = R_PosInf;// best sofar in range
   int i, imin;
   int bi = 0;
   std::vector<int> best_indices;
   
   imin = std::max(0, nd-w);
   
   // Rcout << ndisvec<< " "<< w<< "  "<< ndisvec-w<<"\n";
   for(i = nd-1; i > imin; i--){
      if(disvec[i] < bsfir){
         bsfir = disvec[i];
         bi = i;
         // Rcout << i << "  "<< bsfir<<"\n";
      }
   }
   
   if(debug == 1){
      Rcout << "i: "<<i << " ---- end of initial ---- \n";            
   } 
   
   for(i = imin; i >= 0 ; i--){
      if(debug == 1){
         Rcout << "i: "<<i << " bi: "<<bi << " bsfir: "<< bsfir <<"\n";            
      } 
      
      if(bi-i  >= w){
         if(debug == 1){
            Rcout << "i: "<<i << " ---- adding ---- \n";
            Rcout << "i: "<<i << " bi: "<<bi << " bsfir: "<< bsfir <<"\n";            
            // Rcout << "i: "<<i << "adding... bi: "<<bi << "... bsfir: "<< bsfir <<"\n";
         }
         best_indices.push_back(bi);
         if(disvec[i] == disvec[i]){
            bsfir = disvec[i];
         }else{
            // disvec[i] is NaN
            bsfir = R_PosInf;
         }
         bi = i;
         
      }else if(disvec[i] < bsfir){
         if(debug == 1){
            Rcout << "i: "<<i << " ---- updating ---- \n";
            Rcout << "i: "<<i << " bi: "<<bi << " bsfir: "<< bsfir <<"\n";            
            
            // Rcout << "i: "<<i << "updating... bi: "<<bi << "... bsfir: "<< bsfir <<"\n";
         }
         
         bsfir = disvec[i];
         bi = i;
      }
   }
   // Rcout << "i: "<<i << "... bi: "<<bi << "... disvec[bi]: "<<disvec[bi]<<"\n";
   // sumsavings += bsfir;
   
   if(disvec[bi] == disvec[bi]){
      best_indices.push_back(bi);
   }
   
   return wrap( best_indices );
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
IntegerVector cpp_local_min(const NumericVector &x, int w, int strict){
   int nx = x.size();
   int i, j, imin, imax;
   int counter=0;
   double xmax = x[0];
   std::vector<int> best_indices;
   
   // get maximum of x
   for(i = 1; i< nx; i++){
      if(x[i] > xmax){
         xmax = x[i];
      }
   }
   // +1 so that saddle points at the beginning are still detected
   xmax = xmax + 1;
   
   // create new vector y=x and replace all nan values with xmax + 1
   NumericVector y(nx);
   for(i = 0; i< nx; i++){
      if(x[i] != x[i]){
         y[i] = xmax;
      }else{
         y[i] = x[i];
      }
   }
   
   if(strict == 0){
      // counter adds 1 if y[j] < y[i]
      // and also      if y[j] == y[i]
      for( i = 0; i< nx; i++){
         if(x[i]!=x[i]){
            continue;
         }
         imin = std::max(0, i - w);
         imax = std::min(nx, i + w+1);
         counter = 0;
         
         for(j = imin; j < imax; j++){
            if(y[j] < y[i]){
               break;
            }
            counter += 1;
         }
         if(counter == imax - imin){
            best_indices.push_back(i+1);
         }
      }
      
   } else{
      // counter only adds 1 if y[j] < y[i]
      for( i = 0; i< nx; i++){
         if(x[i]!=x[i]){
            continue;
         }
         imin = std::max(0, i - w);
         imax = std::min(nx, i + w+1);
         counter = 0;
         
         for(j = imin; j < imax; j++){
            if(j != i){
               if(y[j] <= y[i]){
                  break;
               }
            }
            counter += 1;
         }
         if(counter == imax - imin){
            best_indices.push_back(i+1);
         }
      }
   }
   return wrap( best_indices );
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_rundtw(const NumericVector &h,
                const NumericVector &x,
                std::string step_pattern,
                List &kNN_inf_list,
                int ws, double threshold, int overlap_tol = 0, 
                int do_norm = 1, int use_ea = 1, int use_lb = 1, int debug = 0) {
   // scan x with h and get the diff of DL(ss(x)) - DL(ss(x)-h) for each
   // subsegment ss(x)
   // h ... normalized
   // x ... not normalized
   // ws ... window size parameter for sakoe chiba window
   // threshold ... global threshold for early abandoning dtw calc
   //             (if threshold = Inf, then no thresholding takes place,
   //              in that case use_ea should be 0)
   // do_norm, use_ea, use_lb ... either 0 or 1
   
   
   // set step function
   SEXP step_SEXP = selectVecStepRun(step_pattern);
   XPtr<funcPtr_step_run> xpfun_step(step_SEXP);
   funcPtr_step_run step_fun = *xpfun_step;
   
   if(use_lb == 1){
      use_ea = 1; // lb makes no sense if the loc_threshold is not set by bsfiw
   }
   
   if(use_ea == 0){
      threshold = R_PosInf;
   }
   
   int pastcheck = 0;
   int nh = h.size();
   int overlap_size = nh - overlap_tol;
   int nx = x.size();
   int j, k, jmax, jsup, bsfiw_ix;
   int counter_norm_full = 0; // ... 1 for the initial step
   int counter_norm_1step = 0;
   int counter_norm_newEx = 0;
   int counter_cm_full = 0;
   int counter_cm_1step = 0;
   int counter_lb = 0;
   int counter_ea = 0;
   
   int norm_exit_status = 0;
   int last_cm_calc_j = -1;
   
   double bsfiw, loc_threshold; // local threshold
   
   int kNNk = as<int>(kNN_inf_list["nr_looking4"]);
   IntegerVector kNN_ix (kNNk);
   NumericVector kNN_val (kNNk, 0.0);
   kNN_Info kNN_inf;
   
   NumericVector y (nx, 0.0); // normalized version of x
   NumericVector range(2); // range vector for running min-max
   NumericMatrix cm(nh, nx); //cost matrix
   NumericVector ret(nx-nh+1, 0.0); // return vetor
   
   
   // set initial maximum indices
   j = 0;
   jsup = j + nh;
   jmax = jsup -1;
   
   if(do_norm == 1){
      // get rmin and rmax for the start
      range = cpp_range(x, j, jsup);//cpp_range(x, 0, nh); // range = c(min, max)
      
      // normalize the first nh values of x and fill y
      cpp_norm01(x, y, j, nh, range[0], range[1]);
      counter_norm_full += 1;
      
   } else{
      y = x; 
      norm_exit_status = 4;
   }
   
   
   // get initial dtw value
   loc_threshold = threshold;
   double lb = 0;
   
   NumericMatrix tube(nh,2);
   if(use_lb == 1){
      cpp_set_tube(tube, h, ws);
      lb = get_lb(tube, y, j, jsup);
      
   }else{
      lb = -1;
   }
   // double get_lb(const NumericMatrix &tube, const NumericVector &x, 
   //               int j0, int jsup){
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW BEGIN      //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   int iBegin = 0;
   int iEnd = 0;
   int nanCounter_status = 1;
   
   double * p1 = new double [nh];
   double * p2 = new double [nh];
   double * ptmp;
   double mynan;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   
   if(loc_threshold <= lb){
      //  loc_threshold < lb < DTW 
      ret[j] = mynan;
      counter_lb += 1;
      
   }else{
      // fill the initial nh columns of the costmatrix
      cpp_cm(cm, y, h, j, jsup, nh);
      counter_cm_full += 1;
      last_cm_calc_j = j;
      
      
      //initialize a and b with NAN
      for(int i=0; i < nh; i++){
         p1[i] = mynan;
         p2[i] = mynan;
      }
      
      // first column
      p1[0] = cm(0, j);
      if(p1[0] > loc_threshold){
         ret[j] = mynan;
         counter_ea += 1;
         
      }else{
         
         iEnd   = std::min(nh, ws+1);
         for(int i=1; i < iEnd; i++){
            p1[i] = cm(i, j) + p1[i-1];
            if(p1[i] > loc_threshold) p1[i] = mynan;
         }
         
         for(int k = j +1 ; k < jsup; k++){
            nanCounter = 0;
            iBegin = k-ws-j;
            if(iBegin <= 0){
               *p2 = cm(0, k) + *(p1);
               if(*(p2) > loc_threshold){
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
            
            iEnd   = k+ws+1-j;
            if(iEnd >= nh){
               iEnd = nh;
            }else{
               *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
            }
            
            for (int i = iBegin; i < iEnd; i++){
               *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
               if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                  *(p2+i) = mynan;
                  nanCounter ++;
               }
            }
            
            if(nanCounter == nh) {
               nanCounter_status = 0;
               break;
            }
            ptmp=p1;
            p1 = p2;
            p2 = ptmp;
         }
         
         if(nanCounter_status == 1){
            ret[j] = *(p1+nh-1);//p1[nx-1]
         } else{
            ret[j] = mynan;
            counter_ea += 1;
         }
      }
      
   }// lower bounding if clause
   
   // Rcout << loc_threshold <<" ... "<< lb << " ... " << ret[j]<< "\n";
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW END        //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   // ret[0] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
   // ret[0] = cpp_dtw2vec_cm(cm, i, jsup, nh);
   loc_threshold = mymin(ret[0], threshold);
   
   // ret[0] can be NaN => in this case set bsfiw to Inf, for proper comparisons,
   // since comparisons with NaN are always FALSE
   double initial_bsfiw = mymin(ret[0], R_PosInf);
   
   // fill initially kNN_val
   if(kNNk > 0){
      kNN_val[0] = initial_bsfiw;
      kNN_ix[0] = 0;
      
      for(k = 1; k < kNNk; k++){
         kNN_val[k] = R_PosInf;
         kNN_ix[k] = -99;
      }
      
      kNN_inf = fill_kNN_Info(kNN_inf_list);
   }
   
   bsfiw = initial_bsfiw; // best so far in window
   bsfiw_ix = 0; // last index of best so far in window
   
   
   // --------------------------------------start rolling----------------------------
   for(j=1; j < (nx-nh+1); j++){
      if((j - bsfiw_ix) >= overlap_size){
         bsfiw = threshold;
         bsfiw_ix = j;
      }
      
      jsup = j + nh; // the frontrunner index +1
      jmax = jsup -1; // the index of the frontrunner
      
      // update the discretization of the new window
      // onlx_frontrunner_updated = 0;
      if(do_norm == 1){
         if(pastcheck == 1){
            if(x[jmax] < range[0]){
               // set new rmin and update whole xd
               range[0] = x[jmax];
               cpp_norm01(x, y, j, jsup, range[0], range[1]);
               counter_norm_newEx += 1;
               norm_exit_status = 1;
               
            }else if(x[jmax] > range[1]){
               // set new rmax and update whole xd
               range[1] = x[jmax];
               cpp_norm01(x, y, j, jsup, range[0], range[1]);
               counter_norm_newEx += 1;
               norm_exit_status = 1;
               
            }else{
               // only normalize y[jmax]
               y[jmax] = (x[jmax] - range[0])/(range[1] - range[0]) ;
               counter_norm_1step += 1;
               norm_exit_status = 4;
               
            }
         }else{
            //complete reset
            range = cpp_range(x, j, jsup); // range = c(min, max)
            cpp_norm01(x, y, j, jsup, range[0], range[1]);
            counter_norm_full += 1;
            norm_exit_status = 3;
            
         }
         
         // the next time i increases, x[i] is left behind, if x[i] was the
         // min or max, then new range needs to be defined
         if((x[j] == range[0]) || (x[j] == range[1])){
            pastcheck = 0;
         }else{
            pastcheck = 1;
         }
         
         
      }//end if do_norm
      
      
      //  Lower bounding
      if(use_lb == 1){
         lb = get_lb(tube, y, j, jsup);
      }else{
         lb = -1;
      }
      
      //  bsfiw serves as threshold
      if(use_ea == 1){
         loc_threshold = mymin(bsfiw, threshold);
      }else{
         loc_threshold = threshold;
      }
      
      
      if(kNNk > 0){
         // update local threshold, set it to the k-nearest neighbor
         // loc_threshold = mymin(loc_threshold, kNN_val[kNNk-1]);
         loc_threshold = mymin(loc_threshold, kNN_inf.vmax);
      }
      
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW BEGIN      //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      if(loc_threshold <= lb){
         //  loc_threshold < lb < DTW 
         ret[j] = mynan;
         counter_lb += 1;
         
      }else{
         
         if(norm_exit_status < 4 || last_cm_calc_j < j-1 ){
            // full reset of cm
            cpp_cm(cm, y, h, j, jsup, nh);
            counter_cm_full += 1;
         } else {
            // front running column of cm
            cpp_cm(cm, y, h, jmax, jsup, nh);
            counter_cm_1step += 1;
         }
         last_cm_calc_j = j;
         
         
         nanCounter_status = 1;
         
         //initialize a and b with NAN
         for(int i=0; i < nh; i++){
            p1[i] = mynan;
            p2[i] = mynan;
         }
         
         // first column
         p1[0] = cm(0, j);
         if(p1[0] > loc_threshold){
            ret[j] = mynan;
            counter_ea += 1;
            
         }else{
            
            iEnd   = std::min(nh, ws+1);
            for(int i=1; i < iEnd; i++){
               p1[i] = cm(i, j) + p1[i-1];
               if(p1[i] > loc_threshold) p1[i] = mynan;
            }
            
            for(int k = j +1 ; k < jsup; k++){
               nanCounter = 0;
               iBegin = k-ws-j;
               if(iBegin <= 0){
                  *p2 = cm(0, k) + *(p1);
                  if(*(p2) > loc_threshold){
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
               
               iEnd   = k+ws+1-j;
               if(iEnd >= nh){
                  iEnd = nh;
               }else{
                  *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
               }
               
               for (int i = iBegin; i < iEnd; i++){
                  *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
                  if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                     *(p2+i) = mynan;
                     nanCounter ++;
                  }
               }
               
               if(nanCounter == nh) {
                  nanCounter_status = 0;
                  break;
               }
               ptmp=p1;
               p1 = p2;
               p2 = ptmp;
            }
            
            if(nanCounter_status == 1){
               ret[j] = *(p1+nh-1);//p1[nx-1]
               // Rcout << ret[j] << " ";
            } else{
               ret[j] = mynan;
               counter_ea += 1;
            }
         }
         
         
      }// lower bound if clause
      
      // Rcout << loc_threshold <<" ... "<< lb << " ... " << ret[j]<< "\n";
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW END        //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      
      
      // ret[j] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
      // ret[i] = cpp_dtw2vec_cm(cm, j, jsup, nh);
      
      
      
      
      if(kNNk > 0){
         // update set of k nearest neighbors if new value is smaller than 
         // k-best value so far
         
         debug_print_kNN(debug, "before", j, kNN_inf, kNN_val, kNN_ix);
         
         
         
         if(ret[j] == ret[j]){ // if not nan
            if(j < kNN_inf.imax + overlap_size){
               // there is an overlap
               if(ret[j] < kNN_val[kNN_inf.which_imax]){
                  // because of overlap the new value needs to be smaller
                  // than the most recent best-sofar-value => so we replace
                  // this value with the new one
                  kick_imax_kNN(kNN_val, kNN_ix, kNN_inf, ret[j], j);
                  
                  if(debug == 1){
                     Rcout << "---kick imax---\n";
                  }
               }
            }else{
               // no overlap
               if(ret[j] < kNN_inf.vmax){
                  // new value needs to be at least smaller 
                  // than the farthest nearest neighbor
                  // update_kNN(kNN_val, kNN_ix, ret[j], j);   
                  kick_vmax_kNN(kNN_val, kNN_ix, kNN_inf, ret[j], j);
                  
                  if(debug == 1){
                     Rcout << "---kick vmax---\n";
                  }
               }
            }
         }
         
         debug_print_kNN(debug, "after", j, kNN_inf, kNN_val, kNN_ix);
      }
      
      
      if(ret[j] < bsfiw){
         bsfiw = ret[j];
         bsfiw_ix = j;
      }
      
   }//end rolling
   
   
   List list_ret;
   list_ret["dist"] = ret;
   IntegerVector counter(7);
   counter[0] = counter_norm_full;
   counter[1] = counter_norm_newEx;
   counter[2] = counter_norm_1step;
   counter[3] = counter_cm_full;
   counter[4] = counter_cm_1step;
   counter[5] = counter_ea;
   counter[6] = counter_lb;
   list_ret["counter"] = counter;
   
   if(kNNk > 0){
      // final reverse minimum search of ret:
      // up to here the kNN values and indices dealt as additional 
      // local threshold, but since the possibility of protracting
      // a local minimum, a final reverse check is necessary:
      
      list_ret["all_best_indices"] = cpp_kNN_rev(ret, overlap_size, debug);
   }
   delete[] p1;
   delete[] p2;
   
   return list_ret;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// [[Rcpp::export]]
List cpp_rundtw_lot(const NumericVector &h,
                    const NumericVector &x,
                    NumericVector kNN_val_in,
                    IntegerVector kNN_ix_in,
                    IntegerVector kNN_lot_ix_in,
                    List &kNN_inf_list_in, 
                    int lot_ix,
                    std::string step_pattern,
                    int ws, double threshold, int overlap_tol = 0, 
                    int do_norm = 1, int use_ea = 1, int use_lb = 1, int debug = 0) {
   // scan x with h and get the diff of DL(ss(x)) - DL(ss(x)-h) for each
   // subsegment ss(x)
   // h ... normalized
   // x ... not normalized
   // ws ... window size parameter for sakoe chiba window
   // threshold ... global threshold for early abandoning dtw calc
   //             (if threshold = Inf, then no thresholding takes place,
   //              in that case use_ea should be 0)
   // do_norm, use_ea, use_lb ... either 0 or 1
   
   
   NumericVector kNN_val    = clone(kNN_val_in);
   IntegerVector kNN_ix     = clone(kNN_ix_in);
   IntegerVector kNN_lot_ix = clone(kNN_lot_ix_in);
   List kNN_inf_list        = clone(kNN_inf_list_in);
   
   
   
   // set step function
   SEXP step_SEXP = selectVecStepRun(step_pattern);
   XPtr<funcPtr_step_run> xpfun_step(step_SEXP);
   funcPtr_step_run step_fun = *xpfun_step;
   
   if(use_lb == 1){
      use_ea = 1; // lb makes no sense if the loc_threshold is not set by bsfiw
   }
   
   if(use_ea == 0){
      threshold = R_PosInf;
   }
   
   int pastcheck = 0;
   int nh = h.size();
   int overlap_size = nh - overlap_tol;
   int nx = x.size();
   int j, jmax, jsup, bsfiw_ix;
   int counter_norm_full = 0; // ... 1 for the initial step
   int counter_norm_1step = 0;
   int counter_norm_newEx = 0;
   int counter_cm_full = 0;
   int counter_cm_1step = 0;
   int counter_lb = 0;
   int counter_ea = 0;
   
   int norm_exit_status = 0;
   int last_cm_calc_j = -1;
   
   double bsfiw, loc_threshold; // local threshold
   
   int kNNk = as<int>(kNN_inf_list["nr_looking4"]);
   kNN_Info kNN_inf;
   
   NumericVector y (nx, 0.0); // normalized version of x
   NumericVector range(2); // range vector for running min-max
   NumericMatrix cm(nh, nx); //cost matrix
   NumericVector ret(nx-nh+1, 0.0); // return vetor
   
   
   // set initial maximum indices
   j = 0;
   jsup = j + nh;
   jmax = jsup -1;
   
   if(do_norm == 1){
      // get rmin and rmax for the start
      range = cpp_range(x, j, jsup);//cpp_range(x, 0, nh); // range = c(min, max)
      
      // normalize the first nh values of x and fill y
      cpp_norm01(x, y, j, nh, range[0], range[1]);
      counter_norm_full += 1;
      
   } else{
      y = x; 
      norm_exit_status = 4;
   }
   
   
   // get initial dtw value
   loc_threshold = threshold;
   double lb = 0;
   
   NumericMatrix tube(nh,2);
   if(use_lb == 1){
      cpp_set_tube(tube, h, ws);
      lb = get_lb(tube, y, j, jsup);
      
   }else{
      lb = -1;
   }
   // double get_lb(const NumericMatrix &tube, const NumericVector &x, 
   //               int j0, int jsup){
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW BEGIN      //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   int iBegin = 0;
   int iEnd = 0;
   int nanCounter_status = 1;
   
   double * p1 = new double [nh];
   double * p2 = new double [nh];
   double * ptmp;
   double mynan;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   
   if(loc_threshold <= lb){
      //  loc_threshold < lb < DTW 
      ret[j] = mynan;
      counter_lb += 1;
      
   }else{
      // fill the initial nh columns of the costmatrix
      cpp_cm(cm, y, h, j, jsup, nh);
      counter_cm_full += 1;
      last_cm_calc_j = j;
      
      
      //initialize a and b with NAN
      for(int i=0; i < nh; i++){
         p1[i] = mynan;
         p2[i] = mynan;
      }
      
      // first column
      p1[0] = cm(0, j);
      if(p1[0] > loc_threshold){
         ret[j] = mynan;
         counter_ea += 1;
         
      }else{
         
         iEnd   = std::min(nh, ws+1);
         for(int i=1; i < iEnd; i++){
            p1[i] = cm(i, j) + p1[i-1];
            if(p1[i] > loc_threshold) p1[i] = mynan;
         }
         
         for(int k = j +1 ; k < jsup; k++){
            nanCounter = 0;
            iBegin = k-ws-j;
            if(iBegin <= 0){
               *p2 = cm(0, k) + *(p1);
               if(*(p2) > loc_threshold){
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
            
            iEnd   = k+ws+1-j;
            if(iEnd >= nh){
               iEnd = nh;
            }else{
               *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
            }
            
            for (int i = iBegin; i < iEnd; i++){
               *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
               if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                  *(p2+i) = mynan;
                  nanCounter ++;
               }
            }
            
            if(nanCounter == nh) {
               nanCounter_status = 0;
               break;
            }
            ptmp=p1;
            p1 = p2;
            p2 = ptmp;
         }
         
         if(nanCounter_status == 1){
            ret[j] = *(p1+nh-1);//p1[nx-1]
         } else{
            ret[j] = mynan;
            counter_ea += 1;
         }
      }
      
   }// lower bounding if clause
   
   // Rcout << loc_threshold <<" ... "<< lb << " ... " << ret[j]<< "\n";
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW END        //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   // ret[0] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
   // ret[0] = cpp_dtw2vec_cm(cm, i, jsup, nh);
   loc_threshold = mymin(ret[0], threshold);
   
   // ret[0] can be NaN => in this case set bsfiw to Inf, for proper comparisons,
   // since comparisons with NaN are always FALSE
   double initial_bsfiw = mymin(ret[0], R_PosInf);
   
   // fill initially kNN_val
   if(kNNk > 0){
      kNN_inf = fill_kNN_Info(kNN_inf_list);
      initialize_kNN(kNN_inf, kNN_val, kNN_ix, kNN_lot_ix, lot_ix, kNNk,
                  initial_bsfiw, overlap_size);
   }
   
   
   bsfiw = initial_bsfiw; // best so far in window
   bsfiw_ix = 0; // last index of best so far in window
   
   
   // --------------------------------------start rolling----------------------------
   for(j=1; j < (nx-nh+1); j++){
      if((j - bsfiw_ix) >= overlap_size){
         bsfiw = threshold;
         bsfiw_ix = j;
      }
      
      jsup = j + nh; // the frontrunner index +1
      jmax = jsup -1; // the index of the frontrunner
      
      // update the discretization of the new window
      // onlx_frontrunner_updated = 0;
      if(do_norm == 1){
         if(pastcheck == 1){
            if(x[jmax] < range[0]){
               // set new rmin and update whole xd
               range[0] = x[jmax];
               cpp_norm01(x, y, j, jsup, range[0], range[1]);
               counter_norm_newEx += 1;
               norm_exit_status = 1;
               
            }else if(x[jmax] > range[1]){
               // set new rmax and update whole xd
               range[1] = x[jmax];
               cpp_norm01(x, y, j, jsup, range[0], range[1]);
               counter_norm_newEx += 1;
               norm_exit_status = 1;
               
            }else{
               // only normalize y[jmax]
               y[jmax] = (x[jmax] - range[0])/(range[1] - range[0]) ;
               counter_norm_1step += 1;
               norm_exit_status = 4;
               
            }
         }else{
            //complete reset
            range = cpp_range(x, j, jsup); // range = c(min, max)
            cpp_norm01(x, y, j, jsup, range[0], range[1]);
            counter_norm_full += 1;
            norm_exit_status = 3;
            
         }
         
         // the next time i increases, x[i] is left behind, if x[i] was the
         // min or max, then new range needs to be defined
         if((x[j] == range[0]) || (x[j] == range[1])){
            pastcheck = 0;
         }else{
            pastcheck = 1;
         }
         
         
      }//end if do_norm
      
      
      //  Lower bounding
      if(use_lb == 1){
         lb = get_lb(tube, y, j, jsup);
      }else{
         lb = -1;
      }
      
      //  bsfiw serves as threshold
      if(use_ea == 1){
         loc_threshold = mymin(bsfiw, threshold);
      }else{
         loc_threshold = threshold;
      }
      
      
      if(kNNk > 0){
         // update local threshold, set it to the k-nearest neighbor
         // loc_threshold = mymin(loc_threshold, kNN_val[kNNk-1]);
         loc_threshold = mymin(loc_threshold, kNN_inf.vmax);
      }
      
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW BEGIN      //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      if(loc_threshold <= lb){
         //  loc_threshold < lb < DTW 
         ret[j] = mynan;
         counter_lb += 1;
         
      }else{
         
         if(norm_exit_status < 4 || last_cm_calc_j < j-1 ){
            // full reset of cm
            cpp_cm(cm, y, h, j, jsup, nh);
            counter_cm_full += 1;
         } else {
            // front running column of cm
            cpp_cm(cm, y, h, jmax, jsup, nh);
            counter_cm_1step += 1;
         }
         last_cm_calc_j = j;
         
         
         nanCounter_status = 1;
         
         //initialize a and b with NAN
         for(int i=0; i < nh; i++){
            p1[i] = mynan;
            p2[i] = mynan;
         }
         
         // first column
         p1[0] = cm(0, j);
         if(p1[0] > loc_threshold){
            ret[j] = mynan;
            counter_ea += 1;
            
         }else{
            
            iEnd   = std::min(nh, ws+1);
            for(int i=1; i < iEnd; i++){
               p1[i] = cm(i, j) + p1[i-1];
               if(p1[i] > loc_threshold) p1[i] = mynan;
            }
            
            for(int k = j +1 ; k < jsup; k++){
               nanCounter = 0;
               iBegin = k-ws-j;
               if(iBegin <= 0){
                  *p2 = cm(0, k) + *(p1);
                  if(*(p2) > loc_threshold){
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
               
               iEnd   = k+ws+1-j;
               if(iEnd >= nh){
                  iEnd = nh;
               }else{
                  *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
               }
               
               for (int i = iBegin; i < iEnd; i++){
                  *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
                  if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                     *(p2+i) = mynan;
                     nanCounter ++;
                  }
               }
               
               if(nanCounter == nh) {
                  nanCounter_status = 0;
                  break;
               }
               ptmp=p1;
               p1 = p2;
               p2 = ptmp;
            }
            
            if(nanCounter_status == 1){
               ret[j] = *(p1+nh-1);//p1[nx-1]
               // Rcout << ret[j] << " ";
            } else{
               ret[j] = mynan;
               counter_ea += 1;
            }
         }
         
         
      }// lower bound if clause
      
      // Rcout << loc_threshold <<" ... "<< lb << " ... " << ret[j]<< "\n";
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW END        //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      
      
      // ret[j] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
      // ret[i] = cpp_dtw2vec_cm(cm, j, jsup, nh);
      
      
      
      
      if(kNNk > 0){
         // update set of k nearest neighbors if new value is smaller than 
         // k-best value so far
         
         if(debug == 1){
            debug_print_kNN(debug, "before", j, kNN_inf, kNN_val, kNN_ix);
            Rcout << "lot_ix: " << lot_ix << "\n";
            Rcout << kNN_lot_ix << "\n";
         }
         
         if(ret[j] == ret[j]){ // if not nan
            // Rcout << "\nj: "<< j << "is not nan: " << ret[j]<<"\n";
            if(j < kNN_inf.imax + overlap_size){
               // there is an overlap
               // Rcout << "step1";
               if(ret[j] < kNN_val[kNN_inf.which_imax]){
                  // because of overlap the new value needs to be smaller
                  // than the most recent best-sofar-value => so we replace
                  // this value with the new one
                  kick_imax_kNN_lot(kNN_val, kNN_ix, kNN_lot_ix, kNN_inf, ret[j], j, lot_ix);
                  
                  if(debug == 1){
                     Rcout << "---kick imax---\n";
                  }
                  // Rcout << "step2";
               }
            }else{
               // Rcout << "step3";
               // no overlap
               if(ret[j] < kNN_inf.vmax){
                  // new value needs to be at least smaller 
                  // than the farthest nearest neighbor
                  // update_kNN(kNN_val, kNN_ix, ret[j], j);   
                  kick_vmax_kNN_lot(kNN_val, kNN_ix, kNN_lot_ix, kNN_inf, ret[j], j, lot_ix);
                  // Rcout << "step4";
                  if(debug == 1){
                     Rcout << "---kick vmax---\n";
                  }
               }
            }
         }
         
         if(debug == 1){
            debug_print_kNN(debug, "after", j, kNN_inf, kNN_val, kNN_ix);
            Rcout << "lot_ix: " << lot_ix << "\n";
            Rcout << kNN_lot_ix << "\n";
         }
      }
      
      
      if(ret[j] < bsfiw){
         bsfiw = ret[j];
         bsfiw_ix = j;
      }
      
   }//end rolling
   
   
   List list_ret;
   list_ret["dist"] = ret;
   IntegerVector counter(7);
   counter[0] = counter_norm_full;
   counter[1] = counter_norm_newEx;
   counter[2] = counter_norm_1step;
   counter[3] = counter_cm_full;
   counter[4] = counter_cm_1step;
   counter[5] = counter_ea;
   counter[6] = counter_lb;
   list_ret["counter"] = counter;
   
   if(kNNk > 0){
      // final reverse minimum search of ret:
      // up to here the kNN values and indices dealt as additional 
      // local threshold, but since the possibility of protracting
      // a local minimum, a final reverse check is necessary:
      
      list_ret["all_best_indices"] = cpp_kNN_rev(ret, overlap_size, debug);
      // list_ret["kNN_val"] = kNN_val;
      // list_ret["kNN_ix"] = kNN_ix;
      // list_ret["kNN_lot_ix"] = kNN_lot_ix;
      
      kNN_inf_list["imax"]        = kNN_inf.imax;
      kNN_inf_list["nr_detected"] = kNN_inf.nr_detected;
      kNN_inf_list["nr_looking4"] = kNN_inf.nr_looking4;
      kNN_inf_list["vmax"]        = kNN_inf.vmax;
      kNN_inf_list["which_imax"]  = kNN_inf.which_imax;
      kNN_inf_list["which_vmax"]  = kNN_inf.which_vmax;
      
      list_ret["kNN_inf"] = kNN_inf_list;
   }
   delete[] p1;
   delete[] p2;
   
   return list_ret;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_rundtw_mv(const NumericMatrix &h,
                   const NumericMatrix &x,
                   std::string step_pattern, std::string dist_method,
                   List &kNN_inf_list,
                   int ws, double threshold, int overlap_tol = 0, 
                   int do_norm = 1, int use_ea = 1, int use_lb = 1, int debug = 0) {
   
   // scan x with h and get the diff of DL(ss(x)) - DL(ss(x)-h) for each
   // subsegment ss(x)
   // h ... normalized
   // x ... not normalized
   // ws ... window size parameter for sakoe chiba window
   // threshold ... global threshold for early abandoning dtw calc
   //             (if threshold = Inf, then no thresholding takes place,
   //              in that case use_ea should be 0)
   // do_norm, use_ea, use_lb ... either 0 or 1
   
   // set dist_fun
   SEXP cm_SEXP = select_cm(dist_method);
   XPtr<funcPtr_cm> xpfun_cm(cm_SEXP);
   funcPtr_cm cm_fun = *xpfun_cm;
   
   // set step function
   SEXP step_SEXP = selectVecStepRun(step_pattern);
   XPtr<funcPtr_step_run> xpfun_step(step_SEXP);
   funcPtr_step_run step_fun = *xpfun_step;
   
   // set lb function
   SEXP lb_SEXP = select_lb(dist_method);
   XPtr<funcPtr_lb> xpfun_lb(lb_SEXP);
   funcPtr_lb lb_fun = *xpfun_lb;
   
   if(use_lb == 1){
      use_ea = 1; // lb makes no sense if the loc_threshold is not set by bsfiw
   }
   
   if(use_ea == 0){
      threshold = R_PosInf;
   }
   
   int nh = h.nrow();
   int overlap_size = nh - overlap_tol;
   int nx = x.nrow();
   int nc = x.ncol();
   int k, j,c, jmax, jsup, bsfiw_ix;
   int counter_norm_full = 0; // ... 1 for the initial step
   int counter_norm_1step = 0;
   int counter_norm_newEx = 0;
   int counter_cm_full = 0;
   int counter_cm_1step = 0;
   int counter_lb = 0;
   int counter_ea = 0;
   
   int norm_exit_status = 0;
   int last_cm_calc_j = -1;
   // int j_cm0 = 0;
   double bsfiw, loc_threshold; // local threshold
   double lb = 0;
   
   int kNNk = as<int>(kNN_inf_list["nr_looking4"]);
   IntegerVector kNN_ix (kNNk);
   NumericVector kNN_val (kNNk, 0.0);
   kNN_Info kNN_inf;
   
   NumericMatrix y (nx, nc); // normalized version of x
   NumericMatrix range(2, nc); // range vector for running min-max
   NumericMatrix cm(nh, nx); //cost matrix
   NumericVector ret(nx-nh+1, 0.0); // return vetor
   IntegerVector pastcheck(nc, 0);
   
   // set initial maximum indices
   j = 0;
   jsup = j + nh;
   jmax = jsup -1;
   
   
   if(do_norm == 1){				
      // get rmin and rmax for the start
      for(c = 0; c < nc; c++){
         cpp_set_range_mv(range, x, c, j, jsup);//cpp_range(x, 0, nh); // range = c(min, max)
      }
      
      // normalize the first nh values of x and fill y
      for(c = 0; c < nc; c++){
         cpp_norm01_mv(x, y, c, j, jsup, range(0, c), range(1,c));
         counter_norm_full += 1;
      }
      
   } else{
      y= x; 
      // j_cm0 = jmax;
      norm_exit_status = 4;
   }   
   
   // get initial dtw value
   loc_threshold = threshold;
   
   NumericMatrix tube(nh, 2*nc);
   if(use_lb == 1){
      cpp_set_tube_mv(tube, h, ws);
      lb = lb_fun(tube, y, j, jsup, nc);
      
   }else{
      lb = -1;
   }
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW BEGIN      //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   int iBegin = 0;
   int iEnd = 0;
   int nanCounter_status = 1;
   
   double * p1 = new double [nh];
   double * p2 = new double [nh];
   double * ptmp;
   double mynan;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   
   
   if(loc_threshold <= lb){
      //  loc_threshold < lb < DTW 
      ret[j] = mynan;
      counter_lb += 1;
      
   }else{							   
      
      // fill the initial nh columns of the costmatrix
      // cpp_cm1_mv(cm, y, h, j, jsup, nh, nc);
      cm_fun(cm, y, h, j, jsup, nh, nc);
      
      //initialize a and b with NAN
      for(int i=0; i < nh; i++){
         p1[i] = mynan;
         p2[i] = mynan;
      }
      
      // first column
      p1[0] = cm(0, j);
      if(p1[0] > loc_threshold){
         ret[j] = mynan;
         
         counter_ea += 1;
      }else{
         
         iEnd   = std::min(nh, ws+1);
         for(int i=1; i < iEnd; i++){
            p1[i] = cm(i, j) + p1[i-1];
            if(p1[i] > loc_threshold) p1[i] = mynan;
         }
         
         for(int k = j +1 ; k < jsup; k++){
            nanCounter = 0;
            iBegin = k-ws-j;
            if(iBegin <= 0){
               *p2 = cm(0, k) + *(p1);
               if(*(p2) > loc_threshold){
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
            
            iEnd   = k+ws+1-j;
            if(iEnd >= nh){
               iEnd = nh;
            }else{
               *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
            }
            
            for (int i = iBegin; i < iEnd; i++){
               *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
               if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                  *(p2+i) = mynan;
                  nanCounter ++;
               }
            }
            
            if(nanCounter == nh) {
               nanCounter_status = 0;
               break;
            }
            ptmp=p1;
            p1 = p2;
            p2 = ptmp;
         }
         
         if(nanCounter_status == 1){
            ret[j] = *(p1+nh-1);//p1[nx-1]
         } else{
            ret[j] = mynan;
            counter_ea += 1;
            
         }
      }
   }// lower bounding if clause  
   
   
   
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW END        //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   // ret[0] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
   // ret[0] = cpp_dtw2vec_cm(cm, i, jsup, nh);
   loc_threshold = mymin(ret[0], threshold);
   
   //the first bsfiw can be NaN => set it to Inf  if necessary, for proper comparisons 
   double initial_bsfiw = mymin(ret[0], R_PosInf);
   
   bsfiw = initial_bsfiw; // best so far in window
   bsfiw_ix = 0; // last index of of best so far in window
   
   
   // fill initially kNN_val
   if(kNNk > 0){
      kNN_val[0] = initial_bsfiw;
      kNN_ix[0] = 0;
      
      for(k = 1; k < kNNk; k++){
         kNN_val[k] = R_PosInf;
         kNN_ix[k] = -99;
      }
      
      kNN_inf = fill_kNN_Info(kNN_inf_list);
   }
   
   
   // --------------------------------------start rolling----------------------------
   for(j=1; j < (nx-nh+1); j++){
      
      if((j - bsfiw_ix) >= overlap_size){
         bsfiw = threshold;
         bsfiw_ix = j;
      }
      
      jsup = j + nh; // the frontrunner index +1
      jmax = jsup -1; // the index of the frontrunner
      
      // update the normalization of the new window
      // j_cm0 = jmax; // if one dimension disagrees and needs new normalization
      norm_exit_status = 4;
      // from scratch => j_cm0 needs to be set to j
      if(do_norm == 1){
         for(c = 0; c < nc; c++){
            if(pastcheck[c] == 1){
               if(x(jmax, c) < range(0,c)){
                  // set new min and update whole y
                  range(0,c) = x(jmax,c);
                  cpp_norm01_mv(x, y, c, j, jsup, range(0, c), range(1,c));
                  counter_norm_newEx += 1;
                  norm_exit_status = 2;
                  // j_cm0 = j;
                  
               }else if(x(jmax, c) > range(1,c)){
                  // set new max and update whole y
                  range(1,c) = x(jmax,c);
                  cpp_norm01_mv(x, y, c, j, jsup, range(0, c), range(1,c));
                  counter_norm_newEx += 1;
                  norm_exit_status = 2;
                  // j_cm0 = j;
                  
               }else{
                  // only normalize y[jmax]
                  y(jmax,c) = (x(jmax,c) - range(0,c))/(range(1,c) - range(0,c)) ;
                  counter_norm_1step += 1;
                  
               }
            }else{
               //complete reset
               cpp_set_range_mv(range, x, c, j, jsup);
               cpp_norm01_mv(x, y, c, j, jsup, range(0, c), range(1,c));
               counter_norm_full += 1;
               norm_exit_status = 2;
               // j_cm0 = j;
               
            }
         }
         
         // the next time i increases, x[i] is left behind, if x[i] was the
         // min or max, then new range needs to be defined
         for(c = 0; c < nc; c++){
            if((x(j,c) == range(0,c)) || (x(j,c) == range(1,c))){
               pastcheck[c] = 0;
            }else{
               pastcheck[c] = 1;
            }
         }
         
      }//end 	if do_norm
      
      
      //  Lower bounding
      if(use_lb == 1){
         lb = lb_fun(tube, y, j, jsup, nc);
         
      }else{
         lb = -1;
         
      }
      
      //  bsfiw serves as threshold
      if(use_ea == 1){
         loc_threshold = mymin(bsfiw, threshold);
      }else{
         loc_threshold = threshold;
      }
      
      if(kNNk > 0){
         // update local threshold, set it to the k-nearest neighbor
         // loc_threshold = mymin(loc_threshold, kNN_val[kNNk-1]);
         loc_threshold = mymin(loc_threshold, kNN_inf.vmax);
      }
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW BEGIN      //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      
      if(loc_threshold <= lb){
         //  loc_threshold < lb < DTW 
         ret[j] = mynan;
         counter_lb += 1;
         
      }else{
         
         // Rcout << norm_exit_status <<" ";
         // norm_exit_status = 2;
         
         if(norm_exit_status == 2 || last_cm_calc_j < j-1 ){
            // full reset of cm
            cm_fun(cm, y, h, j, jsup, nh, nc);
            counter_cm_full += 1;
         } else {
            // front running column of cm
            cm_fun(cm, y, h, jmax, jsup, nh, nc);
            counter_cm_1step += 1;
         }
         last_cm_calc_j = j;
         
         nanCounter_status = 1;
         
         //initialize a and b with NAN
         for(int i=0; i < nh; i++){
            p1[i] = mynan;
            p2[i] = mynan;
         }
         
         // first column
         p1[0] = cm(0, j);
         if(p1[0] > loc_threshold){
            ret[j] = mynan;
            
            counter_ea += 1;
         }else{
            
            iEnd   = std::min(nh, ws+1);
            for(int i=1; i < iEnd; i++){
               p1[i] = cm(i, j) + p1[i-1];
               if(p1[i] > loc_threshold) p1[i] = mynan;
            }
            
            for(int k = j +1 ; k < jsup; k++){
               nanCounter = 0;
               iBegin = k-ws-j;
               if(iBegin <= 0){
                  *p2 = cm(0, k) + *(p1);
                  if(*(p2) > loc_threshold){
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
               
               iEnd   = k+ws+1-j;
               if(iEnd >= nh){
                  iEnd = nh;
               }else{
                  *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
               }
               
               for (int i = iBegin; i < iEnd; i++){
                  *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
                  if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                     *(p2+i) = mynan;
                     nanCounter ++;
                  }
               }
               
               if(nanCounter == nh) {
                  nanCounter_status = 0;
                  break;
               }
               ptmp=p1;
               p1 = p2;
               p2 = ptmp;
            }
            
            if(nanCounter_status == 1){
               ret[j] = *(p1+nh-1);//p1[nx-1]
               
            } else{
               ret[j] = mynan;
               counter_ea += 1;
               
            }
         }
      }// lower bound if clause
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW END        //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      
      
      if(kNNk > 0){
         // update set of k nearest neighbors if new value is smaller than 
         // k-best value so far
         
         debug_print_kNN(debug, "before", j, kNN_inf, kNN_val, kNN_ix);
         
         if(ret[j] == ret[j]){ // if not nan
            if(j < kNN_inf.imax + overlap_size){
               // there is an overlap
               if(ret[j] < kNN_val[kNN_inf.which_imax]){
                  // because of overlap the new value needs to be smaller
                  // than the most recent best-sofar-value => so we replace
                  // this value with the new one
                  kick_imax_kNN(kNN_val, kNN_ix, kNN_inf, ret[j], j);
                  
                  if(debug == 1){
                     Rcout << "---kick imax---\n";
                  }
               }
            }else{
               // no overlap
               if(ret[j] < kNN_inf.vmax){
                  // new value needs to be at least smaller 
                  // than the farthest nearest neighbor
                  // update_kNN(kNN_val, kNN_ix, ret[j], j);   
                  kick_vmax_kNN(kNN_val, kNN_ix, kNN_inf, ret[j], j);
                  
                  if(debug == 1){
                     Rcout << "---kick vmax---\n";
                  }
               }
            }
         }
         debug_print_kNN(debug, "after", j, kNN_inf, kNN_val, kNN_ix);
         
      }
      
      // ret[j] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
      // ret[i] = cpp_dtw2vec_cm(cm, j, jsup, nh);
      
      if(ret[j] < bsfiw){
         bsfiw = ret[j];
         bsfiw_ix = j;
      }
   }
   
   List list_ret;
   list_ret["dist"] = ret;
   IntegerVector counter(7);
   counter[0] = counter_norm_full;
   counter[1] = counter_norm_newEx;
   counter[2] = counter_norm_1step;
   counter[3] = counter_cm_full;
   counter[4] = counter_cm_1step;
   counter[5] = counter_ea;
   counter[6] = counter_lb;
   
   list_ret["counter"] = counter;
   
   if(kNNk > 0){
      //TODO: final reverse minimum search of ret?!
      // list_ret["knn_values"] = kNN_val;
      // list_ret["knn_indices"] = kNN_ix;
      list_ret["all_best_indices"] = cpp_kNN_rev(ret, overlap_size, debug);
   }
   
   delete[] p1;
   delete[] p2;
   
   return list_ret;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_rundtw_mv_lot(const NumericMatrix &h,
                         const NumericMatrix &x,
                         NumericVector kNN_val_in,
                         IntegerVector kNN_ix_in,
                         IntegerVector kNN_lot_ix_in,
                         List &kNN_inf_list_in, 
                         int lot_ix,
                         std::string step_pattern, std::string dist_method,
                         int ws, double threshold, int overlap_tol = 0, 
                         int do_norm = 1, int use_ea = 1, int use_lb = 1, int debug = 0) {
         
         // scan x with h and get the diff of DL(ss(x)) - DL(ss(x)-h) for each
         // subsegment ss(x)
         // h ... normalized
         // x ... not normalized
         // ws ... window size parameter for sakoe chiba window
         // threshold ... global threshold for early abandoning dtw calc
         //             (if threshold = Inf, then no thresholding takes place,
         //              in that case use_ea should be 0)
         // do_norm, use_ea, use_lb ... either 0 or 1
         
         NumericVector kNN_val    = clone(kNN_val_in);
         IntegerVector kNN_ix     = clone(kNN_ix_in);
         IntegerVector kNN_lot_ix = clone(kNN_lot_ix_in);
         List kNN_inf_list        = clone(kNN_inf_list_in);
   
   
         // set dist_fun
         SEXP cm_SEXP = select_cm(dist_method);
         XPtr<funcPtr_cm> xpfun_cm(cm_SEXP);
         funcPtr_cm cm_fun = *xpfun_cm;
         
         // set step function
         SEXP step_SEXP = selectVecStepRun(step_pattern);
         XPtr<funcPtr_step_run> xpfun_step(step_SEXP);
         funcPtr_step_run step_fun = *xpfun_step;
         
         // set lb function
         SEXP lb_SEXP = select_lb(dist_method);
         XPtr<funcPtr_lb> xpfun_lb(lb_SEXP);
         funcPtr_lb lb_fun = *xpfun_lb;
         
         if(use_lb == 1){
            use_ea = 1; // lb makes no sense if the loc_threshold is not set by bsfiw
         }
         
         if(use_ea == 0){
            threshold = R_PosInf;
         }
         
         int nh = h.nrow();
         int overlap_size = nh - overlap_tol;
         int nx = x.nrow();
         int nc = x.ncol();
         int j,c, jmax, jsup, bsfiw_ix;
         int counter_norm_full = 0; // ... 1 for the initial step
         int counter_norm_1step = 0;
         int counter_norm_newEx = 0;
         int counter_cm_full = 0;
         int counter_cm_1step = 0;
         int counter_lb = 0;
         int counter_ea = 0;
         
         int norm_exit_status = 0;
         int last_cm_calc_j = -1;
         // int j_cm0 = 0;
         double bsfiw, loc_threshold; // local threshold
         double lb = 0;
         
         int kNNk = as<int>(kNN_inf_list["nr_looking4"]);
         kNN_Info kNN_inf;
         
         NumericMatrix y (nx, nc); // normalized version of x
         NumericMatrix range(2, nc); // range vector for running min-max
         NumericMatrix cm(nh, nx); //cost matrix
         NumericVector ret(nx-nh+1, 0.0); // return vetor
         IntegerVector pastcheck(nc, 0);
         
         // set initial maximum indices
         j = 0;
         jsup = j + nh;
         jmax = jsup -1;
         
         
         if(do_norm == 1){				
            // get rmin and rmax for the start
            for(c = 0; c < nc; c++){
               cpp_set_range_mv(range, x, c, j, jsup);//cpp_range(x, 0, nh); // range = c(min, max)
            }
            
            // normalize the first nh values of x and fill y
            for(c = 0; c < nc; c++){
               cpp_norm01_mv(x, y, c, j, jsup, range(0, c), range(1,c));
               counter_norm_full += 1;
            }
            
         } else{
            y= x; 
            // j_cm0 = jmax;
            norm_exit_status = 4;
         }   
         
         // get initial dtw value
         loc_threshold = threshold;
         
         NumericMatrix tube(nh, 2*nc);
         if(use_lb == 1){
            cpp_set_tube_mv(tube, h, ws);
            lb = lb_fun(tube, y, j, jsup, nc);
            
         }else{
            lb = -1;
         }
         
         ////////////////////////////////////////////////////////////////////////////////
         //////////////////////////////     1st DTW BEGIN      //////////////////////////
         ////////////////////////////////////////////////////////////////////////////////
         
         
         int iBegin = 0;
         int iEnd = 0;
         int nanCounter_status = 1;
         
         double * p1 = new double [nh];
         double * p2 = new double [nh];
         double * ptmp;
         double mynan;
         int nanCounter;
         mynan = std::numeric_limits<double>::quiet_NaN();
         
         
         
         if(loc_threshold <= lb){
            //  loc_threshold < lb < DTW 
            ret[j] = mynan;
            counter_lb += 1;
            
         }else{							   
            
            // fill the initial nh columns of the costmatrix
            // cpp_cm1_mv(cm, y, h, j, jsup, nh, nc);
            cm_fun(cm, y, h, j, jsup, nh, nc);
            
            //initialize a and b with NAN
            for(int i=0; i < nh; i++){
               p1[i] = mynan;
               p2[i] = mynan;
            }
            
            // first column
            p1[0] = cm(0, j);
            if(p1[0] > loc_threshold){
               ret[j] = mynan;
               
               counter_ea += 1;
            }else{
               
               iEnd   = std::min(nh, ws+1);
               for(int i=1; i < iEnd; i++){
                  p1[i] = cm(i, j) + p1[i-1];
                  if(p1[i] > loc_threshold) p1[i] = mynan;
               }
               
               for(int k = j +1 ; k < jsup; k++){
                  nanCounter = 0;
                  iBegin = k-ws-j;
                  if(iBegin <= 0){
                     *p2 = cm(0, k) + *(p1);
                     if(*(p2) > loc_threshold){
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
                  
                  iEnd   = k+ws+1-j;
                  if(iEnd >= nh){
                     iEnd = nh;
                  }else{
                     *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
                  }
                  
                  for (int i = iBegin; i < iEnd; i++){
                     *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
                     if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                        *(p2+i) = mynan;
                        nanCounter ++;
                     }
                  }
                  
                  if(nanCounter == nh) {
                     nanCounter_status = 0;
                     break;
                  }
                  ptmp=p1;
                  p1 = p2;
                  p2 = ptmp;
               }
               
               if(nanCounter_status == 1){
                  ret[j] = *(p1+nh-1);//p1[nx-1]
               } else{
                  ret[j] = mynan;
                  counter_ea += 1;
                  
               }
            }
         }// lower bounding if clause  
         
         
         
         
         ////////////////////////////////////////////////////////////////////////////////
         //////////////////////////////     1st DTW END        //////////////////////////
         ////////////////////////////////////////////////////////////////////////////////
         
         
         // ret[0] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
         // ret[0] = cpp_dtw2vec_cm(cm, i, jsup, nh);
         loc_threshold = mymin(ret[0], threshold);
         
         //the first bsfiw can be NaN => set it to Inf  if necessary, for proper comparisons 
         double initial_bsfiw = mymin(ret[0], R_PosInf);
         
         bsfiw = initial_bsfiw; // best so far in window
         bsfiw_ix = 0; // last index of of best so far in window
         
         
         // fill initially kNN_val
         if(kNNk > 0){
            kNN_inf = fill_kNN_Info(kNN_inf_list);
            initialize_kNN(kNN_inf, kNN_val, kNN_ix, kNN_lot_ix, lot_ix, kNNk,
                        initial_bsfiw, overlap_size);
         }
         
         
         // --------------------------------------start rolling----------------------------
         for(j=1; j < (nx-nh+1); j++){
            
            if((j - bsfiw_ix) >= overlap_size){
               bsfiw = threshold;
               bsfiw_ix = j;
            }
            
            jsup = j + nh; // the frontrunner index +1
            jmax = jsup -1; // the index of the frontrunner
            
            // update the normalization of the new window
            // j_cm0 = jmax; // if one dimension disagrees and needs new normalization
            norm_exit_status = 4;
            // from scratch => j_cm0 needs to be set to j
            if(do_norm == 1){
               for(c = 0; c < nc; c++){
                  if(pastcheck[c] == 1){
                     if(x(jmax, c) < range(0,c)){
                        // set new min and update whole y
                        range(0,c) = x(jmax,c);
                        cpp_norm01_mv(x, y, c, j, jsup, range(0, c), range(1,c));
                        counter_norm_newEx += 1;
                        norm_exit_status = 2;
                        // j_cm0 = j;
                        
                     }else if(x(jmax, c) > range(1,c)){
                        // set new max and update whole y
                        range(1,c) = x(jmax,c);
                        cpp_norm01_mv(x, y, c, j, jsup, range(0, c), range(1,c));
                        counter_norm_newEx += 1;
                        norm_exit_status = 2;
                        // j_cm0 = j;
                        
                     }else{
                        // only normalize y[jmax]
                        y(jmax,c) = (x(jmax,c) - range(0,c))/(range(1,c) - range(0,c)) ;
                        counter_norm_1step += 1;
                        
                     }
                  }else{
                     //complete reset
                     cpp_set_range_mv(range, x, c, j, jsup);
                     cpp_norm01_mv(x, y, c, j, jsup, range(0, c), range(1,c));
                     counter_norm_full += 1;
                     norm_exit_status = 2;
                     // j_cm0 = j;
                     
                  }
               }
               
               // the next time i increases, x[i] is left behind, if x[i] was the
               // min or max, then new range needs to be defined
               for(c = 0; c < nc; c++){
                  if((x(j,c) == range(0,c)) || (x(j,c) == range(1,c))){
                     pastcheck[c] = 0;
                  }else{
                     pastcheck[c] = 1;
                  }
               }
               
            }//end 	if do_norm
            
            
            //  Lower bounding
            if(use_lb == 1){
               lb = lb_fun(tube, y, j, jsup, nc);
               
            }else{
               lb = -1;
               
            }
            
            //  bsfiw serves as threshold
            if(use_ea == 1){
               loc_threshold = mymin(bsfiw, threshold);
            }else{
               loc_threshold = threshold;
            }
            
            if(kNNk > 0){
               // update local threshold, set it to the k-nearest neighbor
               // loc_threshold = mymin(loc_threshold, kNN_val[kNNk-1]);
               loc_threshold = mymin(loc_threshold, kNN_inf.vmax);
            }
            
            ////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////     RUN DTW BEGIN      //////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            
            if(loc_threshold <= lb){
               //  loc_threshold < lb < DTW 
               ret[j] = mynan;
               counter_lb += 1;
               
            }else{
               
               // Rcout << norm_exit_status <<" ";
               // norm_exit_status = 2;
               
               if(norm_exit_status == 2 || last_cm_calc_j < j-1 ){
                  // full reset of cm
                  cm_fun(cm, y, h, j, jsup, nh, nc);
                  counter_cm_full += 1;
               } else {
                  // front running column of cm
                  cm_fun(cm, y, h, jmax, jsup, nh, nc);
                  counter_cm_1step += 1;
               }
               last_cm_calc_j = j;
               
               nanCounter_status = 1;
               
               //initialize a and b with NAN
               for(int i=0; i < nh; i++){
                  p1[i] = mynan;
                  p2[i] = mynan;
               }
               
               // first column
               p1[0] = cm(0, j);
               if(p1[0] > loc_threshold){
                  ret[j] = mynan;
                  
                  counter_ea += 1;
               }else{
                  
                  iEnd   = std::min(nh, ws+1);
                  for(int i=1; i < iEnd; i++){
                     p1[i] = cm(i, j) + p1[i-1];
                     if(p1[i] > loc_threshold) p1[i] = mynan;
                  }
                  
                  for(int k = j +1 ; k < jsup; k++){
                     nanCounter = 0;
                     iBegin = k-ws-j;
                     if(iBegin <= 0){
                        *p2 = cm(0, k) + *(p1);
                        if(*(p2) > loc_threshold){
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
                     
                     iEnd   = k+ws+1-j;
                     if(iEnd >= nh){
                        iEnd = nh;
                     }else{
                        *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
                     }
                     
                     for (int i = iBegin; i < iEnd; i++){
                        *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
                        if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                           *(p2+i) = mynan;
                           nanCounter ++;
                        }
                     }
                     
                     if(nanCounter == nh) {
                        nanCounter_status = 0;
                        break;
                     }
                     ptmp=p1;
                     p1 = p2;
                     p2 = ptmp;
                  }
                  
                  if(nanCounter_status == 1){
                     ret[j] = *(p1+nh-1);//p1[nx-1]
                     
                  } else{
                     ret[j] = mynan;
                     counter_ea += 1;
                     
                  }
               }
            }// lower bound if clause
            
            ////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////     RUN DTW END        //////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            
            
            if(kNNk > 0){
               // update set of k nearest neighbors if new value is smaller than 
               // k-best value so far
               
               if(ret[j] == ret[j]){ // if not nan
                  if(j < kNN_inf.imax + overlap_size){
                     // there is an overlap
                     if(ret[j] < kNN_val[kNN_inf.which_imax]){
                        // because of overlap the new value needs to be smaller
                        // than the most recent best-sofar-value => so we replace
                        // this value with the new one
                        kick_imax_kNN_lot(kNN_val, kNN_ix, kNN_lot_ix, kNN_inf, ret[j], j, lot_ix);
                        
                     }
                  }else{
                     // no overlap
                     if(ret[j] < kNN_inf.vmax){
                        // new value needs to be at least smaller 
                        // than the farthest nearest neighbor
                        // update_kNN(kNN_val, kNN_ix, ret[j], j);   
                        kick_vmax_kNN_lot(kNN_val, kNN_ix, kNN_lot_ix, kNN_inf, ret[j], j, lot_ix);
                     }
                  }
               }
            }
            
            if(ret[j] < bsfiw){
               bsfiw = ret[j];
               bsfiw_ix = j;
            }
         }// end rolling
         
         List list_ret;
         list_ret["dist"] = ret;
         IntegerVector counter(7);
         counter[0] = counter_norm_full;
         counter[1] = counter_norm_newEx;
         counter[2] = counter_norm_1step;
         counter[3] = counter_cm_full;
         counter[4] = counter_cm_1step;
         counter[5] = counter_ea;
         counter[6] = counter_lb;
         
         list_ret["counter"] = counter;
         
         if(kNNk > 0){
            // final reverse minimum search of ret:
            // up to here the kNN values and indices dealt as additional 
            // local threshold, but since the possibility of protracting
            // a local minimum, a final reverse check is necessary:
            
            list_ret["all_best_indices"] = cpp_kNN_rev(ret, overlap_size, debug);
            
            kNN_inf_list["imax"]        = kNN_inf.imax;
            kNN_inf_list["nr_detected"] = kNN_inf.nr_detected;
            kNN_inf_list["nr_looking4"] = kNN_inf.nr_looking4;
            kNN_inf_list["vmax"]        = kNN_inf.vmax;
            kNN_inf_list["which_imax"]  = kNN_inf.which_imax;
            kNN_inf_list["which_vmax"]  = kNN_inf.which_vmax;
            list_ret["kNN_inf"] = kNN_inf_list;
         }
         
         delete[] p1;
         delete[] p2;
         
         return list_ret;
      }


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// [[Rcpp::export]]
List cpp_rundtw_znorm(const NumericVector &h,
                      const NumericVector &x,
                      std::string step_pattern,
                      List &kNN_inf_list,
                      int ws, double threshold, int overlap_tol = 0, 
                      int use_ea = 1, int use_lb = 1, int debug = 0) {
   // scan x with h and get the diff of DL(ss(x)) - DL(ss(x)-h) for each
   // subsegment ss(x)
   // h ... normalized
   // x ... not normalized
   // ws ... window size parameter for sakoe chiba window
   // threshold ... global threshold for early abandoning dtw calc
   //             (if threshold = Inf, then no thresholding takes place,
   //              in that case use_ea should be 0)
   // use_ea, use_lb ... either 0 or 1
   
   
   // set step function
   SEXP step_SEXP = selectVecStepRun(step_pattern);
   XPtr<funcPtr_step_run> xpfun_step(step_SEXP);
   funcPtr_step_run step_fun = *xpfun_step;
   
   if(use_lb == 1){
      use_ea = 1; // lb makes no sense if the loc_threshold is not set by bsfiw
   }
   
   if(use_ea == 0){
      threshold = R_PosInf;
   }
   
   int nh = h.size();
   
   int overlap_size = nh - overlap_tol;
   int nx = x.size();
   int j, k, jmax, jsup, bsfiw_ix;
   
   int counter_lb = 0;
   int counter_ea = 0;
   
   
   double bsfiw, loc_threshold; // local threshold
   
   int kNNk = as<int>(kNN_inf_list["nr_looking4"]);
   IntegerVector kNN_ix (kNNk);
   NumericVector kNN_val (kNNk, 0.0);
   kNN_Info kNN_inf;
   
   double mu0 = 0;
   double mu = 0;
   double mu2 = 0;
   double sigma = 0;
   double var_precision = 0.000000001 ; // 1e-9
   double sigma2 = 0;
   double z_j = 0; 
   NumericVector ret(nx-nh+1, 0.0); // return vetor
   
   
   
   // set initial indices: it is possible that the first 
   //                      some entries of ret are NaN
   j = 0;
   while(ret[j] != ret[j] && j < (nx-nh+1)){
      j += 1;
   }
   jsup = j + nh;
   jmax = jsup -1;
   
   // set initial mu and sigma
   mu = get_mean(x, j, jmax);
   mu2 = mu * mu;
   sigma2 = get_sigma2(x, j, jmax, mu2);
   if(sigma2 < var_precision){
      sigma = 1; 
   } else{
      sigma = sqrt(sigma2);
   }
   
   if(debug == 1){
      Rcout << "j = "<< j << " ... mu: "<< mu << " "<< " sigma: "<< sigma << " sigma2: "<< sigma2 << "\n";
   }
   
   
   // get initial dtw value
   loc_threshold = threshold;
   double lb = 0;
   
   NumericMatrix tube(nh,2);
   if(use_lb == 1){
      cpp_set_tube(tube, h, ws);
      lb = get_lb_znorm(tube, x, mu, sigma, loc_threshold, j, jsup);
      
   }else{
      lb = -1;
   }
   // double get_lb(const NumericMatrix &tube, const NumericVector &x, 
   //               int j0, int jsup){
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW BEGIN      //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   int iBegin = 0;
   int iEnd = 0;
   int nanCounter_status = 1;
   
   double * p1 = new double [nh];
   double * p2 = new double [nh];
   double * ptmp;
   double mynan;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   
   if(loc_threshold <= lb){
      //  loc_threshold < lb < DTW 
      ret[j] = mynan;
      counter_lb += 1;
      
   }else{
      // fill the initial nh columns of the costmatrix
      // cpp_cm(cm, y, h, j, jsup, nh);
      // counter_cm_full += 1;
      // last_cm_calc_j = j;
      
      
      //initialize a and b with NAN
      for(int i=0; i < nh; i++){
         p1[i] = mynan;
         p2[i] = mynan;
      }
      
      // first column
      z_j = (x[j] - mu)/sigma;
      p1[0] = abs(z_j - h[0]);
      
      if(p1[0] > loc_threshold){
         ret[j] = mynan;
         counter_ea += 1;
         
      }else{
         
         iEnd = std::min(nh, ws+1);
         for(int i=1; i < iEnd; i++){
            
            p1[i] = abs(z_j - h[i]) + p1[i-1];
            // p1[i] = cm(i, j) + p1[i-1];
            if(p1[i] > loc_threshold) p1[i] = mynan;
         }
         
         for(int k = j +1 ; k < jsup; k++){
            
            z_j = (x[k] - mu)/sigma;
            
            nanCounter = 0;
            iBegin = k-ws-j;
            if(iBegin <= 0){
               
               
               *p2 = abs(z_j - h[0]) + *(p1);
               // *p2 = cm(0, k) + *(p1);
               if(*(p2) > loc_threshold){
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
            
            iEnd   = k+ws+1-j;
            if(iEnd >= nh){
               iEnd = nh;
            }else{
               *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
            }
            
            
            for (int i = iBegin; i < iEnd; i++){
               *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), abs(z_j - h[i]));
               // *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
               if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                  *(p2+i) = mynan;
                  nanCounter ++;
               }
            }
            
            if(nanCounter == nh) {
               nanCounter_status = 0;
               break;
            }
            ptmp=p1;
            p1 = p2;
            p2 = ptmp;
         }
         
         if(nanCounter_status == 1){
            ret[j] = *(p1+nh-1);//p1[nx-1]
         } else{
            ret[j] = mynan;
            counter_ea += 1;
         }
      }
      
   }// lower bounding if clause
   
   // Rcout << loc_threshold <<" ... "<< lb << " ... " << ret[j]<< "\n";
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW END        //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   // ret[0] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
   // ret[0] = cpp_dtw2vec_cm(cm, i, jsup, nh);
   loc_threshold = mymin(ret[0], threshold);
   
   // ret[0] can be NaN => in this case set bsfiw to Inf, for proper comparisons,
   // since comparisons with NaN are always FALSE
   double initial_bsfiw = mymin(ret[0], R_PosInf);
   
   
   
   // fill initially kNN_val
   if(kNNk > 0){
      kNN_val[0] = initial_bsfiw;
      kNN_ix[0] = 0;
      
      for(k = 1; k < kNNk; k++){
         kNN_val[k] = R_PosInf;
         kNN_ix[k] = -99;
      }
      
      kNN_inf = fill_kNN_Info(kNN_inf_list);
   }
   
   bsfiw = initial_bsfiw; // best so far in window
   bsfiw_ix = 0; // last index of best so far in window
   
   
   // --------------------------------------start rolling----------------------------
   int j_prev_norm = 0;
   j = 1;
   while(j < (nx-nh+1)){
      // for(j=1; j < (nx-nh+1); j++){
      
      // for lot-mode: test if the current value should be calculated 
      while(ret[j] != ret[j] && j < (nx-nh+1)){
         j += 1;
      }
      
      jsup = j + nh; // the frontrunner index +1
      jmax = jsup -1; // the index of the frontrunner
      
      
      if(j_prev_norm == j-1){
         // update the discretization of the new window
         mu0 = mu;
         mu += (x[jmax] - x[j-1])/nh;
         sigma2 += (x[jmax] * x[jmax] - x[j-1] * x[j-1])/(nh-1) + 
            (mu0 * mu0 - mu * mu) * nh/(nh-1);
         
      }else{
         //  discretization of the new window from scratch
         mu0 = mu;
         mu = get_mean(x, j, jmax);
         mu2 = mu * mu;
         sigma2 = get_sigma2(x, j, jmax, mu2);
      }
      
      if(sigma2 < var_precision){
         sigma = 1; 
      } else{
         sigma = sqrt(sigma2);
      }
      // Rcout << "\nsigma2: "<< sigma2 << "sigma: "<< sigma << "\n";
      
      j_prev_norm = j;
      
      if(debug == 1){
         Rcout << "j = "<< j << " ... mu: "<< mu << " "<< " sigma: "<< sigma << " sigma2: "<< sigma2 << "\n";
      }
      
      if((j - bsfiw_ix) >= overlap_size){
         bsfiw = threshold;
         bsfiw_ix = j;
      }
      
      
      
      //  bsfiw serves as threshold
      if(use_ea == 1){
         loc_threshold = mymin(bsfiw, threshold);
      }else{
         loc_threshold = threshold;
      }
      
      //  Lower bounding
      if(use_lb == 1){
         lb = get_lb_znorm(tube, x, mu, sigma, loc_threshold, j, jsup);
      }else{
         lb = -1;
      }
      
      if(kNNk > 0){
         // update local threshold, set it to the k-nearest neighbor
         // loc_threshold = mymin(loc_threshold, kNN_val[kNNk-1]);
         loc_threshold = mymin(loc_threshold, kNN_inf.vmax);
      }
      
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW BEGIN      //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      if(loc_threshold <= lb){
         //  loc_threshold < lb < DTW 
         ret[j] = mynan;
         counter_lb += 1;
         
      }else{
         
         nanCounter_status = 1;
         
         //initialize a and b with NAN
         for(int i=0; i < nh; i++){
            p1[i] = mynan;
            p2[i] = mynan;
         }
         
         // first column
         z_j = (x[j] - mu)/sigma;
         p1[0] = abs(z_j - h[0]);
         
         if(p1[0] > loc_threshold){
            ret[j] = mynan;
            counter_ea += 1;
            
         }else{
            
            iEnd   = std::min(nh, ws+1);
            for(int i=1; i < iEnd; i++){
               
               p1[i] = abs(z_j - h[i]) + p1[i-1];
               // p1[i] = cm(i, j) + p1[i-1];
               if(p1[i] > loc_threshold) p1[i] = mynan;
            }
            
            for(int k = j +1 ; k < jsup; k++){
               
               z_j = (x[k] - mu)/sigma;
               
               nanCounter = 0;
               iBegin = k-ws-j;
               if(iBegin <= 0){
                  // *p2 = cm(0, k) + *(p1);
                  *p2 = abs(z_j - h[0]) + *(p1);
                  if(*(p2) > loc_threshold){
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
               
               iEnd   = k+ws+1-j;
               if(iEnd >= nh){
                  iEnd = nh;
               }else{
                  *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
               }
               
               
               for (int i = iBegin; i < iEnd; i++){
                  *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), abs(z_j - h[i]));
                  if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                     *(p2+i) = mynan;
                     nanCounter ++;
                  }
               }
               
               if(nanCounter == nh) {
                  nanCounter_status = 0;
                  break;
               }
               ptmp=p1;
               p1 = p2;
               p2 = ptmp;
            }
            
            if(nanCounter_status == 1){
               ret[j] = *(p1+nh-1);//p1[nx-1]
               // Rcout << ret[j] << " ";
            } else{
               ret[j] = mynan;
               counter_ea += 1;
            }
         }
         
         
      }// lower bound if clause
      
      // Rcout << loc_threshold <<" ... "<< lb << " ... " << ret[j]<< "\n";
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW END        //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      
      
      // ret[j] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
      // ret[i] = cpp_dtw2vec_cm(cm, j, jsup, nh);
      
      
      
      
      if(kNNk > 0){
         // update set of k nearest neighbors if new value is smaller than 
         // k-best value so far
         
         debug_print_kNN(debug, "before", j, kNN_inf, kNN_val, kNN_ix);
         
         if(ret[j] == ret[j]){ // if not nan
            if(j < kNN_inf.imax + overlap_size){
               // there is an overlap
               if(ret[j] < kNN_val[kNN_inf.which_imax]){
                  // because of overlap the new value needs to be smaller
                  // than the most recent best-sofar-value => so we replace
                  // this value with the new one
                  kick_imax_kNN(kNN_val, kNN_ix, kNN_inf, ret[j], j);
                  
                  if(debug == 1){
                     Rcout << "---kick imax---\n";
                  }
               }
            }else{
               // no overlap
               if(ret[j] < kNN_inf.vmax){
                  // new value needs to be at least smaller 
                  // than the farthest nearest neighbor
                  // update_kNN(kNN_val, kNN_ix, ret[j], j);   
                  kick_vmax_kNN(kNN_val, kNN_ix, kNN_inf, ret[j], j);
                  
                  if(debug == 1){
                     Rcout << "---kick vmax---\n";
                  }
               }
            }
         }
         
         debug_print_kNN(debug, "after", j, kNN_inf, kNN_val, kNN_ix);
         
      }
      
      
      if(ret[j] < bsfiw){
         bsfiw = ret[j];
         bsfiw_ix = j;
      }
      
      j += 1;
   }//end rolling
   
   
   List list_ret;
   list_ret["dist"] = ret;
   IntegerVector counter(7);
   counter[0] = 0;
   counter[1] = 0;
   counter[2] = 0;
   counter[3] = 0;
   counter[4] = 0;
   counter[5] = counter_ea;
   counter[6] = counter_lb;
   list_ret["counter"] = counter;
   
   if(kNNk > 0){
      // final reverse minimum search of ret:
      // up to here the kNN values and indices dealt as additional 
      // local threshold, but since the possibility of protracting
      // a local minimum, a final reverse check is necessary:
      
      list_ret["all_best_indices"] = cpp_kNN_rev(ret, overlap_size, debug);
   }
   delete[] p1;
   delete[] p2;
   
   return list_ret;
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// [[Rcpp::export]]
List cpp_rundtw_znorm_lot(const NumericVector &h,
                            const NumericVector &x,
                            NumericVector kNN_val_in,
                            IntegerVector kNN_ix_in,
                            IntegerVector kNN_lot_ix_in,
                            List &kNN_inf_list_in, 
                            int lot_ix,
                            std::string step_pattern,
                            int ws, double threshold, int overlap_tol = 0, 
                            int use_ea = 1, int use_lb = 1, int debug = 0) {
   // scan x with h and get the diff of DL(ss(x)) - DL(ss(x)-h) for each
   // subsegment ss(x)
   // h ... normalized
   // x ... not normalized
   // ws ... window size parameter for sakoe chiba window
   // threshold ... global threshold for early abandoning dtw calc
   //             (if threshold = Inf, then no thresholding takes place,
   //              in that case use_ea should be 0)
   // use_ea, use_lb ... either 0 or 1
   
   NumericVector kNN_val    = clone(kNN_val_in);
   IntegerVector kNN_ix     = clone(kNN_ix_in);
   IntegerVector kNN_lot_ix = clone(kNN_lot_ix_in);
   List kNN_inf_list        = clone(kNN_inf_list_in);
   
   
   // set step function
   SEXP step_SEXP = selectVecStepRun(step_pattern);
   XPtr<funcPtr_step_run> xpfun_step(step_SEXP);
   funcPtr_step_run step_fun = *xpfun_step;
   
   if(use_lb == 1){
      use_ea = 1; // lb makes no sense if the loc_threshold is not set by bsfiw
   }
   
   if(use_ea == 0){
      threshold = R_PosInf;
   }
   
   int nh = h.size();
   
   int overlap_size = nh - overlap_tol;
   int nx = x.size();
   int j, jmax, jsup, bsfiw_ix;
   
   int counter_lb = 0;
   int counter_ea = 0;
   
   
   double bsfiw, loc_threshold; // local threshold
   
   int kNNk = as<int>(kNN_inf_list["nr_looking4"]);
   kNN_Info kNN_inf;
   
   double mu0 = 0;
   double mu = 0;
   double mu2 = 0;
   double sigma = 0;
   double var_precision = 0.000000001 ; // 1e-9
   double sigma2 = 0;
   double z_j = 0; 
   NumericVector ret(nx-nh+1, 0.0); // return vetor
   
   
   // set initial indices: it is possible that the first 
   //                      some entries of ret are NaN
   j = 0;
   while(ret[j] != ret[j] && j < (nx-nh+1)){
      j += 1;
   }
   jsup = j + nh;
   jmax = jsup -1;
   
   // set initial mu and sigma
   mu = get_mean(x, j, jmax);
   mu2 = mu * mu;
   sigma2 = get_sigma2(x, j, jmax, mu2);
   if(sigma2 < var_precision){
      sigma = 1; 
   } else{
      sigma = sqrt(sigma2);
   }
   
   if(debug == 1){
      Rcout << "j = "<< j << " ... mu: "<< mu << " "<< " sigma: "<< sigma << " sigma2: "<< sigma2 << "\n";
   }
   
   
   // get initial dtw value
   loc_threshold = threshold;
   double lb = 0;
   
   NumericMatrix tube(nh,2);
   if(use_lb == 1){
      cpp_set_tube(tube, h, ws);
      lb = get_lb_znorm(tube, x, mu, sigma, loc_threshold, j, jsup);
      
   }else{
      lb = -1;
   }
   // double get_lb(const NumericMatrix &tube, const NumericVector &x, 
   //               int j0, int jsup){
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW BEGIN      //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   int iBegin = 0;
   int iEnd = 0;
   int nanCounter_status = 1;
   
   double * p1 = new double [nh];
   double * p2 = new double [nh];
   double * ptmp;
   double mynan;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   
   if(loc_threshold <= lb){
      //  loc_threshold < lb < DTW 
      ret[j] = mynan;
      counter_lb += 1;
      
   }else{
      // fill the initial nh columns of the costmatrix
      // cpp_cm(cm, y, h, j, jsup, nh);
      // counter_cm_full += 1;
      // last_cm_calc_j = j;
      
      
      //initialize a and b with NAN
      for(int i=0; i < nh; i++){
         p1[i] = mynan;
         p2[i] = mynan;
      }
      
      // first column
      z_j = (x[j] - mu)/sigma;
      p1[0] = abs(z_j - h[0]);
      
      if(p1[0] > loc_threshold){
         ret[j] = mynan;
         counter_ea += 1;
         
      }else{
         
         iEnd = std::min(nh, ws+1);
         for(int i=1; i < iEnd; i++){
            
            p1[i] = abs(z_j - h[i]) + p1[i-1];
            // p1[i] = cm(i, j) + p1[i-1];
            if(p1[i] > loc_threshold) p1[i] = mynan;
         }
         
         for(int k = j +1 ; k < jsup; k++){
            
            z_j = (x[k] - mu)/sigma;
            
            nanCounter = 0;
            iBegin = k-ws-j;
            if(iBegin <= 0){
               
               
               *p2 = abs(z_j - h[0]) + *(p1);
               // *p2 = cm(0, k) + *(p1);
               if(*(p2) > loc_threshold){
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
            
            iEnd   = k+ws+1-j;
            if(iEnd >= nh){
               iEnd = nh;
            }else{
               *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
            }
            
            
            for (int i = iBegin; i < iEnd; i++){
               *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), abs(z_j - h[i]));
               // *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
               if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                  *(p2+i) = mynan;
                  nanCounter ++;
               }
            }
            
            if(nanCounter == nh) {
               nanCounter_status = 0;
               break;
            }
            ptmp=p1;
            p1 = p2;
            p2 = ptmp;
         }
         
         if(nanCounter_status == 1){
            ret[j] = *(p1+nh-1);//p1[nx-1]
         } else{
            ret[j] = mynan;
            counter_ea += 1;
         }
      }
      
   }// lower bounding if clause
   
   // Rcout << loc_threshold <<" ... "<< lb << " ... " << ret[j]<< "\n";
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW END        //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   // ret[0] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
   // ret[0] = cpp_dtw2vec_cm(cm, i, jsup, nh);
   loc_threshold = mymin(ret[0], threshold);
   
   // ret[0] can be NaN => in this case set bsfiw to Inf, for proper comparisons,
   // since comparisons with NaN are always FALSE
   double initial_bsfiw = mymin(ret[0], R_PosInf);
   
   // fill initially kNN_val
   if(kNNk > 0){
      kNN_inf = fill_kNN_Info(kNN_inf_list);
      initialize_kNN(kNN_inf, kNN_val, kNN_ix, kNN_lot_ix, lot_ix, kNNk,
                  initial_bsfiw, overlap_size);
   }
   
   bsfiw = initial_bsfiw; // best so far in window
   bsfiw_ix = 0; // last index of best so far in window
   
   
   // --------------------------------------start rolling----------------------------
   int j_prev_norm = 0;
   j = 1;
   while(j < (nx-nh+1)){
      // for(j=1; j < (nx-nh+1); j++){
      
      // for lot-mode: test if the current value should be calculated 
      while(ret[j] != ret[j] && j < (nx-nh+1)){
         j += 1;
      }
      
      jsup = j + nh; // the frontrunner index +1
      jmax = jsup -1; // the index of the frontrunner
      
      
      if(j_prev_norm == j-1){
         // update the discretization of the new window
         mu0 = mu;
         mu += (x[jmax] - x[j-1])/nh;
         sigma2 += (x[jmax] * x[jmax] - x[j-1] * x[j-1])/(nh-1) + 
            (mu0 * mu0 - mu * mu) * nh/(nh-1);
         
      }else{
         //  discretization of the new window from scratch
         mu0 = mu;
         mu = get_mean(x, j, jmax);
         mu2 = mu * mu;
         sigma2 = get_sigma2(x, j, jmax, mu2);
      }
      
      if(sigma2 < var_precision){
         sigma = 1; 
      } else{
         sigma = sqrt(sigma2);
      }
      j_prev_norm = j;
      
      
      if(debug == 1){
         Rcout << "j = "<< j << " ... mu: "<< mu << " "<< " sigma: "<< sigma << " sigma2: "<< sigma2 << "\n";
      }
      
      if((j - bsfiw_ix) >= overlap_size){
         bsfiw = threshold;
         bsfiw_ix = j;
      }
      
      
      
      //  bsfiw serves as threshold
      if(use_ea == 1){
         loc_threshold = mymin(bsfiw, threshold);
      }else{
         loc_threshold = threshold;
      }
      
      //  Lower bounding
      if(use_lb == 1){
         lb = get_lb_znorm(tube, x, mu, sigma, loc_threshold, j, jsup);
      }else{
         lb = -1;
      }
      
      if(kNNk > 0){
         // update local threshold, set it to the k-nearest neighbor
         // loc_threshold = mymin(loc_threshold, kNN_val[kNNk-1]);
         loc_threshold = mymin(loc_threshold, kNN_inf.vmax);
      }
      
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW BEGIN      //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      if(loc_threshold <= lb){
         //  loc_threshold < lb < DTW 
         ret[j] = mynan;
         counter_lb += 1;
         
      }else{
         
         nanCounter_status = 1;
         
         //initialize a and b with NAN
         for(int i=0; i < nh; i++){
            p1[i] = mynan;
            p2[i] = mynan;
         }
         
         // first column
         z_j = (x[j] - mu)/sigma;
         p1[0] = abs(z_j - h[0]);
         
         if(p1[0] > loc_threshold){
            ret[j] = mynan;
            counter_ea += 1;
            
         }else{
            
            iEnd   = std::min(nh, ws+1);
            for(int i=1; i < iEnd; i++){
               
               p1[i] = abs(z_j - h[i]) + p1[i-1];
               // p1[i] = cm(i, j) + p1[i-1];
               if(p1[i] > loc_threshold) p1[i] = mynan;
            }
            
            for(int k = j +1 ; k < jsup; k++){
               
               z_j = (x[k] - mu)/sigma;
               
               nanCounter = 0;
               iBegin = k-ws-j;
               if(iBegin <= 0){
                  // *p2 = cm(0, k) + *(p1);
                  *p2 = abs(z_j - h[0]) + *(p1);
                  if(*(p2) > loc_threshold){
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
               
               iEnd   = k+ws+1-j;
               if(iEnd >= nh){
                  iEnd = nh;
               }else{
                  *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
               }
               
               
               for (int i = iBegin; i < iEnd; i++){
                  *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), abs(z_j - h[i]));
                  if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                     *(p2+i) = mynan;
                     nanCounter ++;
                  }
               }
               
               if(nanCounter == nh) {
                  nanCounter_status = 0;
                  break;
               }
               ptmp=p1;
               p1 = p2;
               p2 = ptmp;
            }
            
            if(nanCounter_status == 1){
               ret[j] = *(p1+nh-1);//p1[nx-1]
               // Rcout << ret[j] << " ";
            } else{
               ret[j] = mynan;
               counter_ea += 1;
            }
         }
         
         
      }// lower bound if clause
      
      // Rcout << loc_threshold <<" ... "<< lb << " ... " << ret[j]<< "\n";
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW END        //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      
      if(kNNk > 0){
         // update set of k nearest neighbors if new value is smaller than 
         // k-best value so far
         
         if(ret[j] == ret[j]){ // if not nan
            if(j < kNN_inf.imax + overlap_size){
               // there is an overlap
               if(ret[j] < kNN_val[kNN_inf.which_imax]){
                  // because of overlap the new value needs to be smaller
                  // than the most recent best-sofar-value => so we replace
                  // this value with the new one
                  kick_imax_kNN_lot(kNN_val, kNN_ix, kNN_lot_ix, kNN_inf, ret[j], j, lot_ix);
                  
               }
            }else{
               // no overlap
               if(ret[j] < kNN_inf.vmax){
                  // new value needs to be at least smaller 
                  // than the farthest nearest neighbor
                  // update_kNN(kNN_val, kNN_ix, ret[j], j);   
                  kick_vmax_kNN_lot(kNN_val, kNN_ix, kNN_lot_ix, kNN_inf, ret[j], j, lot_ix);
               }
            }
         }
      }
      
      
      if(ret[j] < bsfiw){
         bsfiw = ret[j];
         bsfiw_ix = j;
      }
      
      j += 1;
   }//end rolling
   
   
   List list_ret;
   list_ret["dist"] = ret;
   IntegerVector counter(7);
   counter[0] = 0;
   counter[1] = 0;
   counter[2] = 0;
   counter[3] = 0;
   counter[4] = 0;
   counter[5] = counter_ea;
   counter[6] = counter_lb;
   list_ret["counter"] = counter;
   
   if(kNNk > 0){
      // final reverse minimum search of ret:
      // up to here the kNN values and indices dealt as additional 
      // local threshold, but since the possibility of protracting
      // a local minimum, a final reverse check is necessary:
      
      list_ret["all_best_indices"] = cpp_kNN_rev(ret, overlap_size, debug);
      
      kNN_inf_list["imax"]        = kNN_inf.imax;
      kNN_inf_list["nr_detected"] = kNN_inf.nr_detected;
      kNN_inf_list["nr_looking4"] = kNN_inf.nr_looking4;
      kNN_inf_list["vmax"]        = kNN_inf.vmax;
      kNN_inf_list["which_imax"]  = kNN_inf.which_imax;
      kNN_inf_list["which_vmax"]  = kNN_inf.which_vmax;
      list_ret["kNN_inf"] = kNN_inf_list;
   }
   delete[] p1;
   delete[] p2;
   
   return list_ret;
}



//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_rundtw_znorm_mv(const NumericMatrix &h,
                         const NumericMatrix &x,
                         std::string step_pattern, std::string dist_method,
                         List &kNN_inf_list,
                         int ws, double threshold, int overlap_tol = 0, 
                         int use_ea = 1, int use_lb = 1, int debug = 0) {
   
   // scan x with h and get the diff of DL(ss(x)) - DL(ss(x)-h) for each
   // subsegment ss(x)
   // h ... normalized
   // x ... not normalized
   // ws ... window size parameter for sakoe chiba window
   // threshold ... global threshold for early abandoning dtw calc
   //             (if threshold = Inf, then no thresholding takes place,
   //              in that case use_ea should be 0)
   // do_norm, use_ea, use_lb ... either 0 or 1
   
   // set dist_fun
   SEXP dist_SEXP = select_dist_mv_z(dist_method);
   XPtr<funcPtr_dist_mv_z> xpfun_dist_mv_z(dist_SEXP);
   funcPtr_dist_mv_z dist_fun = *xpfun_dist_mv_z;
   
   // set step function
   SEXP step_SEXP = selectVecStepRun(step_pattern);
   XPtr<funcPtr_step_run> xpfun_step(step_SEXP);
   funcPtr_step_run step_fun = *xpfun_step;
   
   // set lb function
   SEXP lb_SEXP = select_lb_z(dist_method);
   XPtr<funcPtr_lb_z> xpfun_lb_z(lb_SEXP);
   funcPtr_lb_z lb_fun = *xpfun_lb_z;
   
   if(use_lb == 1){
      use_ea = 1; // lb makes no sense if the loc_threshold is not set by bsfiw
   }
   
   if(use_ea == 0){
      threshold = R_PosInf;
   }
   
   int nh = h.nrow();
   int overlap_size = nh - overlap_tol;
   int nx = x.nrow();
   int nc = x.ncol();
   int k, i, j,c, jmax, jsup, bsfiw_ix;
   int counter_lb = 0;
   int counter_ea = 0;
   
   double bsfiw, loc_threshold; // local threshold
   double lb = 0;
   double mu0 = 0;
   
   int kNNk = as<int>(kNN_inf_list["nr_looking4"]);
   IntegerVector kNN_ix (kNNk);
   NumericVector kNN_val (kNNk, 0.0);
   kNN_Info kNN_inf;
   
   double var_precision = 0.000000001 ; // 1e-9
   NumericVector sigma2(nc); 
   NumericVector sigma(nc); 
   NumericVector mu(nc); 
   NumericVector ret(nx-nh+1, 0.0); // return vetor
   
   
   // set initial maximum indices
   j = 0;
   while(ret[j] != ret[j] && j < (nx-nh+1)){
      // if ret[j] is NaN, then ret[j] != ret[j] => these seubsegments 
      // are artifificial and emerge due to time series concatenation in
      // lot_mode
      j += 1;
   }
   jsup = j + nh;
   jmax = jsup -1;
   
   set_mean_sigma_mv(mu, sigma, sigma2, x, j, jmax, nc);
   if(debug == 1){
      for(int c = 0; c<nc; c++){
         Rcout << "j = "<< j << " ... mu: "<< mu[c] << " "<< " sigma: "<< sigma[c] << " sigma2: "<< sigma2[c] << "\n";   
      }
   }
   
   // get initial dtw value
   loc_threshold = threshold;
   
   NumericMatrix tube(nh, 2*nc);
   if(use_lb == 1){
      cpp_set_tube_mv(tube, h, ws);
      lb = lb_fun(tube, x, mu, sigma, loc_threshold, j, jsup, nc);
      
   }else{
      lb = -1;
   }
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW BEGIN      //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   int iBegin = 0;
   int iEnd = 0;
   int nanCounter_status = 1;
   
   double * p1 = new double [nh];
   double * p2 = new double [nh];
   double * ptmp;
   double mynan;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   
   
   if(loc_threshold <= lb){
      //  loc_threshold < lb < DTW 
      ret[j] = mynan;
      counter_lb += 1;
      
   }else{							   
      
      // fill the initial nh columns of the costmatrix
      // cpp_cm1_mv(cm, y, h, j, jsup, nh, nc);
      // cm_fun(cm, y, h, j, jsup, nh, nc);
      
      //initialize a and b with NAN
      for( i=0; i < nh; i++){
         p1[i] = mynan;
         p2[i] = mynan;
      }
      
      // first column
      // p1[0] = cm(0, j);
      p1[0] = dist_fun( h, x, mu, sigma, 0, j, nc);
      if(p1[0] > loc_threshold){
         ret[j] = mynan;
         
         counter_ea += 1;
      }else{
         
         iEnd   = std::min(nh, ws+1);
         for( i=1; i < iEnd; i++){
            p1[i] = dist_fun( h, x, mu, sigma, i, j, nc) + p1[i-1];
            // p1[i] = cm(i, j) + p1[i-1];
            if(p1[i] > loc_threshold) p1[i] = mynan;
         }
         
         for( k = j +1 ; k < jsup; k++){
            nanCounter = 0;
            iBegin = k-ws-j;
            if(iBegin <= 0){
               *p2 = dist_fun( h, x, mu, sigma, 0, k, nc) + *(p1);
               // *p2 = cm(0, k) + *(p1);
               if(*(p2) > loc_threshold){
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
            
            iEnd   = k+ws+1-j;
            if(iEnd >= nh){
               iEnd = nh;
            }else{
               *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
            }
            
            for ( i = iBegin; i < iEnd; i++){
               *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), dist_fun( h, x, mu, sigma, i, k, nc));
               // *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
               if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                  *(p2+i) = mynan;
                  nanCounter ++;
               }
            }
            
            if(nanCounter == nh) {
               nanCounter_status = 0;
               break;
            }
            ptmp=p1;
            p1 = p2;
            p2 = ptmp;
         }
         
         if(nanCounter_status == 1){
            ret[j] = *(p1+nh-1);//p1[nx-1]
         } else{
            ret[j] = mynan;
            counter_ea += 1;
            
         }
      }
   }// lower bounding if clause  
   
   
   
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW END        //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   // ret[0] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
   // ret[0] = cpp_dtw2vec_cm(cm, i, jsup, nh);
   loc_threshold = mymin(ret[0], threshold);
   
   //the first bsfiw can be NaN => set it to Inf  if necessary, for proper comparisons 
   double initial_bsfiw = mymin(ret[0], R_PosInf);
   
   bsfiw = initial_bsfiw; // best so far in window
   bsfiw_ix = 0; // last index of of best so far in window
   
   // fill initially kNN_val
   if(kNNk > 0){
      kNN_val[0] = initial_bsfiw;
      kNN_ix[0] = 0;
      
      for(k = 1; k < kNNk; k++){
         kNN_val[k] = R_PosInf;
         kNN_ix[k] = -99;
      }
      
      kNN_inf = fill_kNN_Info(kNN_inf_list);
   }
   
   
   // --------------------------------------start rolling----------------------------
   int j_prev_norm = 0;
   j = 1;
   while(j < (nx-nh+1)){
      // for(j=1; j < (nx-nh+1); j++){
      
      // for lot-mode: test if the current value should be calculated 
      while(ret[j] != ret[j] && j < (nx-nh+1)){
         j += 1;
      }
      
      jsup = j + nh; // the frontrunner index +1
      jmax = jsup -1; // the index of the frontrunner
      
      if(j_prev_norm == j-1){
         // update mu and sigma
         for( c = 0; c < nc; c++){
            mu0 = mu[c];
            mu[c] += ( x(jmax, c) - x(j-1, c) )/nh;
            sigma2[c] += (x(jmax, c) * x(jmax, c) - x(j-1, c) * x(j-1, c) )/(nh-1) + 
               (mu0 * mu0 - mu[c] * mu[c]) * nh/(nh-1);
            
            if(sigma2[c] < var_precision){
               sigma[c] = 1; 
            } else{
               sigma[c] = sqrt(sigma2[c]);
            }
         }
         
      }else{
         //  mu and sigma from scratch
         set_mean_sigma_mv(mu, sigma, sigma2, x, j, jmax, nc);
      }
      if(debug == 1){
         for(int c = 0; c<nc; c++){
            Rcout << "j = "<< j << " ... mu: "<< mu[c] << " "<< " sigma: "<< sigma[c] << " sigma2: "<< sigma2[c] << "\n";   
         }
      }
      j_prev_norm = j;
      
      
      if((j - bsfiw_ix) >= overlap_size){
         bsfiw = threshold;
         bsfiw_ix = j;
      }
      
      
      //  bsfiw serves as threshold
      if(use_ea == 1){
         loc_threshold = mymin(bsfiw, threshold);
      }else{
         loc_threshold = threshold;
      }
      
      //  Lower bounding
      if(use_lb == 1){
         lb = lb_fun(tube, x, mu, sigma, loc_threshold, j, jsup, nc);
      }else{
         lb = -1;
      }
      
      
      if(kNNk > 0){
         // update local threshold, set it to the k-nearest neighbor
         // loc_threshold = mymin(loc_threshold, kNN_val[kNNk-1]);
         loc_threshold = mymin(loc_threshold, kNN_inf.vmax);
      }
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW BEGIN      //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      
      if(loc_threshold <= lb){
         //  loc_threshold < lb < DTW 
         ret[j] = mynan;
         counter_lb += 1;
         
      }else{
         
         nanCounter_status = 1;
         
         //initialize a and b with NAN
         for( i=0; i < nh; i++){
            p1[i] = mynan;
            p2[i] = mynan;
         }
         
         // first column
         // p1[0] = cm(0, j);
         p1[0] = dist_fun( h, x, mu, sigma, 0, j, nc);
         if(p1[0] > loc_threshold){
            ret[j] = mynan;
            
            counter_ea += 1;
         }else{
            
            iEnd   = std::min(nh, ws+1);
            for( i=1; i < iEnd; i++){
               p1[i] = dist_fun( h, x, mu, sigma, i, j, nc) + p1[i-1];
               // p1[i] = cm(i, j) + p1[i-1];
               if(p1[i] > loc_threshold) p1[i] = mynan;
            }
            
            for( k = j +1 ; k < jsup; k++){
               nanCounter = 0;
               iBegin = k-ws-j;
               if(iBegin <= 0){
                  *p2 = dist_fun(h, x, mu, sigma, 0, k, nc) + *(p1);
                  // *p2 = cm(0, k) + *(p1);
                  if(*(p2) > loc_threshold){
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
               
               iEnd   = k+ws+1-j;
               if(iEnd >= nh){
                  iEnd = nh;
               }else{
                  *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
               }
               
               for ( i = iBegin; i < iEnd; i++){
                  *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), dist_fun(h, x, mu, sigma, i, k, nc));
                  // *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
                  if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                     *(p2+i) = mynan;
                     nanCounter ++;
                  }
               }
               
               if(nanCounter == nh) {
                  nanCounter_status = 0;
                  break;
               }
               ptmp=p1;
               p1 = p2;
               p2 = ptmp;
            }
            
            if(nanCounter_status == 1){
               ret[j] = *(p1+nh-1);//p1[nx-1]
               
            } else{
               ret[j] = mynan;
               counter_ea += 1;
               
            }
         }
      }// lower bound if clause
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW END        //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      
      
      if(kNNk > 0){
         // update set of k nearest neighbors if new value is smaller than 
         // k-best value so far
         
         debug_print_kNN(debug, "before", j, kNN_inf, kNN_val, kNN_ix);
         
         if(ret[j] == ret[j]){ // if not nan
            if(j < kNN_inf.imax + overlap_size){
               // there is an overlap
               if(ret[j] < kNN_val[kNN_inf.which_imax]){
                  // because of overlap the new value needs to be smaller
                  // than the most recent best-sofar-value => so we replace
                  // this value with the new one
                  kick_imax_kNN(kNN_val, kNN_ix, kNN_inf, ret[j], j);
                  
                  if(debug == 1){
                     Rcout << "---kick imax---\n";
                  }
               }
            }else{
               // no overlap
               if(ret[j] < kNN_inf.vmax){
                  // new value needs to be at least smaller 
                  // than the farthest nearest neighbor
                  // update_kNN(kNN_val, kNN_ix, ret[j], j);   
                  kick_vmax_kNN(kNN_val, kNN_ix, kNN_inf, ret[j], j);
                  
                  if(debug == 1){
                     Rcout << "---kick vmax---\n";
                  }
               }
            }
         }
         
         debug_print_kNN(debug, "after", j, kNN_inf, kNN_val, kNN_ix);
         
      }
      
      // ret[j] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
      // ret[i] = cpp_dtw2vec_cm(cm, j, jsup, nh);
      
      if(ret[j] < bsfiw){
         bsfiw = ret[j];
         bsfiw_ix = j;
      }
      
      j += 1;
   }// end rolling
   
   List list_ret;
   list_ret["dist"] = ret;
   IntegerVector counter(7);
   counter[0] = 0;
   counter[1] = 0;
   counter[2] = 0;
   counter[3] = 0;
   counter[4] = 0;
   counter[5] = counter_ea;
   counter[6] = counter_lb;
   
   list_ret["counter"] = counter;
   
   if(kNNk > 0){
      //TODO: final reverse minimum search of ret?!
      // list_ret["knn_values"] = kNN_val;
      // list_ret["knn_indices"] = kNN_ix;
      list_ret["all_best_indices"] = cpp_kNN_rev(ret, overlap_size, debug);
   }
   
   delete[] p1;
   delete[] p2;
   
   return list_ret;
}




//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_rundtw_znorm_mv_lot(const NumericMatrix &h,
                            const NumericMatrix &x,
                            NumericVector kNN_val_in,
                            IntegerVector kNN_ix_in,
                            IntegerVector kNN_lot_ix_in,
                            List &kNN_inf_list_in, 
                            int lot_ix,
                            std::string step_pattern, std::string dist_method,
                            int ws, double threshold, int overlap_tol = 0, 
                            int use_ea = 1, int use_lb = 1, int debug = 0) {
   
   // scan x with h and get the diff of DL(ss(x)) - DL(ss(x)-h) for each
   // subsegment ss(x)
   // h ... normalized
   // x ... not normalized
   // ws ... window size parameter for sakoe chiba window
   // threshold ... global threshold for early abandoning dtw calc
   //             (if threshold = Inf, then no thresholding takes place,
   //              in that case use_ea should be 0)
   // do_norm, use_ea, use_lb ... either 0 or 1
   
   NumericVector kNN_val    = clone(kNN_val_in);
   IntegerVector kNN_ix     = clone(kNN_ix_in);
   IntegerVector kNN_lot_ix = clone(kNN_lot_ix_in);
   List kNN_inf_list        = clone(kNN_inf_list_in);
   
   // set dist_fun
   SEXP dist_SEXP = select_dist_mv_z(dist_method);
   XPtr<funcPtr_dist_mv_z> xpfun_dist_mv_z(dist_SEXP);
   funcPtr_dist_mv_z dist_fun = *xpfun_dist_mv_z;
   
   // set step function
   SEXP step_SEXP = selectVecStepRun(step_pattern);
   XPtr<funcPtr_step_run> xpfun_step(step_SEXP);
   funcPtr_step_run step_fun = *xpfun_step;
   
   // set lb function
   SEXP lb_SEXP = select_lb_z(dist_method);
   XPtr<funcPtr_lb_z> xpfun_lb_z(lb_SEXP);
   funcPtr_lb_z lb_fun = *xpfun_lb_z;
   
   if(use_lb == 1){
      use_ea = 1; // lb makes no sense if the loc_threshold is not set by bsfiw
   }
   
   if(use_ea == 0){
      threshold = R_PosInf;
   }
   
   int nh = h.nrow();
   int overlap_size = nh - overlap_tol;
   int nx = x.nrow();
   int nc = x.ncol();
   int k, i, j,c, jmax, jsup, bsfiw_ix;
   int counter_lb = 0;
   int counter_ea = 0;
   
   double bsfiw, loc_threshold; // local threshold
   double lb = 0;
   double mu0 = 0;
   
   int kNNk = as<int>(kNN_inf_list["nr_looking4"]);
   kNN_Info kNN_inf;
   
   double var_precision = 0.000000001 ; // 1e-9
   NumericVector sigma2(nc); 
   NumericVector sigma(nc); 
   NumericVector mu(nc); 
   NumericVector ret(nx-nh+1, 0.0); // return vetor
   
   
   // set initial maximum indices
   j = 0;
   while(ret[j] != ret[j] && j < (nx-nh+1)){
      // if ret[j] is NaN, then ret[j] != ret[j] => these seubsegments 
      // are artifificial and emerge due to time series concatenation in
      // lot_mode
      j += 1;
   }
   jsup = j + nh;
   jmax = jsup -1;
   
   
   // set initial mu and sigma
   // for( c=0; c<nc; c++){
   //    tmp = 0;
   //    for( i = 0; i < nh; i++){
   //       tmp += x(i,c);
   //    }
   //    mu[c] = tmp/(nh);
   //    mu2 = mu[c] * mu[c];
   //    
   //    tmp = 0;
   //    for( i = 0; i < nh; i++){
   //       tmp += x(i,c) * x(i,c);
   //    }
   //    sigma2[c] = tmp/(nh-1) - mu2 * nh/(nh-1);
   //    sigma[c] = sqrt(sigma2[c]);
   // }
   
   set_mean_sigma_mv(mu, sigma, sigma2, x, j, jmax, nc);
   if(debug == 1){
      for(int c = 0; c<nc; c++){
         Rcout << "j = "<< j << " ... mu: "<< mu[c] << " "<< " sigma: "<< sigma[c] << " sigma2: "<< sigma2[c] << "\n";   
      }
   }
   
   // get initial dtw value
   loc_threshold = threshold;
   
   NumericMatrix tube(nh, 2*nc);
   if(use_lb == 1){
      cpp_set_tube_mv(tube, h, ws);
      lb = lb_fun(tube, x, mu, sigma, loc_threshold, j, jsup, nc);
      
   }else{
      lb = -1;
   }
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW BEGIN      //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   int iBegin = 0;
   int iEnd = 0;
   int nanCounter_status = 1;
   
   double * p1 = new double [nh];
   double * p2 = new double [nh];
   double * ptmp;
   double mynan;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   
   
   if(loc_threshold <= lb){
      //  loc_threshold < lb < DTW 
      ret[j] = mynan;
      counter_lb += 1;
      
   }else{							   
      
      // fill the initial nh columns of the costmatrix
      // cpp_cm1_mv(cm, y, h, j, jsup, nh, nc);
      // cm_fun(cm, y, h, j, jsup, nh, nc);
      
      //initialize a and b with NAN
      for( i=0; i < nh; i++){
         p1[i] = mynan;
         p2[i] = mynan;
      }
      
      // first column
      // p1[0] = cm(0, j);
      p1[0] = dist_fun( h, x, mu, sigma, 0, j, nc);
      if(p1[0] > loc_threshold){
         ret[j] = mynan;
         
         counter_ea += 1;
      }else{
         
         iEnd   = std::min(nh, ws+1);
         for( i=1; i < iEnd; i++){
            p1[i] = dist_fun( h, x, mu, sigma, i, j, nc) + p1[i-1];
            // p1[i] = cm(i, j) + p1[i-1];
            if(p1[i] > loc_threshold) p1[i] = mynan;
         }
         
         for( k = j +1 ; k < jsup; k++){
            nanCounter = 0;
            iBegin = k-ws-j;
            if(iBegin <= 0){
               *p2 = dist_fun( h, x, mu, sigma, 0, k, nc) + *(p1);
               // *p2 = cm(0, k) + *(p1);
               if(*(p2) > loc_threshold){
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
            
            iEnd   = k+ws+1-j;
            if(iEnd >= nh){
               iEnd = nh;
            }else{
               *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
            }
            
            for ( i = iBegin; i < iEnd; i++){
               *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), dist_fun( h, x, mu, sigma, i, k, nc));
               // *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
               if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                  *(p2+i) = mynan;
                  nanCounter ++;
               }
            }
            
            if(nanCounter == nh) {
               nanCounter_status = 0;
               break;
            }
            ptmp=p1;
            p1 = p2;
            p2 = ptmp;
         }
         
         if(nanCounter_status == 1){
            ret[j] = *(p1+nh-1);//p1[nx-1]
         } else{
            ret[j] = mynan;
            counter_ea += 1;
            
         }
      }
   }// lower bounding if clause  
   
   
   
   
   ////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////     1st DTW END        //////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   
   
   // ret[0] = cpp_dtw2vec_cm_ws_ea(cm, ws, loc_threshold, j, jsup, nh);
   // ret[0] = cpp_dtw2vec_cm(cm, i, jsup, nh);
   loc_threshold = mymin(ret[0], threshold);
   
   //the first bsfiw can be NaN => set it to Inf  if necessary, for proper comparisons 
   double initial_bsfiw = mymin(ret[0], R_PosInf);
   
   bsfiw = initial_bsfiw; // best so far in window
   bsfiw_ix = 0; // last index of of best so far in window
   
   // fill initially kNN_val
   if(kNNk > 0){
      kNN_inf = fill_kNN_Info(kNN_inf_list);
      initialize_kNN(kNN_inf, kNN_val, kNN_ix, kNN_lot_ix, lot_ix, kNNk,
                  initial_bsfiw, overlap_size);
   }
   
   
   // --------------------------------------start rolling----------------------------
   int j_prev_norm = 0;
   j = 1;
   while(j < (nx-nh+1)){
      // for(j=1; j < (nx-nh+1); j++){
      
      // for lot-mode: test if the current value should be calculated 
      while(ret[j] != ret[j] && j < (nx-nh+1)){
         j += 1;
      }
      
      jsup = j + nh; // the frontrunner index +1
      jmax = jsup -1; // the index of the frontrunner
      
      if(j_prev_norm == j-1){
         // update mu and sigma
         for( c = 0; c < nc; c++){
            mu0 = mu[c];
            mu[c] += ( x(jmax, c) - x(j-1, c) )/nh;
            sigma2[c] += (x(jmax, c) * x(jmax, c) - x(j-1, c) * x(j-1, c) )/(nh-1) + 
               (mu0 * mu0 - mu[c] * mu[c]) * nh/(nh-1);
            
            if(sigma2[c] < var_precision){
               sigma[c] = 1; 
            } else{
               sigma[c] = sqrt(sigma2[c]);
            }
         }
         
      }else{
         //  mu and sigma from scratch
         set_mean_sigma_mv(mu, sigma, sigma2, x, j, jmax, nc);
      }
      if(debug == 1){
         for(int c = 0; c<nc; c++){
            Rcout << "j = "<< j << " ... mu: "<< mu[c] << " "<< " sigma: "<< sigma[c] << " sigma2: "<< sigma2[c] << "\n";   
         }
      }
      j_prev_norm = j;
      
      
      if((j - bsfiw_ix) >= overlap_size){
         bsfiw = threshold;
         bsfiw_ix = j;
      }
      
      
      //  bsfiw serves as threshold
      if(use_ea == 1){
         loc_threshold = mymin(bsfiw, threshold);
      }else{
         loc_threshold = threshold;
      }
      
      //  Lower bounding
      if(use_lb == 1){
         lb = lb_fun(tube, x, mu, sigma, loc_threshold, j, jsup, nc);
      }else{
         lb = -1;
      }
      
      
      if(kNNk > 0){
         // update local threshold, set it to the k-nearest neighbor
         // loc_threshold = mymin(loc_threshold, kNN_val[kNNk-1]);
         loc_threshold = mymin(loc_threshold, kNN_inf.vmax);
      }
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW BEGIN      //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      
      if(loc_threshold <= lb){
         //  loc_threshold < lb < DTW 
         ret[j] = mynan;
         counter_lb += 1;
         
      }else{
         
         nanCounter_status = 1;
         
         //initialize a and b with NAN
         for( i=0; i < nh; i++){
            p1[i] = mynan;
            p2[i] = mynan;
         }
         
         // first column
         // p1[0] = cm(0, j);
         p1[0] = dist_fun( h, x, mu, sigma, 0, j, nc);
         if(p1[0] > loc_threshold){
            ret[j] = mynan;
            
            counter_ea += 1;
         }else{
            
            iEnd   = std::min(nh, ws+1);
            for( i=1; i < iEnd; i++){
               p1[i] = dist_fun( h, x, mu, sigma, i, j, nc) + p1[i-1];
               // p1[i] = cm(i, j) + p1[i-1];
               if(p1[i] > loc_threshold) p1[i] = mynan;
            }
            
            for( k = j +1 ; k < jsup; k++){
               nanCounter = 0;
               iBegin = k-ws-j;
               if(iBegin <= 0){
                  *p2 = dist_fun(h, x, mu, sigma, 0, k, nc) + *(p1);
                  // *p2 = cm(0, k) + *(p1);
                  if(*(p2) > loc_threshold){
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
               
               iEnd   = k+ws+1-j;
               if(iEnd >= nh){
                  iEnd = nh;
               }else{
                  *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
               }
               
               for ( i = iBegin; i < iEnd; i++){
                  *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), dist_fun(h, x, mu, sigma, i, k, nc));
                  // *(p2+i) =  step_fun(*(p2+i-1), *(p1+i-1), *(p1+i), cm(i,k));
                  if((*(p2+i) > loc_threshold) | (*(p2+i) != *(p2+i))){
                     *(p2+i) = mynan;
                     nanCounter ++;
                  }
               }
               
               if(nanCounter == nh) {
                  nanCounter_status = 0;
                  break;
               }
               ptmp=p1;
               p1 = p2;
               p2 = ptmp;
            }
            
            if(nanCounter_status == 1){
               ret[j] = *(p1+nh-1);//p1[nx-1]
               
            } else{
               ret[j] = mynan;
               counter_ea += 1;
               
            }
         }
      }// lower bound if clause
      
      ////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////     RUN DTW END        //////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      
      
      if(kNNk > 0){
         // update set of k nearest neighbors if new value is smaller than 
         // k-best value so far
         
         if(ret[j] == ret[j]){ // if not nan
            if(j < kNN_inf.imax + overlap_size){
               // there is an overlap
               if(ret[j] < kNN_val[kNN_inf.which_imax]){
                  // because of overlap the new value needs to be smaller
                  // than the most recent best-sofar-value => so we replace
                  // this value with the new one
                  kick_imax_kNN_lot(kNN_val, kNN_ix, kNN_lot_ix, kNN_inf, ret[j], j, lot_ix);
                  
               }
            }else{
               // no overlap
               if(ret[j] < kNN_inf.vmax){
                  // new value needs to be at least smaller 
                  // than the farthest nearest neighbor
                  // update_kNN(kNN_val, kNN_ix, ret[j], j);   
                  kick_vmax_kNN_lot(kNN_val, kNN_ix, kNN_lot_ix, kNN_inf, ret[j], j, lot_ix);
               }
            }
         }
      }
      
      if(ret[j] < bsfiw){
         bsfiw = ret[j];
         bsfiw_ix = j;
      }
      
      j += 1;
   }// end rolling
   
   List list_ret;
   list_ret["dist"] = ret;
   IntegerVector counter(7);
   counter[0] = 0;
   counter[1] = 0;
   counter[2] = 0;
   counter[3] = 0;
   counter[4] = 0;
   counter[5] = counter_ea;
   counter[6] = counter_lb;
   
   list_ret["counter"] = counter;
   
   if(kNNk > 0){
      // final reverse minimum search of ret:
      // up to here the kNN values and indices dealt as additional 
      // local threshold, but since the possibility of protracting
      // a local minimum, a final reverse check is necessary:
      
      list_ret["all_best_indices"] = cpp_kNN_rev(ret, overlap_size, debug);
      
      kNN_inf_list["imax"]        = kNN_inf.imax;
      kNN_inf_list["nr_detected"] = kNN_inf.nr_detected;
      kNN_inf_list["nr_looking4"] = kNN_inf.nr_looking4;
      kNN_inf_list["vmax"]        = kNN_inf.vmax;
      kNN_inf_list["which_imax"]  = kNN_inf.which_imax;
      kNN_inf_list["which_vmax"]  = kNN_inf.which_vmax;
      list_ret["kNN_inf"] = kNN_inf_list;
   }
   
   delete[] p1;
   delete[] p2;
   
   return list_ret;
}




//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double Rcode(double a){
   return a * 1;
}






// lessons learned: 
// cpp_rnudtw_mv does NOT get faster if
//       + I copy paste the code of norm, range in the function directly
//         instead of using function calls or using "inline"
//       + using the function names of cpp_cm1_mv and mystep_symmetric1
//         directly instead of function pointer


/*** R

# test the kNN
if(FALSE){
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- NULL
   noise <- function(nr) cumsum(rnorm(nr))
   goal <- function(hnorm, x, sp, ws = 10, kNNk){
      allv <- cpp_rundtw(hnorm, x, sp, ws = ws, threshold = Inf, 
                         kNNk = 0, do_norm = 1, use_ea = 0, use_lb = 0)
      dd <- allv$dist
      best_i <- integer()
      for(i in 1:kNNk){
         ix <- which.min(dd)
         best_i <- c(best_i, ix-1)
         dd[max(0,(ix-nh)):min((ix+nh-1), length(dd))] <- Inf
      }
      best_i
   }
   
   
   nx <- 500
   nh <- 20
   nfits <- 10
   
   nn <- nx - nfits * nh# nnoise
   nn <- nn/nfits
   
   h <- noise(nh)
   hnorm <- IncDTW::norm( h , type="01")
   x <- noise(0)
   for(i in 1:nfits){
      x <- c(x, noise(nn), h)
   }
   
   # rundtw(hnorm, x,dm,  sp, k = 5, ws = 10, threshold = Inf,1,1,1)# DO norm
   # cpp_rundtw(hnorm, x, sp, ws = 10, threshold = Inf, kNNk = 10, do_norm = 1, use_ea = 1, use_lb = 1)
   
   k <- nfits + 2
   
   
   tmp <- cpp_rundtw(hnorm, x, sp, ws = 10, threshold = Inf, kNNk = k,
                     do_norm = 1, use_ea = 1, use_lb = 1, debug = 0)
   sort(tmp$knn_indices)
   
   sort( goal(hnorm, x, sp, ws=10, kNNk = k) )
   
   
   zz <- file("build_ignore/mySink.Rout", open = "wt"); sink(zz); sink(zz, type = "message")
   cpp_rundtw(hnorm, x, sp, ws = 10, threshold = Inf, kNNk = 3,
              do_norm = 1, use_ea = 1, use_lb = 1, debug = 1)
   
   sink();close(zz);
   closeAllConnections()
   
}


# test the synthetic kNN, where samples are concatenated
if(FALSE){
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- NULL
   noise <- function(nr) cumsum(rnorm(nr))
   goal <- function(hnorm, x, sp, ws = 10, kNNk){
      allv <- cpp_rundtw(hnorm, x, sp, ws = ws, threshold = Inf, 
                         kNNk = 0, do_norm = 1, use_ea = 0, use_lb = 0)
      dd <- allv$dist
      best_i <- integer()
      for(i in 1:kNNk){
         ix <- which.min(dd)
         best_i <- c(best_i, ix-1)
         dd[max(0,(ix-nh)):min((ix+nh-1), length(dd))] <- Inf
      }
      best_i
   }
   
   
   nx <- 50
   nh <- 20
   nfits <- 10
   h <- noise(nh)
   x <- c(h, noise(10), Inf,noise(10), h)
   hnorm <- IncDTW::norm( h , type="01")
   
   k <- nfits + 2
   
   
   tmp <- cpp_rundtw(hnorm, x, sp, ws = 10, threshold = Inf, kNNk = k,
                     do_norm = 1, use_ea = 1, use_lb = 1, debug = 0)
   sort(tmp$knn_indices)
   
   sort( goal(hnorm, x, sp, ws=10, kNNk = k) )
   
   
   zz <- file("build_ignore/mySink.Rout", open = "wt"); sink(zz); sink(zz, type = "message")
   cpp_rundtw(hnorm, x, sp, ws = 10, threshold = Inf, kNNk = 3,
              do_norm = 1, use_ea = 1, use_lb = 1, debug = 1)
   
   sink();close(zz);
   closeAllConnections()
   
}


# test lower bounding and norming multivar
if(FALSE){
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- NULL
   noise <- function(nr, nc) matrix(cumsum(rnorm(nr*nc)), nrow =nr, ncol=nc)
   # deform <- function(x) (x + rnorm(1, 0, 100)) * abs(rnorm(1, 0, 100))
   
   nx <- 500
   nh <- 50
   nc <- 2
   nfits <- 5
   
   nn <- nx - nfits * nh# nnoise
   nn <- nn/nfits
   
   h <- noise(nh, nc)
   hnorm <- IncDTW::norm( h , type="01")
   x <- matrix(numeric(), ncol=nc)
   for(i in 1:nfits){
      x <- rbind(x, noise(nn, nc), h)
   }
   
   cpp_rundtw_mv(hnorm, x, sp, dm, 20, Inf, 1, 1, 1)
   cpp_rundtw_mv(h, x, sp, dm, 20, Inf, 0, 1, 1)
   
   # no window size parameter
   cpp_rundtw_mv(hnorm, x, sp, dm, nh, Inf, 1, 1, 1)
   cpp_rundtw_mv(h, x, sp, dm, nh, Inf, 0, 1, 1)
}


#test lower bounding univ
if(FALSE){
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- NULL
   noise <- function(nr) cumsum(rnorm(nr))
   
   nx <- 500
   nh <- 30
   nfits <- 5
   
   nn <- nx - nfits * nh# nnoise
   nn <- nn/nfits
   
   h <- noise(nh)
   hnorm <- IncDTW::norm( h , type="01")
   x <- noise(0)
   for(i in 1:nfits){
      x <- c(x, noise(nn), h)
   }
   
   cpp_rundtw(hnorm, x, sp, ws = 10, threshold = Inf,1,1,1)# DO norm
   cpp_rundtw(h, x, sp, ws = 10, threshold = Inf,0,1,1)# do NOT norm
}


if(FALSE){
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- NULL
   noise <- function(nr, nc) matrix(cumsum(rnorm(nr*nc)), nrow =nr, ncol=nc)
   # deform <- function(x) (x + rnorm(1, 0, 100)) * abs(rnorm(1, 0, 100))
   
   # speed comparison
   foo_2vec <- function(hnorm, x, ws){
      nh <- nrow(hnorm)
      nx <- nrow(x)
      sapply(1:(nx-nh+1), function(i){
         y <- IncDTW::norm(x[i:(i+nh-1), , drop=F], type = "01")
         IncDTW::dtw2vec(y, hnorm, dist_method = dm,
                         step_pattern = sp, ws = ws)$dist
      })
   }
   
   foo_cm <- function(hnorm, x, ws){
      nh <- nrow(hnorm)
      nx <- nrow(x)
      sapply(1:(nx-nh+1), function(i){
         y <- IncDTW::norm(x[i:(i+nh-1), , drop=F], type = "01")
         cm_tmp <- IncDTW::cm(y, hnorm, dist_method = dm, ws = ws)
         IncDTW::dtw2vec(cm_tmp, C = "cm", step_pattern = sp, ws = ws)$dist
      })
   }
   
   nx <- 500
   nh <- 30
   nfits <- 5
   
   nn <- nx - nfits * nh# nnoise
   nn <- nn/nfits
   
   h <- noise(nh, nc)
   hnorm <- IncDTW::norm( h , type="01")
   x <- matrix(numeric(), ncol=nc)
   for(i in 1:nfits){
      x <- rbind(x, noise(nn, nc), h)
   }
   
   
   microbenchmark::microbenchmark(cpp_rundtw_mv(h, x,sp, dm, 100, Inf, 1),
                                  cpp_rundtw_mv_test(h, x, 100, Inf, 1),
                                  times = 30)
}
*/
