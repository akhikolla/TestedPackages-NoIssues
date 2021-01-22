/*
Grøtan Grotan and Steinar Engen 2007
Masatoshi Katabuchi 2020
*/

#include <RcppNumerical.h>
#include <RcppEigen.h>
using namespace Numer;

double maxf(int x, double mu, double sig)
{
   double d,z;
   z = 0;
   d = 100;
   while (d > 0.00001) {
     if (x - 1 - exp(z) - 1 / sig * (z - mu)  > 0) z = z + d; else z = z - d;
     d = d / 2;
   }
   return(z);
 };

double upper(int x, double m, double mu, double sig)
{
   double d, z, mf;
   mf = (x - 1) * m - exp(m) - 0.5 / sig * ((m - mu) * (m - mu));
   z = m + 20;
   d = 10;
   while (d > 0.000001) {
      if ((x - 1) * z - exp(z) - 0.5 / sig * ((z - mu) * (z - mu)) - mf 
          + log(1000000.0) > 0) 
        z = z + d; else z = z - d;
      d = d / 2;
   }
   return(z);
 };

/* ---------------------------------------------------------------------------*/

double lower(int x, double m, double mu, double sig)
{
   double d, z, mf;
   mf = (x - 1) * m - exp(m) - 0.5 / sig * ((m - mu) * (m - mu));
   z = m - 20;
   d = 10;
   while (d > 0.000001) {
      if ((x - 1) * z - exp(z) - 0.5 / sig * ((z - mu) * (z - mu)) - mf 
          + log(1000000.0) > 0) 
        z = z - d; else z = z + d;
      d = d / 2;
   }
   return(z);
 };


class plnintegrand: public Func
{
private:
  int x;
  double mu;
  double sig;
public:
  plnintegrand(int x_, double mu_, double sig_) : x(x_), mu(mu_), sig(sig_){}
  double operator()(const double& z) const {
    double fac;
    fac = std::lgamma(x + 1);
    return exp(z * x - exp(z) - 0.5 / sig *((z - mu) * (z - mu)) - fac);
  }
};

class plnintegrand2: public Func
{
private:
  int x;
  double mu;
  double sig;
public:
  plnintegrand2(int x_, double mu_, double sig_) : x(x_), mu(mu_), sig(sig_){}
  double operator()(const double& z) const {
    double fac, z2;
    fac = std::lgamma(x + 1);
    z2 = exp(z);
    return exp(z * x - log(exp(z2) - 1) - 0.5 / sig *((z - mu) * (z - mu)) - fac);
  }
};


Rcpp::NumericVector do_dpln(Rcpp::IntegerVector x, double mu, double sig){
  int n = x.size();
  Rcpp::NumericVector out(n);
  double a, b, m;
  for (int i = 0; i < n; i++) {
    m = maxf(x[i], mu, sig);
    a = lower(x[i], m, mu, sig);
    b = upper(x[i], m, mu, sig);
    plnintegrand f(x[i], mu, sig);
    double err_est;
    int err_code;
    out[i] = integrate(f, a, b, err_est, err_code) 
      * (1 / std::sqrt(2 * M_PI * sig));
  }
  return out;
}


double check_diff(double mu, double sig) {
  double a, b, m, a0, b0, m0;
  m = maxf(100, mu, sig);
  a = lower(100, m, mu, sig);
  b = upper(100, m, mu, sig);
  m0 = maxf(0, mu, sig);
  a0 = lower(0, m0, mu, sig);
  b0 = upper(0, m0, mu, sig);
  plnintegrand f100(100, mu, sig);
  plnintegrand f0(0, mu, sig);
  plnintegrand2 f100_2(100, mu, sig);
  double err_est;
  int err_code;
  double out100_2 = integrate(f100_2, a, b, err_est, err_code) 
    * (1 / std::sqrt(2 * M_PI * sig)); 
  double out100 = integrate(f100, a, b, err_est, err_code) 
    / (std::sqrt(2 * M_PI * sig) - integrate(f0, a0, b0, err_est, err_code)); 
  
  return out100_2 / out100;
}

Rcpp::NumericVector do_dpln2(Rcpp::IntegerVector x, double mu, double sig){
  int n = x.size();
  Rcpp::NumericVector out(n);
  double a, b, m, a0, b0, m0;
  double diff = check_diff(mu, sig);
  for (int i = 0; i < n; i++) {
    m = maxf(x[i], mu, sig);
    a = lower(x[i], m, mu, sig);
    b = upper(x[i], m, mu, sig);
    /*b needs to be smaller than log(709)*/
    if (b <= 6.563856) {
      plnintegrand2 f2(x[i], mu, sig);
      double err_est;
      int err_code;
      out[i] = integrate(f2, a, b, err_est, err_code) 
        * (1 / std::sqrt(2 * M_PI * sig)); 
    } else {
      double err_est;
      int err_code;
      m0 = maxf(0, mu, sig);
      a0 = lower(0, m0, mu, sig);
      b0 = upper(0, m0, mu, sig);
      plnintegrand f(x[i], mu, sig);
      plnintegrand f0(0, mu, sig);
      double tmp =  integrate(f, a, b, err_est, err_code);
      double tmp0 =  integrate(f0, a0, b0, err_est, err_code);
      out[i] = tmp / (std::sqrt(2 * M_PI * sig) - tmp0) * diff;
    }
  }
  return out;
}


// [[Rcpp::export]]
Rcpp::NumericVector do_dztpln(Rcpp::IntegerVector x, double mu, double sig){
  int n = x.size();
  double sig2 = sig * sig;
  Rcpp::IntegerVector zero (1);
  Rcpp::NumericVector p(n), p00(1), p0(n), lik(n); 
  p = do_dpln(x, mu, sig2);
  p00 = do_dpln(zero, mu, sig2);
  //p00 = do_dpln(0, mu, sig2);
  p0 = Rcpp::rep(p00, n);
  lik = p / (1.0 - p0);
  return lik;
}

// [[Rcpp::export]]
Rcpp::NumericVector do_dztpln2(Rcpp::IntegerVector x, double mu, double sig){
  double sig2 = sig * sig;
  return do_dpln2(x, mu, sig2);
}

// [[Rcpp::export]]
Eigen::VectorXd do_dztplnm(Rcpp::IntegerVector x, Rcpp::NumericVector mu,
    Rcpp::NumericVector sigma, Eigen::VectorXd theta){
    Rcpp::NumericMatrix lik(x.size(), mu.size());
     for ( int j = 0; j < mu.size(); j++ ) {
       lik.column(j) = do_dztpln(x, mu(j), sigma(j));;
    }
  Eigen::Map<Eigen::MatrixXd> lik2(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(lik));
  return  lik2 * theta;
}


// [[Rcpp::export]]
Eigen::VectorXd do_dztplnm2(Rcpp::IntegerVector x, Rcpp::NumericVector mu,
    Rcpp::NumericVector sigma, Eigen::VectorXd theta){
    Rcpp::NumericMatrix lik(x.size(), mu.size());
     for ( int j = 0; j < mu.size(); j++ ) {
       lik.column(j) = do_dztpln2(x, mu(j), sigma(j));;
    }
  Eigen::Map<Eigen::MatrixXd> lik2(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(lik));
  return  lik2 * theta;
}
