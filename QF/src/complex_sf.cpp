#include <complex>
extern "C"
{
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
}


// [[Rcpp::export]]
std::complex<double> pochhammer_complex(std::complex<double> x, std::complex<double> y){

  gsl_sf_result r1;
  gsl_sf_result i1;
  gsl_sf_result r2;
  gsl_sf_result i2;

  gsl_sf_lngamma_complex_e(real(x+y), imag(x+y), &r1, &i1);
  gsl_sf_lngamma_complex_e(real(x), imag(x), &r2, &i2);

  std::complex<double> ln_g_1(r1.val,i1.val);
  std::complex<double> ln_g_2(r2.val,i2.val);

  std::complex<double> pochhammer_def = std::exp(ln_g_1-ln_g_2);

  return(pochhammer_def);
}




// [[Rcpp::export]]
std::complex<double> beta_complex(std::complex<double> a, std::complex<double> b){

  gsl_sf_result r1;
  gsl_sf_result i1;
  gsl_sf_result r2;
  gsl_sf_result i2;
  gsl_sf_result r3;
  gsl_sf_result i3;

  gsl_sf_lngamma_complex_e(real(a), imag(a), &r1, &i1);
  gsl_sf_lngamma_complex_e(real(b), imag(b), &r2, &i2);
  gsl_sf_lngamma_complex_e(real(a+b), imag(a+b), &r3, &i3);

  std::complex<double> ln_g_1(r1.val,i1.val);
  std::complex<double> ln_g_2(r2.val,i2.val);
  std::complex<double> ln_g_3(r3.val,i3.val);

  std::complex<double> beta_def=std::exp(ln_g_1+ln_g_2-ln_g_3);

  return(beta_def);
}


// [[Rcpp::export]]
std::complex<double> gamma_complex(std::complex<double> n){

  gsl_sf_result r1;
  gsl_sf_result i1;

  gsl_sf_lngamma_complex_e(real(n), imag(n), &r1, &i1);

  std::complex<double> ln_g_1(r1.val,i1.val);

  std::complex<double> gamma_def=std::exp(ln_g_1);

  return(gamma_def);
}



