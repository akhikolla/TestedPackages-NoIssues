#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

std::vector<double> rQF_c(int n, std::vector<double> lambdas, std::vector<double> etas) {
  int r = lambdas.size();
  std::vector<double> sample(n,0.0);
  for(int i=0; i<n; i++){
    for(int j=0; j<r; j++){
      sample[i] += lambdas[j]*R::rnchisq(1.0, etas[j]);
    }
  }
  return(sample);
}

// [[Rcpp::export]]

std::vector<double> rQF_ratio_c(int n, std::vector<double> lambdas_1, std::vector<double> lambdas_2,
                              std::vector<double> etas_1, std::vector<double> etas_2) {
  int r_1 = lambdas_1.size(), r_2 = lambdas_2.size();
  double sample_num, sample_den;
  std::vector<double> sample(n);

  for(int i=0; i<n; i++){
    sample_num = 0.0;
    sample_den = 0.0;
    for(int j=0; j<r_1; j++){
      sample_num += lambdas_1[j]*R::rnchisq(1.0, etas_1[j]);
    }
    for(int j=0; j<r_2; j++){
      sample_den += lambdas_2[j]*R::rnchisq(1.0, etas_2[j]);
    }
    sample[i] = sample_num / sample_den;
  }
  return(sample);
}


