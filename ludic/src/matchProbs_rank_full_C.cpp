#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//'Compute the matching probabilities for each pair of observations
//'
//'C++ version: for each observations in \code{(1:n)}, all the matching probabilities are computed
//'for the \code{p} possible pairs.
//'
//'@param computed_dist a \code{n x p} matrix of computed distances used for ranking.
//'@param prop_match a priori proportion of matches ("rho_1")
//'
//'@return a \code{n x p} matrix containing the matching probabilities for each pair
//'
//'@export
// [[Rcpp::export]]
NumericMatrix matchProbs_rank_full_C(NumericMatrix computed_dist, 
                                     double prop_match){
  
  mat compD = as<mat>(computed_dist);
  int p = compD.n_cols;
  int n = compD.n_rows;
  
  vec logexptrick_const = max(compD, 1);
  mat prob0 = mat(n, p);
  
  
  for(int i=0; i<n; i++){
    double logexptrick_const_rowi = logexptrick_const(i);
    rowvec compD_rowi = compD.row(i);
    double normconst = exp(-logexptrick_const_rowi) + sum(exp(compD_rowi - logexptrick_const_rowi + log(prop_match))); //exponential trick
    prob0.row(i) = exp(compD_rowi + log(prop_match) - logexptrick_const_rowi) / normconst;
  //  exp(computed_dist[id0,] + log(prop_match) - logexptrick_const)/(exp(-logexptrick_const) + sum(exp(computed_dist[id0,] - logexptrick_const + log(prop_match))))
  }
  return(wrap(prob0));
}

//rank_match_1 <- t(sapply(rownames(dist_all), FUN=matchProbs_rank_complete, computed_dist=dist_all, 
//                         prop_match=prop_match_1way))
