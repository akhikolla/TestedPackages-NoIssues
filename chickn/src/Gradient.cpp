// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <iostream>

using namespace arma;
using namespace std;
using namespace Rcpp;


//' @title G_fun_cpp
//' @description Function for gradient computation. 
//' @param x is a vector
//' @param y is a vector
//' @param W is a frequency matrix
//' @return A vector
//' @keywords internal
//' @export
// [[Rcpp::export]]
arma::rowvec G_fun_cpp(arma::vec x, arma::vec y, arma::mat W){
  int m = W.n_rows;
  arma::vec sk_y(m);
  for(int i = 0; i<m; i++){
    sk_y[i] = x[i]*y[m+i] - x[m+i]*y[i];
  }
  return (sk_y.t()*W);
}
//'@title Objective function 
//'@description Objective function of the minimization with respect to the cluster centroids
//'
//'@param c is a cluster centroid 1xn
//'@param W is a frequency matrix mxn
//'@param residue is a vector
//'@return The objective function value 
//'@details The residue vector is equal to \eqn{(Sketch(Data, W) - \sum_{k=1}^K \alpha_k *Sketch(c_k, W))}. 
//'@export
//'@keywords internal
// [[Rcpp::export]]
double ObjFun_OMP_cpp(arma::vec c, arma::mat W, arma::vec residue){
  int m = W.n_rows;
  arma::vec prodWc = W*c;
  arma::vec SK_c(2*m);
  SK_c.subvec(0, m-1) = cos(prodWc);
  SK_c.subvec(m, 2*m-1) = sin(prodWc);
  double norm_SK_c = as_scalar(sqrt(SK_c.t()*SK_c));
  SK_c = SK_c/norm_SK_c;
  return as_scalar(-residue.t()*SK_c);
}
//' @title Gradient OMP
//'@description Gradient with respect to cluster centroid vector
//'
//'@param c is a cluster centroid
//'@param W is a frequency matrix mxs
//'@param residue is a vector  
//'@return The gradient vector
//'@keywords internal 
//'@export
// [[Rcpp::export]]
arma::rowvec Gradient_OMP_cpp(arma::vec c, arma::mat W, arma::vec residue){
  int m = W.n_rows;
  arma::vec prodWc = W*c;
  arma::vec SK_c(2*m);
  SK_c.subvec(0, m-1) = cos(prodWc);
  SK_c.subvec(m, 2*m-1) = sin(prodWc);
  //vec SK_c =  join_cols(cos(prodWc), sin(prodWc));
  double norm_SK_c = as_scalar(sqrt(SK_c.t()*SK_c));
  SK_c = SK_c/norm_SK_c;
  double val = as_scalar(-residue.t()*SK_c);
  vec y = val*SK_c + residue;
  rowvec G  = G_fun_cpp(SK_c, y, W);
  return -G/norm_SK_c;
}

//'@title Gradient and objective function
//'@description Global objective function and gradient computations with respect to cluster centroid vectors and their weights.
//'
//'@param x is a data vector. Its first K*n components are cluster centroids and its last K components are the centroid weights. 
//'@param SK is a data sketch vector. 
//'@param W is a frequency matrix.
//'@param K is a number of cluster centroids. 
//'@return \itemize{
//'           \item \code{gradient} is a gradient vector
//'           \item \code{objective} is an objective function value \eqn{\|SK - \sum_{k=1}^K \alpha_k Sketch(c_k, W)\|}} 
//' @keywords internal
//'@export
// [[Rcpp::export]]
List Gradient_cpp(arma::rowvec x, arma::vec SK, arma::mat W, int K) { 
  int m = W.n_rows;
  int d = W.n_cols;
  
  arma::mat C(x.subvec(0,K*d-1));
  C.reshape(d, K);
  arma::rowvec A(x.subvec(K*d, K*(d+1)-1));
  
  arma::mat prodWC = W*C;
  
  arma::mat P = join_cols(cos(prodWC), sin(prodWC));
  arma::mat phial(P);
  arma::vec diff(2*m);
  
   for(int i = 0; i< 2*m;i++){
     phial.row(i) = P.row(i)%A;
     diff(i) = sum(phial.row(i)) -SK(i);
   }
  
  List res;
  res["objective"] = diff.t()*diff;
  rowvec grad(d*K);
  for(int i = 0; i< K; i++){
  grad.subvec(i*d, (i+1)*d -1) = G_fun_cpp(phial.col(i), diff, W);
  }
  rowvec grad_weight(K);
  grad_weight = diff.t()*P;
  res["gradient"] = 2*join_rows(grad, grad_weight);
  return res;
}
