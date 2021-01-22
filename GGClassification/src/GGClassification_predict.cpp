// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <vector>
#include <cmath>

using std::vector;
using Rcpp::as;
using Rcpp::List;
using Rcpp::wrap;
using Rcpp::NumericMatrix;
using Rcpp::IntegerVector;
using Eigen::ArrayXXd;
using Eigen::ArrayXd;
using Eigen::ArrayXi;

typedef Eigen::Map<ArrayXXd> MapArrayXXd;

//IntegerVector

// [[Rcpp::export]]
IntegerVector predict(List model, NumericMatrix& X){

  // Receives the model created by model function and an unlabeled matrix of data.
  // Then it classifies the unlabeled data using the model parameters.

  MapArrayXXd X_array(as<MapArrayXXd>(X)); // Mapping R variables to C++/Eigen variables types.

  ArrayXXd array_midpoints = model["Midpoints"];
  ArrayXXd array_w = model["W"];
  vector<double> vector_bias =  as<vector<double> >(model["Bias"]);
  vector<int> labels = as<vector<int> >(model["Labels"]);

  int nrows = X_array.rows(); // Getting the dimensions of data
  int ncols = X_array.cols();
  int number_of_midpoints = vector_bias.size();

  ArrayXd decision_maker(nrows); // More auxiliary variables
  ArrayXd vector_class(nrows);
  ArrayXXd array_class(nrows, number_of_midpoints);
  ArrayXXd array_class_abs(nrows, number_of_midpoints);
  ArrayXXd array_distances(nrows, number_of_midpoints);
  ArrayXd vector_sigma_squared(nrows);
  ArrayXi vector_y(nrows);
  ArrayXd output(number_of_midpoints);

  // Here begins the settings of classification.
  for(int i=0; i<number_of_midpoints; i++){
    decision_maker = (X_array.rowwise()*array_w.row(i)).rowwise().sum();
    decision_maker = (decision_maker-vector_bias[i]);

    for(int j=0; j<nrows; j++){
      if( tanh(decision_maker(j)) > 0){
        vector_class(j) = 1;
      }
      else{
        vector_class(j) = -1;
      }
    }

    array_distances.col(i) = (X_array.rowwise() - array_midpoints.row(i)).matrix().rowwise().lpNorm<2>(); //
    array_class.col(i) = (vector_class*array_distances.col(i));
  }


  array_class_abs = array_class.abs();

  vector_sigma_squared = (array_distances.rowwise().maxCoeff()).square();

  for(int m=0; m<nrows; m++){

    array_class_abs.row(m) = 1/array_class_abs.row(m);
    array_class_abs.row(m) = -vector_sigma_squared(m)*array_class_abs.row(m);
    array_class_abs.row(m) = array_class_abs.row(m).exp();
    array_class_abs.row(m) = (1/array_class_abs.row(m) + 1e-16);
    array_class_abs.row(m) = array_class_abs.row(m)/(array_class_abs.row(m)).sum();
    output = array_class_abs.row(m)*(array_class.row(m).sign());

    if( output.sum() >= 0 ){ // Classification
      vector_y(m) = labels[1]; // Class 2
    }
    else{
      vector_y(m) = labels[0]; // Class 1
    }
  }

  return(wrap(vector_y));
}
