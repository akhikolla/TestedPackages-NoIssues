// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <limits>

using std::numeric_limits;
using Rcpp::as;
using Rcpp::wrap;
using Rcpp::NumericMatrix;
using Rcpp::IntegerMatrix;
using Eigen::ArrayXXd;
using Eigen::ArrayXXi;
using Eigen::Ref;

typedef Eigen::Map<ArrayXXd> MapArrayXXd;

// [[Rcpp::export]]
IntegerMatrix GabrielGraph(const NumericMatrix X){

  // This function computes the Gabriel Graph of an input data.
  //
  // It first calculates the squared distance array(SDA) containing the squared euclidean distance
  // between the data points.
  //
  // Then it returns the graph, or array of adjacency.

  MapArrayXXd data(Rcpp::as<MapArrayXXd>(X));
  int n = data.rows();
  ArrayXXd fourth_power_distance_array(n, n);

  for(int i=0; i<n; i++){
    fourth_power_distance_array.col(i) = ((data.rowwise() - data.row(i)).matrix().rowwise().squaredNorm()).array().square(); // vectorized fourth power distance of col(i)
    fourth_power_distance_array(i, i) = numeric_limits<double>::infinity(); // distance between same point its infinity ( convention used in this code )
  }

  ArrayXXi array_of_adjacency = ArrayXXi::Zero(n, n);
  double min_sum_of_distances;

  for(int i=0; i<(n-1); i++){ // No need to iterate the last row.
    for(int j=(i+1); j<n; j++){ // No need to iterate over j <= i
      min_sum_of_distances = ( fourth_power_distance_array.row(i) + fourth_power_distance_array.row(j) ).minCoeff(); // get the sum of the minimum distances of the points to i and j.
      if( fourth_power_distance_array(i,j) <= min_sum_of_distances ){
        // if the sum of the minimum distances between other points to i and j isn't
        // less or equal to the distance between i and j, then (i,j) is an edge that belongs to the graph.
        array_of_adjacency(i,j) = 1;
        array_of_adjacency(j,i) = 1;
      }
    }
  }

  return wrap(array_of_adjacency);
}
