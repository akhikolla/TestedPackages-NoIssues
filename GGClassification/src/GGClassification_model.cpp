// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <vector>
#include <limits>

using std::numeric_limits;
using std::vector;
using Rcpp::as;
using Rcpp::wrap;
using Rcpp::stop;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::NumericMatrix;
using Rcpp::IntegerVector;
using Eigen::ArrayXXd;
using Eigen::ArrayXXi;
using Eigen::ArrayXd;
using Eigen::ArrayXi;
using Eigen::VectorXi;
using Eigen::Ref;

typedef Eigen::Map<ArrayXXd> MapArrayXXd;

void RemoveArrayElementsByIndex(const ArrayXXd& data, int nrows, Ref<ArrayXXd> updated_data, vector<int> index_of_element_to_remove){

  // Receives an array of data and an vector containing the indexes of rows to remove.
  // An empty array is parsed as reference to receive the data without the rows removed.

  int iter = 0;
  for(int i=0; i<nrows; i++){
    if( i == index_of_element_to_remove[0] && int(index_of_element_to_remove.size()) != 0 ){
      index_of_element_to_remove.erase(index_of_element_to_remove.begin());
      continue;
    }
    updated_data.row(iter) = data.row(i);
    iter++;
  }

}

void RemoveVectorElementsByIndex(vector<int>& vector_of_classes, vector<int> index_of_element_to_remove){

  // Receives an vector of labels and an vector containing the indexes of elements to remove.
  // The vector of labels is parsed as reference and modified inside de function.

  index_of_element_to_remove.push_back(-1);
  index_of_element_to_remove.push_back(vector_of_classes.size());
  sort(index_of_element_to_remove.begin(), index_of_element_to_remove.end());
  vector<int>::iterator last = vector_of_classes.begin();
  for(size_t i = 1; i != index_of_element_to_remove.size(); ++i){
    size_t range_begin = index_of_element_to_remove[i - 1] + 1;
    size_t range_end = index_of_element_to_remove[i];
    copy(vector_of_classes.begin() + range_begin, vector_of_classes.begin() + range_end, last);
    last += range_end - range_begin;
  }
  vector_of_classes.erase(last, vector_of_classes.end());

}

void GabrielGraph(const ArrayXXd& data, int nrows, Ref<ArrayXXi> array_of_adjacency){

  // This function computes the Gabriel Graph of an input data.
  //
  // It first calculates the squared distance array(SDA) containing the squared euclidean distance
  // between the data points.
  //
  // Then, using the SDA, computes the gabriel graph and stores in an array parsed by reference.

  ArrayXXd fourth_power_distance_array(nrows, nrows);

  for(int i=0; i<nrows; i++){
    fourth_power_distance_array.col(i) = ((data.rowwise() - data.row(i)).matrix().rowwise().squaredNorm()).array().square(); // vectorized fourth power distance of col(i)
    fourth_power_distance_array(i, i) = numeric_limits<double>::infinity(); // distance between same point its infinity ( convention used in this code )
  }

  double min_sum_of_distances;

  for(int i=0; i<(nrows-1); i++){ // No need to iterate the last row.
    for(int j=(i+1); j<nrows; j++){ // No need to iterate over j <= i
      min_sum_of_distances = ( fourth_power_distance_array.row(i) + fourth_power_distance_array.row(j) ).minCoeff(); // get the sum of the minimum distances of the points to i and j.
      // if the sum of the minimum distances between other points to i and j isn't
      // less or equal to the distance between i and j, then (i,j) is an edge that belongs to the graph.
      if( fourth_power_distance_array(i,j) <= min_sum_of_distances ){
        array_of_adjacency(i,j) = 1;
        array_of_adjacency(j,i) = 1;
      }
    }
  }

}

vector<int> FilterGraph(const ArrayXXi& array_of_adjacency, const vector<int>& vector_of_classes, const vector<int> labels){

  // This function receives the graph and its labels and filter the noise using an
  // vertex quality rule. It returns an vector containing the indexes of the noisy data points.
  //
  // The quality of the vertex is the number of connections(degree) between same label vertexes
  // and number of all connections of the vertex.
  //
  // The threshold that say if the vertex is an noise is the sum of all qualities of the vertexes
  // belonging to the same class, divided by the length of the current class.
  // ( this will give us two thresholds )

  vector<int> index_of_element_to_remove;

  int nrows = array_of_adjacency.rows();
  int size_of_class_1 = 0, size_of_class_2 = 0;
  double threshold_class_1 = 0, threshold_class_2 = 0;
  double degree_of_vertex, degree_of_vertex_same_class, quality_of_vertex;

  ArrayXd array_of_vertex_quality(nrows);
  VectorXi mask_of_class_1 = VectorXi::Zero(nrows);
  VectorXi mask_of_class_2 = VectorXi::Zero(nrows);

  for(int i=0; i<nrows; i++){
    if( vector_of_classes[i] == labels[0] ){
      size_of_class_1++;
      mask_of_class_1(i) = 1;
    }
    else{
      size_of_class_2++;
      mask_of_class_2(i) = 1;
    }
  }

  bool flag_noise = false; // Flag that indicates if the dataset has noise.

  for(int i=0; i<nrows; i++){ // Getting the quality of the vertexes and the thresholds of each class.

    degree_of_vertex = array_of_adjacency.row(i).sum();
    if(vector_of_classes[i] == labels[0]){
      degree_of_vertex_same_class = array_of_adjacency.row(i).matrix().dot(mask_of_class_1);
      quality_of_vertex = (degree_of_vertex_same_class/degree_of_vertex);
      threshold_class_1 += quality_of_vertex;
    }
    else{
      degree_of_vertex_same_class = array_of_adjacency.row(i).matrix().dot(mask_of_class_2);
      quality_of_vertex = (degree_of_vertex_same_class/degree_of_vertex);
      threshold_class_2 += quality_of_vertex;
    }

    if( quality_of_vertex == 0 ){ // Verifying if the current vertex is a noise.
      flag_noise = true;
    }

    array_of_vertex_quality(i) = quality_of_vertex;

  }

  if(flag_noise){ // If dataset contains noise, then apply filtering process.
    threshold_class_1 /= size_of_class_1; // Dividing by class length
    threshold_class_2 /= size_of_class_2;

    for(int i=0; i<nrows; i++){ // Filtering graph based on thresolds and qualities.
      if( vector_of_classes[i] == labels[0] ){
        if( array_of_vertex_quality(i) < threshold_class_1 ){
          index_of_element_to_remove.push_back(i);
        }
      }
      else{
        if( array_of_vertex_quality(i) < threshold_class_2 ){
          index_of_element_to_remove.push_back(i);
        }
      }
    }

    return(index_of_element_to_remove); // vector containing the indexes of elements to remove.
  }
  else{
    return(index_of_element_to_remove); // Empty vector, zero elements to remove.
  }

}

List GetModelParams(const ArrayXXi& array_of_adjacency, const ArrayXXd& data, int nrows, int ncols, const vector<int>& vector_of_classes, const vector<int> labels){

  // Extracts the parameters needed to classification of new data points, that will be
  // done by the predict function.
  //
  // Receives the already filtered graph, the data points and the labels. It will extract
  // the paremeters midpoints, w and bias.
  //
  //The return is an Named List that will be the main return of the program.


  int current_vertex_class_1, current_vertex_class_2;

  vector<int> indexes_of_class_1, indexes_of_class_2;

  for(int i=0; i<nrows; i++){
    if( vector_of_classes[i] == labels[0] ){
      indexes_of_class_1.push_back(i);
    }
    else{
      indexes_of_class_2.push_back(i);
    }
  }

  int size_of_class_1 = indexes_of_class_1.size();
  int size_of_class_2 = indexes_of_class_2.size();

  const int size_to_reserve = size_of_class_1*size_of_class_2;

  vector<ArrayXd> vector_midpoints;
  vector<ArrayXd> vector_w;
  vector<double> vector_bias;

  vector_midpoints.reserve(size_to_reserve);
  vector_w.reserve(size_to_reserve);
  vector_bias.reserve(size_to_reserve);

  ArrayXXd vertex_data(2, ncols);
  ArrayXd midpoint(ncols);
  ArrayXd w(ncols);
  int sig;
  double bias;

  for(int i=0; i<size_of_class_1; i++){
    current_vertex_class_1 = indexes_of_class_1[i];
    for(int j=0; j<size_of_class_2; j++){
      current_vertex_class_2 = indexes_of_class_2[j];
      if(array_of_adjacency(current_vertex_class_1, current_vertex_class_2) == 1){

        vertex_data.row(0) = data.row(current_vertex_class_1);
        vertex_data.row(1) = data.row(current_vertex_class_2);

        midpoint = ( vertex_data.row(0) + vertex_data.row(1) )/2;
        vector_midpoints.push_back(midpoint);

        w = ( vertex_data.row(1) - vertex_data.row(0) );
        vector_w.push_back(w);

        bias = (midpoint*w).sum();
        vector_bias.push_back(bias);
      }
    }
  }

  // Converting vectors to arrays for better return structure.
  int midpoints_length = vector_midpoints.size();
  ArrayXXd array_midpoints(midpoints_length, ncols);
  ArrayXXd array_w(midpoints_length, ncols);
  for(int k=0; k<midpoints_length; k++){
    array_midpoints.row(k) = vector_midpoints[k];
    array_w.row(k) = vector_w[k];
  }

  return List::create(
    Named("Midpoints") = wrap(array_midpoints),
    Named("W") = wrap(array_w),
    Named("Bias") = wrap(vector_bias),
    Named("Labels") = wrap(labels) // Auxiliar variable used in the predict function.
  );
}

vector<int> VerificationOfParameters(const ArrayXXd& X, const vector<int> Y){

  // Verifies if the number of data points and the length of the labels is equivalent.
  // Also, test the labels contains two different types only and return this two types.
  // Any break of the conditions above an error is raised, returning to R enviroment.

  int nrows_X = X.rows();
  int length_Y = Y.size();

  if( nrows_X != length_Y ){
    stop("Error: The number of data elements(x) should be equal to the number of class labels(y).");
  }

  vector<int> labels = Y;
  sort(labels.begin(), labels.end());
  vector<int>::iterator last = unique(labels.begin(), labels.end());
  labels.erase(last, labels.end());

  int number_of_labels = labels.size();
  if( number_of_labels  != 2 ){
    stop("Error: The current classifier model only supports two labels of classes.");
  }

  return(labels);
}

// [[Rcpp::export]]
List model(NumericMatrix& X, IntegerVector& y, bool normalize=false){

  // Main function of the model.
  // It invokes all of the steps necessary to obtain the model, and returns it
  // as an list.

  MapArrayXXd X_array(as<MapArrayXXd>(X)); // Conversion of R data to C++/Eigen data type.
  vector<int> Y = as<vector<int> >(y);

  vector<int> labels = VerificationOfParameters(X_array, Y);

  if(normalize){ // Normalization of data if user desiries.
    X_array = ( ( X_array-X_array.mean() )/( (X_array-X_array.mean()).abs().maxCoeff() ) );
  }

  int nrows = X_array.rows(); // data dimensions
  int ncols = X_array.cols();

  ArrayXXi gabriel_graph_unfiltered = ArrayXXi::Zero(nrows,nrows); // instantiating the graph
  GabrielGraph(X_array, nrows, gabriel_graph_unfiltered); // unfiltered graph

  vector<int> index_of_element_to_remove = FilterGraph(gabriel_graph_unfiltered, Y, labels); // obtaining indexes of noises.

  List final_model;

  if(index_of_element_to_remove.empty()){ // No noise to filter.

    final_model = GetModelParams(gabriel_graph_unfiltered, X_array, nrows, ncols, Y, labels); //extracting model params

  }
  else{ // Data has noise and have to be filtered.

    int new_nrows = nrows - index_of_element_to_remove.size();
    ArrayXXd X_filtered(new_nrows, ncols);
    RemoveArrayElementsByIndex(X_array, nrows, X_filtered, index_of_element_to_remove); // filtering the data
    RemoveVectorElementsByIndex(Y, index_of_element_to_remove);

    ArrayXXi gabriel_graph_filtered = ArrayXXi::Zero(new_nrows, new_nrows);
    GabrielGraph(X_filtered, new_nrows, gabriel_graph_filtered); // obtaining filtered graph

    final_model = GetModelParams(gabriel_graph_filtered, X_filtered, new_nrows, ncols, Y, labels); //extracting model parameters

  }

  return(final_model);
}
