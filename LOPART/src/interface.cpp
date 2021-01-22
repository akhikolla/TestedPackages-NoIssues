#include "LOPART.h"
#include <Rcpp.h>

//' Low-level interface to LOPART C code
//'
//' Avoid using this function and instead use the LOPART function.
//' @title Labeled Optimal Partitioning interface
//' @param input_data numeric vector of N data to segment
//' @param input_label_start integer vector of label start positions
//'   in 0, ..., N-2
//' @param input_label_end integer vector of label end positions in 1,
//'   ..., N-1
//' @param input_label_changes integer vector of 0/1, number of
//'   labeled changes
//' @param n_updates number of dynamic programming updates to perform,
//'   usually should be number of input_data N, but can be less if you
//'   want to analyze/plot the cost/candidates at previous data.
//' @param penalty_unlabeled non-negative numeric scalar (bigger for
//'   fewer changes in unlabeled regions, smaller for more changes)
//' @param penalty_labeled non-negative numeric scalar (penalty for
//'   each change in a positive label).
//' @return data frame with four columns: cost_candidates is the cost
//'   of each last segment start considered (from 1 to N) for the
//'   computation of the optimal cost up to the last data point (Inf
//'   means infeasible); cost_optimal is the optimal cost vector
//'   computed using dynamic programming; mean is the last segment
//'   mean of the optimal model ending at that data point; last_change
//'   is the optimal changepoints (negative numbers are not used).
//' @author Toby Dylan Hocking
// [[Rcpp::export]]
Rcpp::DataFrame LOPART_interface
(Rcpp::NumericVector input_data,
 Rcpp::IntegerVector input_label_start,
 Rcpp::IntegerVector input_label_end,
 Rcpp::IntegerVector input_label_changes,
 int n_updates,
 double penalty_unlabeled,
 double penalty_labeled = 0
 ) {
  int n_data = input_data.size();
  int n_labels = input_label_changes.size();
  if(input_label_start.size() != n_labels){
    Rcpp::stop("input_label_start and input_label_changes sizes must match");
  }
  if(input_label_end.size() != n_labels){
    Rcpp::stop("input_label_end and input_label_changes sizes must match");
  }
  Rcpp::NumericVector out_cumsum(n_updates);
  Rcpp::IntegerVector out_change_candidates(n_updates);
  Rcpp::NumericVector out_cost_candidates(n_updates);
  Rcpp::NumericVector out_cost(n_updates);
  Rcpp::NumericVector out_mean(n_updates);
  Rcpp::IntegerVector out_last_change(n_updates);
  int status = LOPART
    (&input_data[0],
     n_data,
     &input_label_start[0],
     &input_label_end[0],
     &input_label_changes[0],
     n_labels,
     penalty_unlabeled,
     penalty_labeled,
     n_updates,
     //inputs above, outputs below.
     &out_cumsum[0],
     &out_change_candidates[0],
     &out_cost_candidates[0],
     &out_cost[0],
     &out_mean[0],
     &out_last_change[0]);
  if(status == ERROR_PENALTY_MUST_BE_NON_NEGATIVE){
    Rcpp::stop("penalty must be non-negative"); 
  }
  if(status == ERROR_EACH_LABEL_START_MUST_BE_LESS_THAN_ITS_END){
    Rcpp::stop("each label start must be less than its end");
  }
  if(status == ERROR_LABELED_NUMBER_OF_CHANGES_MUST_BE_0_OR_1){
    Rcpp::stop("labeled number of changes must be 0 or 1");
  }
  if(status == ERROR_EACH_LABEL_START_MUST_BE_ON_OR_AFTER_PREVIOUS_END){
    Rcpp::stop("each label start must be on or after previous end");
  }
  if(status == ERROR_LABEL_START_MUST_BE_ZERO_OR_LARGER){
    Rcpp::stop("label start must be zero or larger");
  }
  if(status == ERROR_LABEL_END_MUST_BE_LESS_THAN_N_DATA){
    Rcpp::stop("label end must be less than n data");
  }
  if(status == ERROR_NO_DATA){
    Rcpp::stop("no data");
  }
  if(status == ERROR_DATA_MUST_BE_FINITE){
    Rcpp::stop("data must be finite");
  }
  return Rcpp::DataFrame::create
    (
     Rcpp::Named("cost_candidates", out_cost_candidates),
     Rcpp::Named("cost_optimal", out_cost),
     Rcpp::Named("mean", out_mean),
     Rcpp::Named("last_change", out_last_change)
     ) ;
}

