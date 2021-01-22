#' RScelestial: An R wrapper for scelestial algorithm for single-cell lineage tree reconstruction
#'   through an approximation algorithm based on Steiner tree problem
#'
#' This package provides a wrapper for the scelestial which is implemented in C++.
#' The package contains function \code{scelestial} for running the algorithm and 
#' \code{synthesis} for tumor simulation for providing synthetic data.
#' 
#' @name RScelestial
#' 
#' @useDynLib RScelestial
#' @importFrom Rcpp evalCpp
#' @importFrom "utils" "read.table"
#' @docType package
#' @exportPattern "^[^._][[:alpha:]]*"
NULL