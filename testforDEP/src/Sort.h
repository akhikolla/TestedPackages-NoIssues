#pragma once
#include <vector>
#include <Rcpp.h>
#include <R.h>

using namespace Rcpp;

class MergeSort {
private:
  std::vector<double>* tempA;
  std::vector<double>* tempB;
	std::vector<double>::iterator sa;
	std::vector<double>::iterator sb;
	void splitMerge(NumericVector::iterator a, NumericVector::iterator b, int from, int to);
	void merge(NumericVector::iterator a, NumericVector::iterator b, int from, int mid, int to);
public:
	MergeSort(NumericVector::iterator a, NumericVector::iterator b, int length);
	~MergeSort();
};

NumericMatrix SortByX(NumericVector& x, NumericVector& y);
