

#include <Rcpp.h>
#include <R.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include "Sort.h"

using namespace Rcpp;

MergeSort::MergeSort(NumericVector::iterator a, NumericVector::iterator b, int length) {

  tempA = new std::vector<double>(length);
  tempB = new std::vector<double>(length);
	sa = tempA->begin();
	sb = tempB->begin();
	splitMerge(a, b, 0, length);
}

MergeSort::~MergeSort(){
  delete tempA;
  delete tempB;
}

void MergeSort::merge(NumericVector::iterator a, NumericVector::iterator b, int from, int mid, int to) {

	int l = from;
	int r = mid;
	int k = from;

	while (l<mid && r<to) {
		if (a[l] < a[r]) {
			sa[k] = a[l];
			sb[k++] = b[l++];
		}
		else {
			sa[k] = a[r];
			sb[k++] = b[r++];
		}
	}
	while (l < mid) {
		sa[k] = a[l];
		sb[k++] = b[l++];
	}
	while (r < to) {
		sa[k] = a[r];
		sb[k++] = b[r++];
	}
	for (int i = from; i < to; i++) {
		a[i] = sa[i];
	}
	for (int i = from; i < to; i++) {
		b[i] = sb[i];
	}


}

void MergeSort::splitMerge(NumericVector::iterator a, NumericVector::iterator b, int from, int to) {


	if ((to - from) >1) {
		int mid = (from + to) / 2;
		splitMerge(a, b, from, mid);
		splitMerge(a, b, mid, to);
		merge(a, b, from, mid, to);
	}
}

// [[Rcpp::export]]
NumericMatrix SortByX(NumericVector& x, NumericVector& y){
  NumericVector::iterator a = x.begin();
  NumericVector::iterator b = y.begin();
  int length = x.length();
  MergeSort(a,b,length);
  return Rcpp::cbind(x,y);
}

