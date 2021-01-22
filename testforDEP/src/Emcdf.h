#pragma once
#include <vector>
#include <Rcpp.h>
#include <R.h>

using namespace Rcpp;

class Emcdf {
public:
	Emcdf(NumericVector& x, NumericVector& y, bool is_tie);
	~Emcdf();
	int n;
	bool tie;
	double** table_;
	NumericMatrix* out;
  std::vector<int>* uniqueX;
  std::vector<int>* uniqueY;
  NumericMatrix& getTable();
};

