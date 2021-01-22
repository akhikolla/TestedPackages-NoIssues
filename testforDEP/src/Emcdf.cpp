#include <Rcpp.h>
#include <R.h>
#include "Emcdf.h"
#include "Sort.h"
#include <algorithm>
#include <vector>
using namespace Rcpp;

Emcdf::Emcdf(NumericVector& x, NumericVector& y, bool is_tie)
{
  tie = is_tie;
  n = x.length();
  NumericVector sortY(n);
  for(int i=0; i<n; ++i)
    sortY[i] = y[i];

  std::sort(sortY.begin(), sortY.end());

  table_ = new double*[n];
  for(int i=0; i<n; ++i)
    table_[i] = new double[n];
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      table_[i][j] = 0;

	MergeSort mg2(x.begin(), y.begin(), n);

//compute empirical cdf
	double curY;
	for (int i = 0; i < n; ++i) {
		curY = sortY[i];
		if (y[0] <= curY)
			table_[i][0] = 1;

		for (int j = 1; j < n; ++j) {
			if (y[j] <= curY)
				table_[i][j]= table_[i][j-1] + 1;
			else
				table_[i][j] = table_[i][j-1];
		}
	}

	uniqueX = 0;
	uniqueY = 0;

	if(is_tie){
  	//get unique x, y's (ignore repeated values)
    uniqueX = new std::vector<int>;
    uniqueY = new std::vector<int>;

    double pre = x[n-1];
    uniqueX->push_back(n-1);
    for(int i=n-2; i >= 0; --i){
      if(pre != x[i]){
        pre = x[i];
        uniqueX->push_back(i);
      }
    }

    pre = sortY[n-1];
    uniqueY->push_back(n-1);
    for(int i=n-2; i >= 0; --i){
      if(pre != sortY[i]){
        pre = sortY[i];
        uniqueY->push_back(i);
      }
    }
	}


}

Emcdf::~Emcdf()
{
  delete uniqueX;
  delete uniqueY;
  for(int i=0; i<n; ++i)
  delete[] table_[i];
  delete[] table_;
  delete out;
}

NumericMatrix& Emcdf::getTable() {

  if(tie){
    int nrow = uniqueY->size();
    int ncol = uniqueX->size();
    int a = 0;
    int b = 0;

  	out = new NumericMatrix(nrow, ncol);

  	for (int i = 0; i < nrow; ++i)
  		for (int j = 0; j < ncol; ++j){
  		  a = uniqueX->at(ncol-1-j);
  		  b = uniqueY->at(nrow-1-i);
  		  out->at(i, j) = table_[b][a];
  		}


  	return *out;
  }
  else{
    out = new NumericMatrix(n, n);

  	for (int i = 0; i < n; ++i)
  		for (int j = 0; j < n; ++j)
  		  out->at(i, j) = table_[i][j];

  	return *out;
  }


}

