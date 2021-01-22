// Various smaller utility functions implemented in C++ for efficiency and reduced memory overhead.


// Efficient nonzero count and non-negativity check for dense and sparse matrices,
// based on counting the signs (pos, 0, neg) of values in a numeric vector.
// [[Rcpp::export]]
NumericVector CPP_signcount(NumericVector x) {
  int pos = 0, neg = 0, zero = 0;
  int n = x.size();
  NumericVector::iterator _x = x.begin();
  for (int i = 0; i < n; i++) {
    if (_x[i] > 0.0)
      pos++;
    else if (_x[i] < 0.0)
      neg++;
    else
      zero++;
  }
  return NumericVector::create(pos, zero, neg);
}

// Same for integer vector to avoid memory overhead from automatic type conversion.
// [[Rcpp::export]]
NumericVector CPP_signcount_int(IntegerVector x) {
  int pos = 0, neg = 0, zero = 0;
  int n = x.size();
  IntegerVector::iterator _x = x.begin();
  for (int i = 0; i < n; i++) {
    if (_x[i] > 0.0)
      pos++;
    else if (_x[i] < 0.0)
      neg++;
    else
      zero++;
  }
  return NumericVector::create(pos, zero, neg);
}


// TODO: implement integer version (to avoid conversion)
// TODO: R wrapper function
//  - is.integer
//  - is.double
//  - Matrix classes: dgeMatrix (dense), dgCMatrix, dgRMatrix, dgTMatrix (with adjusted zero count) -- all should have slot @x
//  - Matrix doesn't seem to really support integer matrices yet, so we don't worry abpout those
