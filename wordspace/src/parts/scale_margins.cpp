/*
 *  Scale rows and columns of a dense or sparse matrix
 */

// [[Rcpp::export]]
NumericMatrix CPP_scale_margins_dense(NumericMatrix M, NumericVector rows, NumericVector cols, bool duplicate = true) {
  int nr = M.nrow(), nc = M.ncol();
  if (nr != rows.size() || nc != cols.size())
    stop("internal error -- row/column weights not conformable with matrix");
  
  NumericMatrix res = M;
  if (duplicate)
    res = clone(M);
  
  for (int j = 0; j < nc; j++) {
    double col_weight = cols[j];
    NumericMatrix::Column v = res(_, j);
    v = v * rows * col_weight;
  }
  
  return res;
}

// [[Rcpp::export]]
S4 CPP_scale_margins_sparse(S4 M, NumericVector rows, NumericVector cols, bool duplicate = true) {
  if (!M.is("dgCMatrix"))
    stop("internal error -- not a canonical sparse matrix");
  IntegerVector dims = M.slot("Dim");
  int nr = dims[0], nc = dims[1];
  if (nr != rows.size() || nc != cols.size())
    stop("internal error -- row/column weights not conformable with matrix");

  if (duplicate)
    M = clone(M);

  IntegerVector p = M.slot("p");
  IntegerVector::iterator _p = p.begin();
  IntegerVector row_of = M.slot("i");
  IntegerVector::iterator _row_of = row_of.begin();
  NumericVector x = M.slot("x");
  NumericVector::iterator _x = x.begin();
  NumericVector::iterator _rows = rows.begin();

  for (int col = 0; col < nc; col++) {
    double col_weight = cols[col];
    for (int i = _p[col]; i < _p[col+1]; i++) {
      _x[i] *= _rows[_row_of[i]] * col_weight;
    }
  }

  return M;
}
