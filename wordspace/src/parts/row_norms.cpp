/*
 *  Compute different row norms for dense or sparse matrix
 */

/* internal codes for norms:
 *  0 = euclidean
 *  1 = maximum
 *  2 = manhattan
 *  3 = minkowski  (with exponent *p_norm)
 */

void check_norm(int norm_code, double p) {
  if (norm_code < 0 || norm_code > 3)
    stop("internal error -- invalid norm code");
  if (norm_code == 3 && (!R_FINITE(p) || p < 0.0))
    stop("internal error -- Minkowski p-parameter out of range [0, Inf)");  
}

// [[Rcpp::export]]
NumericVector CPP_row_norms_dense(NumericMatrix x, int norm_code, double p_norm = 2.0) {
  check_norm(norm_code, p_norm);

  int nr = x.nrow(), nc = x.ncol();
  NumericVector norms(nr, 0.0);

  NumericMatrix::iterator _x = x.begin();
  NumericVector::iterator _norms = norms.begin();
  
  int i = 0;
  for (int col = 0; col < nc; col++) {
    for (int row = 0; row < nr; row++) {
      if      (norm_code == 0) _norms[row] += _x[i] * _x[i];
      else if (norm_code == 1) { if (fabs(_x[i]) > _norms[row]) _norms[row] = fabs(_x[i]); }
      else if (norm_code == 2) _norms[row] += fabs(_x[i]);
      else if (norm_code == 3) {
        if (p_norm > 0)
          _norms[row] += pow(fabs(_x[i]), p_norm);
        else
          _norms[row] += (_x[i] != 0);
      }
      i++;
    }
  }
  
  if      (norm_code == 0)                 norms = sqrt(norms);
  else if (norm_code == 3 && p_norm > 1.0) norms = pow(norms, 1.0 / p_norm);
  /* no adjustment needed for Maximum and Manhattan norms */
  
  return norms;
}

// [[Rcpp::export]]
NumericVector CPP_row_norms_sparse(int nr, int nc, IntegerVector p, IntegerVector row_of, NumericVector x, int norm_code, double p_norm = 2.0) {
  check_norm(norm_code, p_norm);

  NumericVector norms(nr, 0.0);

  NumericVector::iterator _x = x.begin();
  IntegerVector::iterator _p = p.begin();
  IntegerVector::iterator _row_of = row_of.begin();
  NumericVector::iterator _norms = norms.begin();

  for (int col = 0; col < nc; col++) {
    for (int i = _p[col]; i < _p[col+1]; i++) {
      int row = _row_of[i];
      if      (norm_code == 0) _norms[row] += _x[i] * _x[i];
      else if (norm_code == 1) { if (fabs(_x[i]) > _norms[row]) _norms[row] = fabs(_x[i]); }
      else if (norm_code == 2) _norms[row] += fabs(_x[i]);
      else if (norm_code == 3) {
        if (p_norm > 0)
          _norms[row] += pow(fabs(_x[i]), p_norm);
        else
          _norms[row] += (_x[i] != 0);
      }
    }
  }

  if      (norm_code == 0)                 norms = sqrt(norms);
  else if (norm_code == 3 && p_norm > 1.0) norms = pow(norms, 1.0 / p_norm);
  /* no adjustment needed for Maximum and Manhattan norms */
  
  return norms;
}

// [[Rcpp::export]]
NumericVector CPP_col_norms_dense(NumericMatrix x, int norm_code, double p_norm = 2.0) {
  check_norm(norm_code, p_norm);

  // int nr = x.nrow();
  int nc = x.ncol();
  NumericVector norms(nc, 0.0);

  // NumericVector::iterator _norms = norms.begin();

  for (int col = 0; col < nc; col++) {
    NumericMatrix::Column v = x(_, col);
    if      (norm_code == 0) norms[col] = sum(v * v);
    else if (norm_code == 1) norms[col] = max(abs(v));
    else if (norm_code == 2) norms[col] = sum(abs(v));
    else if (norm_code == 3) {
      if (p_norm > 0) 
        norms[col] = sum(pow(abs(v), p_norm));
      else
        norms[col] = sum(v != 0);
    }
  }

  if      (norm_code == 0)                 norms = sqrt(norms);
  else if (norm_code == 3 && p_norm > 1.0) norms = pow(norms, 1.0 / p_norm);
  /* no adjustment needed for Maximum and Manhattan norms */
  
  return norms;
}

// [[Rcpp::export]]
NumericVector CPP_col_norms_sparse(int nr, int nc, IntegerVector p, IntegerVector row_of, NumericVector x, int norm_code, double p_norm = 2.0) {
  check_norm(norm_code, p_norm);

  NumericVector norms(nc, 0.0);

  NumericVector::iterator _x = x.begin();
  IntegerVector::iterator _p = p.begin();
  // IntegerVector::iterator _row_of = row_of.begin();
  NumericVector::iterator _norms = norms.begin();

  for (int col = 0; col < nc; col++) {
    double accum = 0.0;
    for (int i = _p[col]; i < _p[col+1]; i++) {
      if      (norm_code == 0) accum += _x[i] * _x[i];
      else if (norm_code == 1) { if (fabs(_x[i]) > accum) accum = fabs(_x[i]); }
      else if (norm_code == 2) accum += fabs(_x[i]);
      else if (norm_code == 3) {
        if (p_norm > 0)
          accum += pow(fabs(_x[i]), p_norm);
        else
          accum += (_x[i] != 0);
      }
    }
    if      (norm_code == 0)                 _norms[col] = sqrt(accum);
    else if (norm_code == 3 && p_norm > 1.0) _norms[col] = pow(accum, 1.0 / p_norm);
    else    /* other norms */                _norms[col] = accum;
  }

  return norms;
}
