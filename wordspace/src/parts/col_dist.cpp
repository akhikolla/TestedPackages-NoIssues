/*
 *  Compute distances between columns of two dense or sparse matrices
 */

/* Internal codes for metric / distance:
 *  0 = euclidean
 *  1 = maximum
 *  2 = manhattan
 *  3 = minkowski   (*param1 = exponent p)
 *  4 = canberra
 *  5 = jaccard
 *  6 = minimum     (asymmetric jaccard without normalization)
 *
 * Reference implementations are in CPP_col_dist_dense(), with comments if necessary.
 * Sparse matrix code in CPP_col_dist_sparse() is often more complex and should be
 * read after consulting the reference implementation.
 */

/* make symmetric matrix from right upper triangle */
void mk_symmetric_matrix(NumericMatrix x) {
  //  int nr = x.nrow();
  int nc = x.ncol();
  for (int c = 0; c < nc; c++)
    for (int r = 0; r < c; r++)
      x(c, r) = x(r, c);
  //      x[nr * r + c] = x[nr * c + r]; /* x[c, r] = x[r, c] */
}

/* validate metric code and parameter(s) */
void check_metric(int metric_code, double p1) {
  if (metric_code < 0 || metric_code > 6)
    stop("internal error -- invalid metric code");
  if (metric_code == 3 && (!R_FINITE(p1) || p1 < 0.0))
    stop("internal error -- Minkowski p-parameter out of range [0, Inf)");  
}

// [[Rcpp::export]]
NumericMatrix CPP_col_dist_dense(NumericMatrix x, NumericMatrix y, int metric_code, double param1, bool symmetric) {
  check_metric(metric_code, param1);
  int nr = x.nrow(), nc1 = x.ncol(), nc2 = y.ncol();
  if (nr != y.nrow()) stop("internal error -- matrices are not conformable");

  NumericMatrix dist(nc1, nc2);

#pragma omp parallel for \
        if (openmp_threads > 1 && (nc1 + 0.0) * (nc2 + 0.0) * (nr + 0.0) > 100e6) \
        num_threads(openmp_threads) \
        shared(dist)
  for (int col2 = 0; col2 < nc2; col2++) {
    NumericVector tmp(nr);
    double accum, denom;
    int col1_max = (symmetric) ? col2 + 1 : nc1;
    for (int col1 = 0; col1 < col1_max; col1++) {
      NumericMatrix::Column vx = x(_, col1);  // column <col1> of first matrix
      NumericMatrix::Column vy = y(_, col2);  // column <col2> of second matrix
      switch (metric_code) {
      case 0: 
        accum = sum((vx - vy) * (vx - vy));
        dist(col1, col2) = sqrt(accum);
        break;
      case 1:
        dist(col1, col2) = max(abs(vx - vy));
        break;
      case 2:
        dist(col1, col2) = sum(abs(vx - vy));
        break;
      case 3:
        accum = sum(pow(abs(vx - vy), param1));
        if (param1 > 1.0)
          dist(col1, col2) = pow(accum, 1.0 / param1);
        else
          dist(col1, col2) = accum;
        break;
      case 4:
        tmp = abs(vx) + abs(vy); // denominator |x_i| + |y_i|
        dist(col1, col2) = sum(ifelse(tmp > 0, abs(vx - vy) / tmp, 0.0));
        break;
      case 5:
        accum = sum(pmin(vx, vy)); /* x_i and y_i must be non-negative */
        denom = sum(pmax(vx, vy));
        dist(col1, col2) = (denom > 0) ? accum / denom : 1.0; /* special case: J(0, 0) = 1 */
        break;
      case 6:
        dist(col1, col2) = sum(pmin(vx, vy)); /* overlap of probability distributions (if ||x||_1 = ||y||_1 = 1) */
        break;
      }
    }
  }
  
  if (symmetric) mk_symmetric_matrix(dist);
  return dist;
}

// [[Rcpp::export]]
NumericMatrix CPP_col_dist_sparse(int nc1, IntegerVector xp, IntegerVector xrow, NumericVector x, int nc2, IntegerVector yp, IntegerVector yrow, NumericVector y, int metric_code, double param1, bool symmetric) {
  check_metric(metric_code, param1);
  NumericVector::iterator _x = x.begin();
  NumericVector::iterator _y = y.begin();
  IntegerVector::iterator _xrow = xrow.begin();
  IntegerVector::iterator _yrow = yrow.begin();

  NumericMatrix dist(nc1, nc2);

#ifdef _OPENMP
  /* average number of entries scanned when comparing two columns, used to decide whether to try parallelization */
  double avg_nr = (xp[nc1] - xp[0] + 0.0) / nc1 + (yp[nc2] - yp[0] + 0.0) / nc2; 
#endif

#pragma omp parallel for \
        if (openmp_threads > 1 && (nc1 + 0.0) * (nc2 + 0.0) * avg_nr > 40e6) \
        num_threads(openmp_threads) \
        shared(dist, nc1, xp, _xrow, _x, nc2, yp, _yrow, _y, metric_code, param1)
  for (int col2 = 0; col2 < nc2; col2++) {
    int col1_max = (symmetric) ? col2 + 1 : nc1;
    int yi_max = yp[col2 + 1];

    for (int col1 = 0; col1 < col1_max; col1++) {
      int xi_max = xp[col1 + 1];
      int xi = xp[col1];
      int yi = yp[col2];
      int xrow_curr = (xi < xi_max) ? _xrow[xi] : INT_MAX;
      int yrow_curr = (yi < yi_max) ? _yrow[yi] : INT_MAX;
      
      double accum = 0.0, denom = 0.0;
      double x_curr, y_curr;
      double d_xy, x_plus_y;
      while (xi < xi_max || yi < yi_max) {

        if (xrow_curr < yrow_curr) {
          x_curr = _x[xi]; y_curr = 0.0;
          xi++;
          xrow_curr = (xi < xi_max) ? _xrow[xi] : INT_MAX;
        }
        else if (xrow_curr == yrow_curr) {
          x_curr = _x[xi]; y_curr = _y[yi];
          xi++; yi++;
          xrow_curr = (xi < xi_max) ? _xrow[xi] : INT_MAX;
          yrow_curr = (yi < yi_max) ? _yrow[yi] : INT_MAX;
        }
        else /* xrow_curr > yrow_curr */ {
          x_curr = 0; y_curr = _y[yi];
          yi++;
          yrow_curr = (yi < yi_max) ? _yrow[yi] : INT_MAX;          
        }
        
        switch (metric_code) {
        case 0:
          d_xy = x_curr - y_curr;
          accum += d_xy * d_xy;
          break;
        case 1:
          d_xy = fabs(x_curr - y_curr);
          if (d_xy > accum) accum = d_xy;
          break;
        case 2:
          d_xy = fabs(x_curr - y_curr);
          accum += d_xy;
          break;
        case 3:
          d_xy = fabs(x_curr - y_curr);
          accum += pow(d_xy, param1);
          break;
        case 4:
          x_plus_y = fabs(x_curr) + fabs(y_curr);
          d_xy = fabs(x_curr - y_curr);
          if (x_plus_y > 0) accum += d_xy / x_plus_y;
          break;
        case 5:
          if (x_curr >= y_curr) {
            accum += y_curr;  /* min(x_i, y_i) */
            denom += x_curr; /* max(x_i, y_i) */
          }
          else {
            accum += x_curr;
            denom += y_curr;
          }
          break;
        case 6:
          accum += (x_curr >= y_curr) ? y_curr : x_curr; /* min(x_i, y_i) */
        }
      } /* while (xi, yi) */

      switch (metric_code) {
      case 0:
        dist(col1, col2) = sqrt(accum);
        break;
      case 3:
        if (param1 > 1.0)
          dist(col1, col2) = pow(accum, 1.0 / param1);
        else
          dist(col1, col2) = accum;
        break;
      case 5:
        dist(col1, col2) = (denom > 0) ? accum / denom : 1.0; /* (SUM min(x_i, y_i)) / (SUM max(x_i, y_i)) */
        break;
      case 1:
      case 2:
      case 4:
      case 6:
        dist(col1, col2) = accum;
        break;
      }

    } /* for (col1) */
  } /* for (col2) */
  
  if (symmetric) mk_symmetric_matrix(dist);
  return dist;
}
