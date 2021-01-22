/*
 *  Amalgamation file of Rcpp/C functions for "wordspace" package
 *  (source code merged into single file in order to avoid slow re-loading of Rcpp headers)
 */

#include <Rcpp.h>
using namespace Rcpp;

/*
 *  Compute different association measures from frequency signatures
 */

#ifndef wordspace_am_h
#define wordspace_am_h

typedef double (*am_func)(double f, double f1, double f2, double N, int sparse); 

extern int am_table_entries; /* number of AM function pointers in am_func table */
extern am_func am_table[];

double transform(double x, int method);

#endif /* wordspace_am_h */
#ifndef wordspace_globals_h
#define wordspace_globals_h

extern int openmp_threads;

#endif /* wordspace_globals_h */
/*
 *  Compute different association measures from frequency signatures
 */

double am_frequency(double f, double f1, double f2, double N, int sparse) {
  return f;
}

double am_simple_ll(double f, double f1, double f2, double N, int sparse) {
  double O = f, E = f1 * f2 / N;
  double ll = 2 * ( ((O > 0) ? O * log(O / E) : 0) - (O - E) );
  if (sparse)
    return (O > E) ? ll : 0;
  else
    return (O >= E) ? ll : -ll;
}

double am_t_score(double f, double f1, double f2, double N, int sparse) {
  double O = f, E = f1 * f2 / N;
  if (sparse) 
    return (O > E) ? (O - E) / sqrt(O) : 0;
  else
    return (O - E) / sqrt(O + 1); /* "discounted" t-score for O == 0 */
}

double am_z_score(double f, double f1, double f2, double N, int sparse) {
  double O = f, E = f1 * f2 / N;
  if (sparse)
    return (O > E) ? (O - E) / sqrt(E) : 0;
  else
    return (O - E) / sqrt(E); /* E == 0 should never happen */
}

double am_Dice(double f, double f1, double f2, double N, int sparse) {
  return 2 * f / (f1 + f2);
}

double am_MI(double f, double f1, double f2, double N, int sparse) {
  double O = f, E = f1 * f2 / N;
  if (sparse)
    return (O > E) ? log2(O / E) : 0;
  else
    return log2(O / E); /* not clear how to avoid the -Inf result here */
}

double am_tf_idf(double f, double f1, double f2, double N, int sparse) {
  /* f1 = dummy, f2 = df, N = total document count (set to 1 if f2 holds relative df) */
  return (f2 > 0) ? f * log(N / f2) : 0; /* avoid division by zero if f2 == 0 */
}

double am_log_likelihood(double f, double f1, double f2, double N, int sparse) {
  double R1 = f1, R2 = N - f1, C1 = f2, C2 = N - f2;
  double O11 = f, O12 = R1 - f, O21 = C1 - f, O22 = C2 - O12;
  double E11 = R1 * C1 / N, E12 = R1 * C2 / N, E21 = R2 * C1 / N, E22 = R2 * C2 / N;
  double G2 =
    ((O11 > 0) ? O11 * log(O11 / E11) : 0) +
    ((O12 > 0) ? O12 * log(O12 / E12) : 0) +
    ((O21 > 0) ? O21 * log(O21 / E21) : 0) +
    ((O22 > 0) ? O22 * log(O22 / E22) : 0);
  if (sparse)
    return (O11 > E11) ? 2 * G2 : 0;
  else
    return (O11 >= E11) ? 2 * G2 : -2 * G2;
}

double am_chi_squared(double f, double f1, double f2, double N, int sparse) {
  double R1 = f1, R2 = N - f1, C1 = f2, C2 = N - f2;
  double O11 = f, O12 = R1 - f, O21 = C1 - f, O22 = C2 - O12;
  double E11 = R1 * C1 / N;
  double yates = fabs(O11 * O22 - O12 * O21) - N / 2;
  double X2 = N * yates * yates / (R1 * R2 * C1 * C2);
  if (sparse)
    return (O11 > E11) ? X2 : 0;
  else
    return (O11 >= E11) ? X2 : -X2;
}

double transform(double x, int method) {
  switch (method) {
    case 0:       /* 0 = none */
      return x;
    case 1:       /* 1 = signed log */
      return R::sign(x) * log(fabs(x) + 1);
    case 2:       /* 2 = signed square root */
      return R::sign(x) * sqrt(fabs(x));
    case 3:       /* 3 = sigmoid (tanh) */
      return tanh(x);
    default:
      stop("internal error -- invalid score transformation code");
      return 0.0; /* just to keep clang from bitching */
  }
}

int am_table_entries = 9;

am_func am_table[] = {
  &am_frequency,
  &am_simple_ll,
  &am_t_score,
  &am_z_score,
  &am_Dice,
  &am_MI,
  &am_tf_idf,
  &am_log_likelihood,
  &am_chi_squared
};

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
/*
 * check OpenMP availability and set desired number of threads
 */

#ifdef _OPENMP
#include <omp.h>
#endif

int openmp_threads = 1;

// [[Rcpp::export]]
DataFrame CPP_get_openmp_threads() {
  int num_threads = openmp_threads;
#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
#else
  int max_threads = 0;
#endif
  DataFrame res =
    DataFrame::create(_["available"] = max_threads > 0,
                      _["max"] = max_threads, 
                      _["threads"] = num_threads);
  res.attr("row.names") = "OpenMP";
  return res;
}

// [[Rcpp::export]]
void CPP_set_openmp_threads(int n) {
  if (n < 1) stop("internal error -- number of threads must be >= 1");
#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  if (n > max_threads) n = max_threads;
  openmp_threads = n;
#else
  if (n > 1) Rf_warning("OpenMP support not available");
  openmp_threads = 1;
#endif
}

/*
 *  Memory-friendly random indexing of a sparse matrix
 */

/* Memory-friendly implementation of random indexing uses a very simple algorithm
 * that iterates through the dimensions of the original matrix M.  For each dimensions,
 * it generates a cross-section of the random basis vectors for this dimension, i.e.
 * a column of the random basis matrix Q, and updates the projected vectors with the
 * outer product of the column of Q and the corresponding column of M.
 * Since both M and Q are sparse, the outer product requires are relatively small number
 * of updates to the projected vectors that can be generated very efficiently.
 * Note that this algorithm isn't cache-friendly at all and will be relatively slow,
 * but it needs only a minimal amount of extra working memory.
 */

// [[Rcpp::export]]
NumericMatrix CPP_random_indexing_sparse(int nr, int nc, IntegerVector p, IntegerVector row_of, NumericVector x, int n_ri, double rate, bool verbose = true) {
  int target_fill = n_ri * rate;                    /* expected number of nonzero entries in each column of Q */
  int max_fill = 2 * target_fill + 1;               /* maximal number of nonzero entries */
  if (max_fill > nc) max_fill = nc;
  std::vector<double> Q_x(max_fill);                /* nonzero entries of current column of Q */
  std::vector<int> Q_row_of(max_fill);              /* row offsets of nonzero entries in current column of Q */
  int nnzero = 0;                                   /* number of nonzero entries in current column of Q */
  std::vector<double> Q_norms(n_ri);                /* collect Euclidean norms of dynamically generated basis vectors */
  double n_updates = 0.0;                           /* performance statistics for <verbose> mode */

  NumericMatrix res(nr, n_ri);                      /* allocate matrix of projected vectors (should be 0-initialised!) */
  for (int i = 0; i < n_ri; i++) Q_norms[i] = 0.0;  /* initialise array of basis vector norms */
  if (rate < 1.0 / n_ri) rate = 1.0 / n_ri; // make sure at least one nonzero entry is expected in each column

  for (int col = 0; col < nc; col++) {
    /* generate column of Q as random vector with fill rate <rate> and values +1 / -1 with equal probability */
    nnzero = 0;
    int rd = R::rgeom(rate);                        /* generate offset of first nonzero entry from geometric distribution */
    while (rd >= n_ri) rd = R::rgeom(rate);         /* make sure there's at least one nonzero entry */
    while (rd < n_ri && nnzero < max_fill) {
      Q_row_of[nnzero] = rd;
      Q_x[nnzero] = (R::unif_rand() >= 0.5) ? 1.0 : -1.0;
      Q_norms[rd] += 1.0;
      nnzero++;
      rd += R::rgeom(rate) + 1;
    }

    /* now compute sparse outer product of matching columns of M and Q, and update <res>, i.e. res += Q[,col] %*% t(M[,col]) */
    for (int i = p[col]; i < p[col+1]; i++) {
      for (int j = 0; j < nnzero; j++) {
        /* k = M_row_of[i], m = Q_row_of[j], M[k, col] = x[i], Q[m, col] = Q_x[j] */
        res(row_of[i], Q_row_of[j]) += Q_x[j] * x[i];   /* res[k, m] += Q[m, col] * M[k, col] */
        n_updates++;
      }
    }
    
    if (verbose && ((col+1) % 100000) == 0)
      Rprintf("%6.0fk columns processed (%.1fG memory updates)\n", (col+1.0)/1000, n_updates / 1e9);
  }
  if (verbose) Rprintf("%.1fG memory updates complete, rescaling RI dimensions\n", n_updates / 1e9);
  
  /* rescale columns of <res>, i.e. random dimensions, with Euclidean norm of basis vectors */
  for (int rd = 0; rd < n_ri; rd++) {
    if (Q_norms[rd] > 0) {
      double factor = 1.0 / sqrt(Q_norms[rd]);
      NumericMatrix::Column v = res(_, rd);
      v = v * factor;
    } // else this dimension is all zeroes (both basis vector and projection)
  }
  
  return res;
}
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
/*
 *  Compute association scores for dense or sparse matrix
 */

// [[Rcpp::export]]
NumericMatrix CPP_dsm_score_dense(NumericMatrix f, NumericVector f1, NumericVector f2, double N, int am_code, int sparse, int transform_code) {
  if (am_code < 0 || am_code >= am_table_entries)
    stop("internal error -- invalid AM code");
  am_func AM = am_table[am_code]; /* selected association measure */

  int nr = f.nrow(), nc = f.ncol();
  if (am_code != 0 && (nr != f1.size() || nc != f2.size()))
    stop("internal error -- marginal vectors f1 and f2 not conformable with matrix f");

  NumericMatrix scores(nr, nc);

  NumericMatrix::iterator _f = f.begin();
  NumericVector::iterator _f1 = f1.begin();
  NumericVector::iterator _f2 = f2.begin();
  NumericMatrix::iterator _scores = scores.begin();

  int i = 0;
  for (int col = 0; col < nc; col++) {
    for (int row = 0; row < nr; row++) {
      /* frequeny measure (am_code == 0) is a special case, since marginals may not be available (in "reweight" mode) */
      double score = (am_code == 0) ? _f[i] : AM(_f[i], _f1[row], _f2[col], N, sparse);
      _scores[i] = (transform_code) ? transform(score, transform_code) : score;
      i++;
    }
  }

  return scores;
}

// [[Rcpp::export]]
NumericVector CPP_dsm_score_sparse(int nr, int nc, IntegerVector p, IntegerVector row_of, NumericVector f, NumericVector f1, NumericVector f2, double N, int am_code, int sparse, int transform_code) {
  if (am_code < 0 || am_code >= am_table_entries)
    stop("internal error -- invalid AM code");
  am_func AM = am_table[am_code]; /* selected association measure */

  // -- don't check whether sparse=TRUE, so power users can compute non-sparse AMs for nonzero entries of the sparse matrix
  //  if (!sparse) stop("only sparse association scores can be used with sparse matrix representation");
  
  int n_items = f.size();
  NumericVector scores(n_items);
  if (am_code != 0 && (nr != f1.size() || nc != f2.size()))
    stop("internal error -- marginal vectors f1 and f2 not conformable with matrix f");

  IntegerVector::iterator _p = p.begin();
  IntegerVector::iterator _row_of = row_of.begin();
  NumericVector::iterator _f = f.begin();
  NumericVector::iterator _f1 = f1.begin();
  NumericVector::iterator _f2 = f2.begin();
  NumericVector::iterator _scores = scores.begin();
  
  for (int col = 0; col < nc; col++) {
    for (int i = _p[col]; i < _p[col+1]; i++) {
      /* frequeny measure (*am_code == 0) is a special case, since marginals may not be available ("reweight" mode) */
      double score = (am_code == 0) ? _f[i] : AM(_f[i], _f1[_row_of[i]], _f2[col], N, sparse);
      _scores[i] = (transform_code) ? transform(score, transform_code) : score;
    }
  }

  return scores;
}
/*
 *  Transform similarity values to distances (optionally as in-place operation on dense matrix)
 */

/* internal codes for type of transformation
 *  0 = cosine -> angle
 *  1 = sim -> d = 1 - sim
*/

// [[Rcpp::export]]
NumericMatrix CPP_similarity_to_distance(NumericMatrix M, int opcode, double tol, bool duplicate = true) {
  if (!R_FINITE(opcode) || opcode < 0 || opcode > 1)
    stop("internal error -- invalid transformation method code");

  unsigned int n_items = M.length();
  NumericMatrix res = M;
  if (duplicate)
    res = clone(M);
  NumericMatrix::iterator _res = res.begin();

  unsigned int n_clamped = 0;
  for (unsigned int i = 0; i < n_items; i++) {
    double x = _res[i];
    switch (opcode) {
    case 0:
      if (x < -(1-tol)) {
        if (x < -(1+tol)) n_clamped++;
        x = -1;
      }
      else if (x > (1-tol)) {
        if (x > (1+tol)) n_clamped++;
        x = 1;
      }
      x = acos(x) * 180 / M_PI; 
      break;
    case 1:
      x = 1 - x;
      break;
    }
    _res[i] = x;
  }
  
  if (n_clamped > 0)
    Rf_warning("angular distance may be inaccurate (some cosine values out of range)");

  return res;
}
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
