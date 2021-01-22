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
