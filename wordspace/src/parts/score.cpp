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
