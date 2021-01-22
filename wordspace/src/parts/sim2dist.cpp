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
