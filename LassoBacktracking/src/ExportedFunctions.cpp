#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// Create a new matrix X that is a scaled and centred version of x. Modify scales and centres
// If there are columns in x that are too close to the intercept these are removed in X
// p is passed as a vector so it can be modified in place (and reduced if there are constant columns)
// var_names_main maps the new columns in X to the old columns in x
// [[Rcpp::export]]
NumericMatrix scale_cen(NumericMatrix x, NumericVector scales, NumericVector centres,
    IntegerVector p, int p_max, IntegerVector var_names_main) {
  double thresh = 1e-10; // threshold for the variance of a constant columns
  NumericMatrix X(x.nrow(), p_max);
  NumericVector centres_temp(p[0]);
  vector<int> var_names;
  var_names.reserve(p[0]);
  unsigned int k=0;
  // Determine centres and scales first
  for (int j=0; j<x.ncol(); j++) {
    for (int i=0; i<x.nrow(); i++) {
      centres_temp[j] += x(i,j);
    }
    centres_temp[j] /= x.nrow();
  }
  for (int j=0; j<x.ncol(); j++) {
    for (int i=0; i<x.nrow(); i++) {
      double temp = x(i,j) - centres_temp[j];
      scales[k] += temp*temp;
    }
    if (abs(scales[k]) < thresh) {
      scales[k]=0;
    } else {
      var_names.push_back(j);
      scales[k] = sqrt(scales[k]);
      for (int i=0; i<X.nrow(); i++) {
        X(i, k) = (x(i,j) - centres_temp[j])/scales[k];
      }
      k++;
    }
  }
  p[0]=k;
  for (k=0; k<var_names.size(); k++) {
    centres[k] = centres_temp[var_names[k]];
    var_names_main[k] = var_names[k] + 1;
  }
  return X;
}

// [[Rcpp::export]]
int change_dim(IntegerVector dim, int p_eff) {
  dim[0] = p_eff;
  return 0;
}

// beta[, l_start:ncol(beta)] <- 0
// [[Rcpp::export]]
int zero(NumericMatrix beta, int l_start) {
  for (int j=l_start-1; j<beta.ncol(); j++) {
    for (int i=0; i<beta.nrow(); i++) {
      beta(i, j) = 0;
    }
  }
  return 0;
}

// return components of a that are TRUE in b
// Same as a[b[a]]
// return components of a that are TRUE in b
// [[Rcpp::export]]
IntegerVector in_log(IntegerVector a, LogicalVector b) {
  vector<int> out;
  out.reserve(a.size());
  for (int j=0; j<a.size(); j++) {
    if (b[a[j]-1]) out.push_back(a[j]);
  }
  return wrap(out);
}

// any(violations[1:p_eff])
// [[Rcpp::export]]
bool any_indmax(LogicalVector violations, int p_eff) {
  for (int j=0; j<p_eff; j++) {
    if (violations[j]) return true;
  }
  return false;
}

// which(violations[1:p_eff])
// [[Rcpp::export]]
IntegerVector which_indmax(LogicalVector violations, int p_eff) {
  vector<int> out;
  for (int j=0; j<p_eff; j++) {
    if (violations[j]) out.push_back(j+1);
  }
  return wrap(out);
}

// returns true if (inter1, inter2) is "less than" interactions[, col]
bool check_less(int inter1, int inter2, int col, IntegerMatrix interactions) {
  if (inter1 < interactions(0, col)) return true;
  if ((inter1 == interactions(0, col)) & (inter2 < interactions(1, col))) return true;
  return false;
}
bool check_more(int inter1, int inter2, int col, IntegerMatrix interactions) {
  if (inter1 > interactions(0, col)) return true;
  if ((inter1 == interactions(0, col)) & (inter2 > interactions(1, col))) return true;
  return false;
}

// returns true if (inter1, inter2) is NOT in interactions
bool check_inter(int inter1, int inter2, IntegerMatrix interactions) {
  // interactions will be sorted by first row and then by second row
  // first row always smaller
  if (interactions.ncol() == 0) return true;
  if (inter1 > inter2) {
    int temp = inter2;
    inter2 = inter1;
    inter1 = temp;
  }
  int lower=0;
  int upper=interactions.ncol()-1;
  int middle=upper/2; // floor

  if (check_less(inter1, inter2, 0, interactions)) return true;
  if (check_more(inter1, inter2, upper, interactions)) return true;
  
  while (middle - lower > 1) {
    if (check_less(inter1, inter2, middle, interactions)) {
      upper=middle;
    } else {
      lower=middle;
    }
    middle=(upper+lower)/2; // floor
  }
  if ((inter1 == interactions(0, middle)) && (inter2 == interactions(1, middle))) return false;
  return true;
}

// [[Rcpp::export]]
int add_inter_orig(NumericMatrix X, NumericVector scales, NumericVector centres, int p_eff,
      IntegerMatrix inter_orig) {
  double scale_prod;
  for (int inter_ind=0; inter_ind<inter_orig.ncol(); inter_ind++) {
    for (int i=0; i<X.nrow(); i++) {
      X(i, p_eff) = X(i, inter_orig(0, inter_ind)) * X(i, inter_orig(1, inter_ind));
      centres[p_eff] += X(i, p_eff);
    }
    centres[p_eff] /= X.nrow();
    for (int i=0; i<X.nrow(); i++) {
      X(i, p_eff) -= centres[p_eff];
      scales[p_eff] += X(i, p_eff)*X(i, p_eff);
    }
    scale_prod = scales[inter_orig(0, inter_ind)]*scales[inter_orig(1, inter_ind)];
    centres[p_eff] *= scale_prod;
    centres[p_eff] -= centres[inter_orig(0, inter_ind)]*centres[inter_orig(1, inter_ind)];
    scales[p_eff] = sqrt(scales[p_eff]);
    for (int i=0; i<X.nrow(); i++) {
      X(i, p_eff) /= scales[p_eff];
    }
    scales[p_eff] *= scale_prod;
    p_eff++;
  }
  return p_eff;
}

// [[Rcpp::export]]
int add_inter(NumericMatrix X, IntegerVector active_old, IntegerVector active_new,
      NumericVector scales, NumericVector centres, int p_eff, IntegerMatrix interactions, int p,
      IntegerMatrix inter_orig) {
  int int_eff = p_eff - p;
  IntegerVector active_old0 = active_old - 1;
  IntegerVector active_new0 = active_new - 1;
  double scale_prod;
  for (int j=0; j<active_old.size(); j++) {
    for (int k=0; k<active_new.size(); k++) {
      if (check_inter(active_old[j], active_new[k], inter_orig)) {
        interactions(0, int_eff) = active_old[j];
        interactions(1, int_eff) = active_new[k];
        for (int i=0; i<X.nrow(); i++) {
          X(i, p_eff) = X(i, active_old0[j]) * X(i, active_new0[k]);
          centres[p_eff] += X(i, p_eff);
        }
        centres[p_eff] /= X.nrow();
        for (int i=0; i<X.nrow(); i++) {
          X(i, p_eff) -= centres[p_eff];
          scales[p_eff] += X(i, p_eff)*X(i, p_eff);
        }
        scale_prod = scales[active_old0[j]]*scales[active_new0[k]];
        centres[p_eff] *= scale_prod;
        centres[p_eff] -= centres[active_old0[j]]*centres[active_new0[k]];
        scales[p_eff] = sqrt(scales[p_eff]);
        for (int i=0; i<X.nrow(); i++) {
          X(i, p_eff) /= scales[p_eff];
        }
        scales[p_eff] *= scale_prod;
        p_eff++;
        int_eff++;
      }
    }
  }
  for (int j=0; j<active_new.size(); j++) {
    for (int k=0; k<j; k++) {
      if (check_inter(active_new[j], active_new[k], inter_orig)) {
        interactions(0, int_eff) = active_new[j];
        interactions(1, int_eff) = active_new[k];
        for (int i=0; i<X.nrow(); i++) {
          X(i, p_eff) = X(i, active_new0[j]) * X(i, active_new0[k]);
          centres[p_eff] += X(i, p_eff);
        }
        centres[p_eff] /= X.nrow();
        for (int i=0; i<X.nrow(); i++) {
          X(i, p_eff) -= centres[p_eff];
          scales[p_eff] += X(i, p_eff)*X(i, p_eff);
        }
        scale_prod = scales[active_new0[j]]*scales[active_new0[k]];
        centres[p_eff] *= scale_prod;
        centres[p_eff] -= centres[active_new0[j]]*centres[active_new0[k]];
        scales[p_eff] = sqrt(scales[p_eff]);
        for (int i=0; i<X.nrow(); i++) {
          X(i, p_eff) /= scales[p_eff];
        }
        scales[p_eff] *= scale_prod;
        p_eff++;
        int_eff++;
      }
    }
  }
  return p_eff;
}

// abs(t(x[, I]) %*% y[, l])
// [[Rcpp::export]]
NumericVector inner_prod_abs_comp(NumericMatrix x, IntegerVector I, NumericMatrix y, int l0) {
  NumericVector out(I.size());
  for (int k=0; k<I.size(); k++) {
    int j = I[k]-1;
    for (int i=0; i<x.nrow(); i++) {
      out[k] += x(i, j) * y(i, l0);
    }
    out[k] = abs(out[k]);
  }
  return out;
}

// any(abs(as.numeric(t(x[, (p_eff_old+1):p_eff]) %*% y)) > lambda[l])
bool violates(NumericMatrix x, NumericMatrix y, int p_eff, int p_eff_old,
    int l0, NumericVector lambda) {
  for (int k=p_eff_old; k<p_eff; k++) {
    double inn_prod=0;
    for (int i=0; i<x.nrow(); i++) {
      inn_prod += x(i, k) * y(i, l0);
    }
    if (abs(inn_prod) > lambda[l0]) return true;
  }
  return false;
}

// abs(as.numeric(t(x[, 1:p_eff][, I[1:p_eff]]) %*% y[, l]))
// [[Rcpp::export]]
NumericVector inner_prod_abs_comp2(NumericMatrix x, LogicalVector I, NumericMatrix y, int p_eff, int l0) {
  vector<double> out;
  out.reserve(I.size());
  for (int k=0; k<p_eff; k++) {
    if (I[k]) {
      double temp=0;
      for (int i=0; i<x.nrow(); i++) {
        temp += x(i, k) * y(i, l0);
      }
      out.push_back(abs(temp));
    }
  }
  return wrap(out);
}

// any(v[indices])
// [[Rcpp::export]]
bool any_ind(LogicalVector v, IntegerVector indices) {
  for (int i=0; i<indices.size(); i++) {
    if (v[indices[i]-1]) return true;
  }
  return false;
}

// 0 indexed
double inner_prod_ind(NumericMatrix x, NumericVector y, int index) {
  // deepest loop
  double out=0;
  for (int i=0; i<x.nrow(); i++) {
    out += x(i, index)*y[i];
  }
  return out;
}

inline double soft_thresh(double u, double t) {
  if (u > 0) {
    if (u > t) {
      return u-t;
    }
  } else {
    if (u < -t) {
      return u+t;
    }
  }
  return 0;
}

// [[Rcpp::export]]
int beta_active(NumericMatrix x, NumericMatrix beta, NumericMatrix resid_cur,
                IntegerVector active_set, int l, double thresh, int maxit,
                double lam, double alpha_lam, double alpha, double alpha_lam_div,
                double null_dev) {
                  // l will be 0 indexed
  double new_coef;
  double old_coef;
  double diff;
  double sq_diff;
  double abs_diff;
  double rel_err;
  int iter=0;
  double obj_change;
  double resid_x_inn_prod;
  NumericVector temp_resid(x.nrow());
  IntegerVector active_set0 = active_set-1; // Rcpp sugar
  do {
    rel_err=0;
    for (int k=0; k<active_set0.size(); k++) {
      // inner_prod_ind is at deepest loop level
      for (int i=0; i<x.nrow(); i++) {
        temp_resid[i] = resid_cur(i, l) + x(i, active_set0[k])*beta(active_set0[k], l);
      }
      resid_x_inn_prod = inner_prod_ind(x, temp_resid, active_set0[k]);
      new_coef = soft_thresh(resid_x_inn_prod/ alpha_lam_div, alpha_lam);
      old_coef = beta(active_set0[k], l);
      diff = old_coef - new_coef;
      abs_diff = abs(new_coef) - abs(old_coef);
      sq_diff = -(old_coef + new_coef)*diff;
      obj_change = 2*resid_x_inn_prod*diff + alpha_lam_div*sq_diff + alpha_lam*abs_diff;
      if (abs(diff) > 0) {
        // Perhaps dangerous may need machine epsilon
        // We do not use sugar here for increased speed
        for (int i=0; i<x.nrow(); i++) {
          resid_cur(i, l) += x(i, active_set0[k])*diff;
        }
        rel_err = max(obj_change / null_dev, rel_err);
        beta(active_set0[k], l) = new_coef;
      }
    }
    iter++;
  } while (rel_err > thresh && iter < maxit);
  return 0;
}

// [[Rcpp::export]]
int find_l0(NumericMatrix X, int p_eff_old, int p_eff, NumericMatrix resid_cur, int l0,
      NumericVector lambda) {
  // first check current l0
  if (!violates(X, resid_cur, p_eff, p_eff_old, l0, lambda)) return l0+1;
  // then check l0=0.  If this fails we return 0
  if (violates(X, resid_cur, p_eff, p_eff_old, 0, lambda)) return 0;
  // then bisection search
  // lower never violates
  int lower=0;
  int upper=l0;
  int middle=l0/2; // floor
  while (middle - lower > 1) {
    if (violates(X, resid_cur, p_eff, p_eff_old, middle, lambda)) {
      upper=middle;
    } else {
      lower=middle;
    }
    middle=(upper+lower)/2;
  }
  return middle;
}
