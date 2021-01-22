#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
const double eps = DBL_EPSILON;
double pcbinom(double x, double n, double p){
  if (x < eps) return(R_NegInf);
  if (x > n + 1 - eps) return(0);
  return(R::pbeta(p, x, n - x + 1, 0, 1));
}

// [[Rcpp::export]]
NumericVector pcbinomC(NumericVector q, NumericVector sz, NumericVector prob,
  bool logp){
  // assumes lower tail
  int n = int(q.size());
  n = std::max(n, int(sz.size()));
  n = std::max(n, int(prob.size()));
  NumericVector p(n);
  LogicalVector qnan = is_nan(q), qna = is_na(q);
  LogicalVector sznan = is_nan(sz), szna = is_na(sz);
  LogicalVector probnan = is_nan(prob), probna = is_na(prob);
  for (int i = 0; i < n; i++){
    if (qna[i%q.size()] || szna[i%sz.size()] || probna[i%prob.size()]){
      if ((qna[i%q.size()] && !qnan[i%q.size()]) ||
          (szna[i%sz.size()] && !sznan[i%sz.size()]) ||
          (probna[i%prob.size()] && !probnan[i%prob.size()])){
        p[i] = NA_REAL;
      } else {
        p[i] = R_NaN;
      }
      continue;
    }
    if (prob[i%prob.size()] < 0 || prob[i%prob.size()] > 1){
      p[i] = R_NaN;
      continue;
    }
    if (sz[i%sz.size()] < 0 || sz[i%sz.size()] == R_PosInf){
      p[i] = R_NaN;
      continue;
    }
    if (q[i%q.size()] > sz[i%sz.size()] + 1){
      if (logp){
        p[i] = 0;
      } else {
        p[i] = 1;
      }
      continue;
    }
    if (q[i%q.size()] < 0){
      if (logp){
        p[i] = R_NegInf;
      } else {
        p[i] = 0;
      }
      continue;
    }
    p[i] = R::pbeta(prob[i%prob.size()], q[i%q.size()],
      sz[i%sz.size()] - q[i%q.size()] + 1, 0, logp);
  }
  return(p);
}

// [[Rcpp::export]]
NumericVector qcbinomC(NumericVector p, NumericVector m, NumericVector g,
  bool rcb = false){
  // assumes probabilities (p) are on logscale and calculates lower tail
  // this level of preprocessing should be done in R
  int n = int(p.size()), ctr = 0;
  if (!rcb){//if random, length = length(p); otherwise, length = maxlen(p, m, g)
    n = std::max(n, int(m.size()));
    n = std::max(n, int(g.size()));
  }
  double tol = sqrt(eps);
  NumericVector x(n);
  double* pp; // local copy of p
  pp = new double[n];
  double x0, x1, x2, f0, f1, f2, sz, prob;
  double maxx;
  double rtt = sqrt(eps);
  double hwid, dx, d0;
  double tmp[4];
  LogicalVector gnan = is_nan(g), gna = is_na(g);
  LogicalVector mnan = is_nan(m), mna = is_na(m);
  LogicalVector pnan = is_nan(p), pna = is_na(p);
    //calculate x[i]
  for (int i = 0; i < n; i++){ 
    // pre-process and error-check

    // NAs and NaNs:
    if (gna[i%g.size()] || mna[i%m.size()] || pna[i%p.size()]){
      if ((pna[i%p.size()] && !pnan[i%p.size()]) ||
          (mna[i%m.size()] && !mnan[i%m.size()]) ||
          (gna[i%g.size()] && !gnan[i%g.size()])){
        x[i] = NA_REAL;
        continue;
      }
      else {
        x[i] = R_NaN;
        continue;
      }
    }
    // bad numeric parameters
    if (g[i%g.size()] > 1 || g[i%g.size()] < 0 ||
        m[i%m.size()] < 0 || p[i%p.size()] > 0 || m[i%m.size()] == R_PosInf){
      x[i] = R_NaN;
      continue;
    }
    // special numeric parameters
    if (g[i%g.size()] >= 1 - 2 * eps){
      x[i] = m[i%m.size()] + 1.0;
      continue;
    }
    if (g[i%g.size()] < 2 * eps){
      x[i] = 0;
      continue;
    }
    if (p[i%p.size()] == R_NegInf){
      x[i] = 0;
      continue;
    }
    sz = m[i%m.size()];
    pp[i] = p[i%p.size()];
    prob = g[i%g.size()];
    if (pp[i] == 0){ // i.e., p = 1
      x[i] = sz + 1;
      continue;
    }

    // find initial bracket points
    maxx = sz + 1;
    // normal approximation for x0
    x0 = R::qnorm(pp[i], sz * prob + 0.5, sqrt(sz*(prob*(1 - prob))), 1, 1);
    if (x0 < 0){
      x0 = eps;
      f0 = pcbinom(x0, sz, prob) - pp[i];
      if (f0 > 0){
        x[i] = 0;
        continue;
      }
      x1 = 0.5;
      f1 = pcbinom(x1, sz, prob) - pp[i];
      ctr = 0;
      while (f1 < 0){
        x1 = fmin(x1 + 0.2, maxx);
        f1 = pcbinom(x1, sz, prob) - pp[i];
        ctr++;
        if (ctr > 100){
          x[i] = R_NaN;
          break;
        }
      }
      if (ctr > 100){
        continue;
      }
    } else { // x0 >= 0
      f0 = pcbinom(x0, sz, prob) - pp[i];
      if (f0 == 0){
        x[i] = x0;
        continue;
      }
      if (f0 < 0){ //then increment x1 upward to bracket x
        x1 = fmin(x0 + 0.2, maxx);
        f1 = pcbinom(x1, sz, prob) - pp[i];
        ctr = 0;
        while(f1 < 0){
          x0 = x1;
          f0 = f1;
          x1 = fmin(x1 + 0.2, maxx);
          f1 = pcbinom(x1, sz, prob) - pp[i];
          ctr++;
          if (ctr > 100){
            x[i] = R_NaN;
            break;
          }
        }
        if (ctr > 100){
          continue;
        }
      } else { // f0 > 0 and increment x0 downward to bracket x
        x1 = x0;
        f1 = f0;
        x0 = fmax(eps, x1 - 0.2);
        f0 = pcbinom(x0, sz, prob) - pp[i];
        ctr = 0;
        while(f0 > 0){
          x1 = x0;
          f1 = f0;
          x0 = fmax(eps, x0 - 0.2);
          f0 = pcbinom(x0, sz, prob) - pp[i];
          ctr++;
          if (ctr > 100){
            x[i] = R_NaN;
            break;
          }
        }
        if (ctr > 100){
          continue;
        }
      }
    }
    // Brent's algorithm
    x2 = x0;
    f2 = f0;
    dx = x1 - x0;
    d0 = dx;
    while(1){
      if (fabs(f2) < fabs(f1)){
        x0 = x1;
        x1 = x2;
        x2 = x0;
        f0 = f1;
        f1 = f2;
        f2 = f0;
      }

      tol = 2.0 * eps * fabs(x1) + rtt;
      hwid = 0.5 * (x2 - x1);

      if (fabs(hwid) <= tol || f1 == 0.0 ){
        break;
      }

      if (fabs(dx) < tol || fabs(f0) <= fabs(f1)){
        dx = hwid;
        d0 = dx;
      } else {
        tmp[3] = f1/f0;
        if (x0 == x2){
          tmp[0] = 2.0 * hwid * tmp[3];
          tmp[1] = 1.0 - tmp[3];
        } else {
          tmp[1] = f0/f2;
          tmp[2] = f1/f2;
          tmp[0] = tmp[3] * (2.0*hwid*tmp[1]*(tmp[1] - tmp[2]) -
              (x1 - x0)*(tmp[2] - 1.0));
          tmp[1] = (tmp[1] - 1.0) * (tmp[2] - 1.0) * (tmp[3] - 1.0);
        }
        if (0.0 < tmp[0]){
          tmp[1] = -tmp[1];
        } else {
          tmp[0] = -tmp[0];
        }
        tmp[3] = dx;
        dx = d0;
        if (2.0 * tmp[0] < 3.0 * hwid * tmp[1] - fabs(tol * tmp[1]) &&
            tmp[0] < fabs(0.5 * tmp[3] * tmp[1])){
          d0 = tmp[0]/tmp[1];
        } else {
          dx = hwid;
          d0 = dx;
        }
      }
      x0 = x1;
      f0 = f1;
      if (tol < fabs(d0)){
        x1 = x1 + d0;
      } else if (0.0 < hwid){
        x1 = x1 + tol;
      } else {
        x1 = x1 - tol;
      }

      f1 =  pcbinom(x1, sz, prob) - pp[i];
      if ((0.0 < f1 && 0.0 < f2) || (f1 <= 0.0 && f2 <= 0.0)){
        x2 = x0;
        f2 = f0;
        dx = x1 - x0;
        d0 = dx;
      }
    }
    x[i] = x1;
  }
  return(x);
}

// [[Rcpp::export]]
NumericMatrix dcblp(NumericVector x, NumericVector m, NumericVector g){
// intermediate calculations for dcbinom (called from R)
  int n = x.size();
  double xi, mi, gi;
  n = std::max(n, int(m.size()));
  n = std::max(n, int(g.size()));
  NumericMatrix f(n, 3);
  double h = 1e-6;
  for (int i = 0; i < n; i++){
    xi = x[i%x.size()];
    mi = m[i%m.size()];
    gi = g[i%g.size()];
    if (xi < 0) {
      f(i, 0) = R_NegInf;
      f(i, 1) = R_NegInf;
      f(i, 2) = h;
    } else if (xi <= eps){
      xi = eps * 1.0000001; // --> density at zero = density at eps
      f(i, 0) = pcbinom(xi + h, mi, gi);
      f(i, 1) = pcbinom(xi, mi, gi);
      f(i, 2) = h;
    } else if (xi > mi + 1) {
      f(i, 0) = 0;
      f(i, 1) = 0;
      f(i, 2) = h;
    } else if (xi >= h && xi <= mi + 1 - h){
      f(i, 0) = pcbinom(xi + h, mi, gi);
      f(i, 1) = pcbinom(xi - h, mi, gi);
      f(i, 2) = 2 * h;
    } else if (xi <= h) {
      f(i, 0) = pcbinom(xi + h, mi, gi);
      f(i, 1) = pcbinom(xi, mi, gi);
      f(i, 2) = h;
    } else {
      f(i, 0) = pcbinom(xi, mi, gi);
      f(i, 1) = pcbinom(xi - h, mi, gi);
      f(i, 2) = h;
    }
  }
  return(f);
}



