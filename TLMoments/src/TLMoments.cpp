#include <Rcpp.h>
using namespace Rcpp;

/////////////////////////////////////////////////////////////////////
// COMPUTATION.METHOD = PWM

// declare pwm_C
double pwm_C(NumericVector x, int r);

// [[Rcpp::export]]
double z_C(int r, int k, int s, int t) {
  double out = Rcpp::internal::factorial(r) * Rcpp::internal::factorial(r+s+t+1) / ((r+1) * Rcpp::internal::factorial(r+s) * Rcpp::internal::factorial(r+t)) * std::pow(-1.0, s+r+k) * R::choose(r+t, k-s) * R::choose(r+k, r);
  return out;
}

// [[Rcpp::export]]
NumericMatrix Z_C(int maxr, int s, int t) {
  NumericMatrix out(maxr, maxr+s+t);

  for (int r = 1; r <= maxr; r++) {
    for (int k = s; k <= r-1+s+t; k++) {
      out(r-1, k) = z_C(r-1, k, s, t);
    }
  }
  return out;
}

// [[Rcpp::export]]
double TLMoment_PWM(NumericVector x, int r, int s, int t) {
  double sum = 0;
  for (int k = s; k <= r-1+s+t; k++) {
    sum += z_C(r-1, k, s, t) * pwm_C(x, k);
  }

  return sum;
}

// [[Rcpp::export]]
NumericVector PWM_to_TLMoments(NumericVector pwm, int s, int t) {
  int maxr = pwm.size() - s - t;
  NumericVector out(maxr);
  double sum;

  for (int r = 1; r <= maxr; r++) {
    sum = 0;
    for (int k = s; k <= r-1+s+t; k++) {
      sum += z_C(r-1, k, s, t) * pwm[k];
    }
    out[r-1] = sum;
  }

  return out;
}


/////////////////////////////////////////////////////////////////////
// COMPUTATION.METHOD = DIRECT

NumericVector w_direct(int r, int n, double s, double t) {
  NumericVector out(n);
  double sum;
  double add;

  if ((s == 0.0) & (t == 0.0)) { // L-Momente

    for (int i = 1; i <= n; i++) {
      sum = 0;
      for (int k = 0; k < r; k++) {
        sum += std::pow(-1.0, k) * R::choose(r-1, k) * R::choose(i-1, r+s-k-1) * R::choose(n-i, t+k);
      }
      out[i-1] = sum / r / R::choose(n, r+s+t);
    }

  } else {

    for (int i = 1; i <= n; i++) {
      sum = 0;
      for (int k = 0; k < r; k++) {
        add = std::pow(-1.0, k) / r / i / (n-i+1) / R::beta(k+1, r-k) / R::beta(r+s-k, i-r-s+k+1) / R::beta(t+k+1, n-i-t-k+1);
        if (NumericVector::is_na(add)) {
          add = 0;
        }
        sum += add;
      }
      out[i-1] = sum / r * (n+1) * R::beta(r+s+t+1, n-r-s-t+1);
    }

  }

  // NaN zu Null und außerhalb der Grenzen zu Null (damit Trimmung eingehalten wird)
  for (int i = 1; i <= n; i++) {
    if ((i < s) | (i > n-t+1) | NumericVector::is_na(out[i-1])) {
      out[i-1] = 0;
    }
  }

  return out;
}

// [[Rcpp::export]]
double TLMoment_direct(NumericVector x, int r, double s, double t) {
  double out;
  int n = x.size();
  NumericVector xs(clone(x));
  xs = xs.sort();

  NumericVector weights = w_direct(r, n, s, t);
  out = sum(xs * weights);
  return out;
}

/////////////////////////////////////////////////////////////////////
// COMPUTATION.METHOD = RECURSIVE

NumericMatrix W_recursive(int maxr, int n, double s, double t) {
  NumericMatrix out(maxr, n);

  if ((s == 0.0) & (t == 0.0)) { // L-Momente


    for (int r = 1; r <= maxr; r++) {
      if (r == 1) {
        for (int k = 1; k <= n; k++) {
          out(0, k-1) = 1;
        }
      } else if (r == 2) {
        IntegerVector k = seq_len(n);
        out(1, _) = (2*(NumericVector)k-n-1)/(n-1);
      } else {
        for (int k = 1; k <= n; k++) {
          out(r-1, k-1) = ((2*((double)r-1)-1)*(2*(double)k-n-1)*out(r-2, k-1) - ((double)r-2)*(n+(double)r-2)*out(r-3, k-1)) / (((double)r-1) * (n-(double)r+1));
        }
      }
    }
    for (int r = 1; r <= maxr; r++) {
      out(r-1, _) = out(r-1, _) / (double)n;
    }

  } else { // TL-Momente


    for (int r = 1; r <= maxr; r++) {
      if (r == 1) {
        for (int k = 1; k <= n; k++) {
          // Verwende entweder Gamma-Funktionen oder die Beta-Funktion
          //           if (k-s > 0 & n-k-t > 0) {
          //             out(0, k-1) = (n+1) / k / (n-k+1) / R::beta(s+1, k-s) / R::beta(t+1, n-k-t+1) * R::beta(s+t+2, n-s-t);
          //           } else {
          out(0, k-1) = R::gammafn(k)/(R::gammafn(s+1)*R::gammafn(k-s)) * R::gammafn(n-k+1)/(R::gammafn(t+1)*R::gammafn(n-k-t+1)) / (R::gammafn(n+1)/(R::gammafn(s+t+2)*R::gammafn(n-s-t)));
          // }
        }
      } else if (r == 2) {
        for (int k = 1; k <= n; k++) {
          out(1, k-1) = (2+s+t) * ((2+s+t)*k - (1+s)*(n+1)) / (2 * (1+s) * (1+t) * (n-s-t-1)) * out(0, k-1);
        }
      } else {
        double A = (r-1+s+t)*(r-1+s)*(n-r+1-s-t) / ((2*(r-1)+s+t-1)*(2*(r-1)+s+t));
        double C = (r-2)*(r-2+t)*(n+r-2) / ((2*(r-1)+s+t-2)*(2*(r-1)+s+t-1));
        for (int k = 1; k <= n; k++) {
          out(r-1, k-1) = ((k-s-1-A-C) * out(r-2, k-1) - (r-2)*(r-1+s+t)/((r-1)*(r-2+t)) * C * out(r-3, k-1)) / ((r*(r-1+t))/((r-1)*(r+s+t)) * A);
        }
      }
    }

  }

  // NaN zu Null und außerhalb der Grenzen zu Null (damit Trimmung eingehalten wird)
  for (int r = 1; r <= maxr; r++) {
    for (int i = 1; i <= n; i++) {
      if ((i < s) | (i > n-t+1) | NumericVector::is_na(out(r-1,i-1))) {
        out(r-1, i-1) = 0;
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector TLMoments_recursive(NumericVector x, int maxr, double s, double t) {
  NumericVector out(maxr);
  int n = x.size();
  NumericVector xs(clone(x));
  xs = xs.sort();

  NumericMatrix weights = W_recursive(maxr, n, s, t);
  for (int r = 1; r <= maxr; r++) {
    out[r-1] = sum(xs * weights(r-1, _));
  }

  return out;
}


/////////////////////////////////////////////////////////////////////
// COMPUTATION.METHOD = TRIMMING-DEGREE RECURRENCE

// [[Rcpp::export]]
NumericVector TLMoments_recurrence(NumericVector x, int maxr, double s, double t) {
  double fs = floor(s);
  double ft = floor(t);
  NumericVector l = TLMoments_recursive(x, maxr+(int)fs+(int)ft, s-fs, t-ft);

  for (double k = 1+t-ft; k <= t; k++) {
    for (int i = 1; i <= maxr+s+t-k; i++) {
      l[i-1] = (((double)i+k)*l[i-1] - ((double)i+1)*l[i]) / (2*(double)i + k - 1);
    }
  }
  for (double k = 1+s-fs; k <= s; k++) {
    for (int i = 1; i <= maxr+s-k; i++) {
      l[i-1] = (((double)i+k+t)*l[i-1] + 1/(double)i * ((double)i+1)*((double)i+t)*l[i]) / (2*(double)i + k + t - 1);
    }
  }

  return l[Range(0, maxr-1)];
}
