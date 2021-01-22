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

