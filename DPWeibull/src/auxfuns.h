// The following includes 4 functions
// logsWeibloglambda
// logsWeib
// sWeibloglambda
// sWeib
// Footnote 1

double logsWeibloglambda(double t, double alpha, double loglambda);

double logsWeib(double t, double alpha, double lambda);

double sWeibloglambda(double t, double alpha, double loglambda);

double sWeib(double t, double alpha, double lambda);
// The following includes 4 functions
// logdWeibloglambda
// logdWeib
// dWeibloglambda
// dWeib
// Footnote2

double logdWeibloglambda(double t, double alpha, double loglambda);

double logdWeib(double t, double alpha, double lambda);

double dWeibloglambda(double t, double alpha, double loglambda);

double dWeib(double t, double alpha, double lambda);
// The following includes 2 functions
// dWeib
// logdWeib

double pWeib(double t, double alpha, double lambda);

double logpWeib(double t, double alpha, double lambda);

double F1(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta);
double F1v2(double t, double p, double alpha, double lambda, double xbeta);
double logf1(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta);
double logf1v2(double t, double p, double alpha, double lambda, double xbeta);
double f1(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta);


double F2(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta1,
 NumericVector beta2);

double F2v2(double t, double p, double alpha, double lambda, double xbeta1, double xbeta2);
double logf2(double t, double p, double alpha, double lambda,
 NumericVector x, NumericVector beta1, NumericVector beta2);

double logf2v2(double t, double p, double alpha, double lambda, double xbeta1, double xbeta2);

double f2(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta1,
 NumericVector beta2);

double f1v2(double t, double p, double alpha, double lambda, double xbeta);

double f2v2(double t, double p, double alpha, double lambda, double xbeta1, double xbeta2);

double logF1(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta);

double logF2(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta1,
 NumericVector beta2);

double logScomp(double t, double alpha1, double lambda1,
double alpha2, double lambda2, double xbeta1, double xbeta2, double p);

double timedWeibloglambda(double t, double ts, double alpha, double loglambda, double beta);

double timesWeibloglambda(double t, double ts, double alpha, double loglambda, double beta);

