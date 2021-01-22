void compreg_sampleloglambda1(double* lambda,
double alpha0, double lambda01,
 double alpha1, 
 double alpha2,double lambda2,  double p,
 double* t,
 int* delta, int size, double* xbeta1, double* xbeta2);

void compreg_sampleloglambda2(double* lambda,
double alpha0, double lambda02,
 double alpha2, 
 double alpha1,double lambda1,  double p,
 double* t,
 int* delta, int size, double* xbeta1, double* xbeta2);

void compreg_samplealpha1(double* alpha, double lambda1,
double alpha2, double lambda2,
double p,
double alphaalpha, double alphalambda,
 double* t,
 int* delta, int size, double* xbeta1, double* xbeta2);

void compreg_samplealpha2(double* alpha, double lambda2,
double alpha1, double lambda1,
double p,
double alphaalpha, double alphalambda,
 double* t,
 int* delta, int size, double* xbeta1, double *xbeta2);

void compreg_samplep(double* p, double gamma0, double gamma1,double *t, int* delta,
 double alpha1, double lambda1,double alpha2, double lambda2, int size,
 double* xbeta1, double *xbeta2);

void compreg_samplebeta1(double* beta, double betasl,
 double alpha1,double lambda1,
 double alpha2, double lambda2,
 double* t,
 int* delta, int size, double*x, double *xbeta1left, double* xbeta2, double p) ;


void compreg_samplebeta2(double* beta, double betasl, double alpha1,double lambda1,
 double alpha2, double lambda2,
 double* t,
 int* delta, int size, double*x, double *xbeta2left, double* xbeta1, double p) ;

