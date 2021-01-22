void compnoreg_sampleloglambda(double* lambda,
double alpha0, double lambda01,
 double alpha1, 
 double alpha2,double lambda2,  double p,
 double* t,
 int* delta, int size);

void compnoreg_samplealpha(double* alpha, double lambda1,
double alpha2, double lambda2,
double p,
double alphaalpha, double alphalambda,
 double* t,
 int* delta, int size);

void compnoreg_samplep(double* p, double gamma0, double gamma1,double *t, int* delta,
 double alpha1, double lambda1,double alpha2, double lambda2, int size);
