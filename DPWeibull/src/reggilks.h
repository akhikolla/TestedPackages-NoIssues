void reg_sampleloglambda(double* lambda, double alpha,
double alpha0, double lambda0,
 double* tl, double* tr,
 int* delta, int* pi, int size, double* xbeta);
 


void reg_samplealpha(double* alpha, double lambda,
double alphaalpha, double alphalambda,
 double* tl, double* tr,
 int* delta, int* pi, int size, double* xbeta);
 


void reg_samplebeta(double* beta, double betasl, double alpha,
 double* tl, double* tr,
 int* delta, int* pi, int size, double*x, double *loglambda) ;
	

void reg_samplebeta2(double* beta, double betasl, double alpha,
 double* tl, double* tr,
 int* delta, int* pi, int size, double *loglambda) ;
