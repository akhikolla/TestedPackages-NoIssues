#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(cpp)]]
double truncgamma(double a,double b,double lower){

double low=lower/b;
double x=0;
double u;
double u2;
if(a<1){
GetRNGstate();
x=low+Rf_rexp(1);
u=unif_rand();
PutRNGstate();
while((log(x)+log(u)/(1-a))>log(a)){
									GetRNGstate();
                                    x=low+Rf_rexp(1);
                                    u=unif_rand();
									PutRNGstate();
}
}

if(a==1){GetRNGstate();x=low+Rf_rexp(1);PutRNGstate();}
if(a>1){

     if(low>(a-1)){
		 GetRNGstate();
          u=Rf_rexp(1);
          u2=Rf_rexp(0.5);
		 PutRNGstate();
          x=low+u/(1-(a-1)/low);
          while((x/low-1+log(low/x))>(u2/(a-1))){
			  GetRNGstate();
                  u=Rf_rexp(1);
                  u2=Rf_rexp(1);
			  PutRNGstate();
                  x=low+u/(1-(a-1)/low);
          }
     }
     if(low<=(a-1)){
		 GetRNGstate();
          x=Rf_rgamma(a,1);
		 PutRNGstate();
          while(x<=low){
			  GetRNGstate();
                  x=Rf_rgamma(a,1);
			  PutRNGstate();
          }
     }
     
}
x=x*b;

return (x);
}

