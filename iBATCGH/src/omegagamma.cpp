#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(cpp)]]
void omegagamma(arma::mat distance,arma::imat xi,double disfix,double alpha,arma::colvec *gammaex,arma::colvec *omega1ex,arma::colvec *omega2ex,arma::colvec *sumex)
{
#define DEN (exp(1.0)-1.0)

int nOss=xi.n_rows;
int m=xi.n_cols;

//Calculate the results
arma::colvec gamma(m);
arma::colvec omega1(m);
arma::colvec omega2(m);
arma::colvec sum(m+1);
sum.zeros();
for (int i = 1; i < m; i++){
    arma::uvec stemp=find(xi.col(i)==xi.col(i-1));
    sum[i]=stemp.n_rows*((exp(1-(distance[i-1]/disfix))-1)/DEN)/nOss;
}

for (int i = 0; i < m; i++){
    omega1[i]=sum[i]/(alpha+sum[i]+sum[i+1]);
    omega2[i]=sum[i+1]/(alpha+sum[i]+sum[i+1]);
    gamma[i]=1-(omega1[i]+omega2[i]);
}

*gammaex=gamma;
*omega1ex=omega1;
*omega2ex=omega2;
*sumex=sum;

}

