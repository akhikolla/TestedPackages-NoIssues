#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(cpp)]]
arma::vec sm(arma::mat distance,arma::imat xi,double disfix)
{
#define DEN (exp(1.0)-1.0)

int nOss=xi.n_rows;
int m=xi.n_cols;

//calculate the results
arma::vec sum(m+1);
sum.zeros();
for (int i = 1; i < m; i++){
    arma::uvec stemp=find(xi.col(i)==xi.col(i-1));
    sum[i]=stemp.n_rows*((exp(1-(distance[i-1]/disfix))-1)/DEN)/nOss;
}

return sum;

}

