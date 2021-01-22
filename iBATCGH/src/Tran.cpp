#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]

arma::mat Tran(arma::mat xi){

int m=xi.n_cols;
int s=xi.n_rows;


arma::mat out(4,4);
out.zeros();
for (int i=0; i<s; i++){
    for(int j=0; j<(m-1);j++){
            out(xi(i,j)-1,xi(i,j+1)-1)++;
    }
}

out.row(0)=out.row(0)/sum(out.row(0));
out.row(1)=out.row(1)/sum(out.row(1));
out.row(2)=out.row(2)/sum(out.row(2));
out.row(3)=out.row(3)/sum(out.row(3));

return out;
}
