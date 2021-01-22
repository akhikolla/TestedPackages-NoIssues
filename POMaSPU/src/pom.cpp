# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;




// [[Rcpp::export]]
double count_if(arma::uvec x) {
    double counter = 0;
    for(int i = 0; i < x.size(); i++) {
        if(x[i] == TRUE) {
            counter++;
        }
    }
    return counter;
}







// [[Rcpp::export]]
arma::vec rankC(NumericVector x) {
//   if (is_true(any(duplicated(x)))) {
//     Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
//   }
  NumericVector sorted = clone(x).sort();
  return as<arma::vec>(wrap(match(x, sorted)));
}






// [[Rcpp::export]]
List permhigh(arma::mat Ts, arma::mat rsum, arma::mat X, arma::mat po, int nperm) {
    
    const int npo = po.n_rows;
    arma::vec Pspu(npo);
    arma::vec minp0(nperm);

    for (int q = 0; q < npo; q++){
        double w = po(q,0);
        arma::vec T0s(nperm);
        T0s.fill(0);

         for (int i = 0; i < nperm; i++){
            arma::mat U0 = rsum.t() * shuffle(X, 0);
            if (w == 0) {
                T0s(i) = abs(U0).max();
            } else {
                T0s(i) = accu(pow(U0, w));
            }
        }
        Pspu(q) = count_if(std::abs(Ts(q)) < abs(T0s))/nperm;
        arma::vec P0s = (nperm - rankC(as<NumericVector>(wrap(abs(T0s)))))/(nperm-1);
        if (q==0){
            minp0 = P0s;
        } else {
            minp0(find(minp0 > P0s)) = P0s(find(minp0 > P0s));
        }
    }
    double PaSPU = count_if(minp0 < Pspu.min())/nperm;
    List A ;
    A["P_SPU"] = Pspu;
    A["P_aSPU"] = PaSPU;
    return(A) ;
}

