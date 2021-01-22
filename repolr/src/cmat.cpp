// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

arma::mat alphpow(double x, arma::mat mat){
  // construct matrix of powers of x
  arma::mat alphpowmat = arma::exp(std::log(x) * mat);
  return(alphpowmat);
}


// [[Rcpp::export]]

List cmat(Rcpp::NumericVector ctimes, double alpha, Rcpp::String corrmod, 
                   Rcpp::String diffmeth, double h){
    
   // input data
   arma::vec times = Rcpp::as<arma::vec>(ctimes);
    
   // sparse matrices 
   unsigned int ntimes = times.n_elem;
   double xntimes(ntimes);
   arma::sp_mat icmat(ntimes, ntimes), gicmat(ntimes, ntimes), ggicmat(ntimes, ntimes); 
   double ntimes1 = std::pow(xntimes - 1, -1.0);

   if(ntimes > 1){

   // two or more time points
   double lalpha, ualpha, lphi, uphi;
   lalpha = alpha; ualpha = alpha; lphi = alpha; uphi = alpha;

    if(corrmod == "ar1"){
     
     // construct icmat
     arma::sp_mat mlag1(ntimes, ntimes), mlag2(ntimes, ntimes);
     arma::mat diff1(1, ntimes), diff2(1, ntimes), lag1, lag2, a1, a2;
     arma::mat la1, ua1, la2, ua2;
     mlag1.diag(0) = -1 * arma::ones(ntimes);
     mlag1.diag(-1) = arma::ones(ntimes - 1);
     diff1 = times.t() * mlag1;
     lag1 = diff1.cols(0, ntimes - 2);
     a1 = arma::pow(1 - alphpow(alpha, 2 * lag1), -1.0);
     if(diffmeth == "numeric"){
      lphi = std::log(alpha) - std::log(1 - alpha) - h;
      uphi = std::log(alpha) - std::log(1 - alpha) + h;
      lalpha = std::exp(lphi) / (1 + std::exp(lphi));
      ualpha = std::exp(uphi) / (1 + std::exp(uphi));
      la1 = arma::pow(1 - alphpow(lalpha, 2 * lag1), -1.0);
      ua1 = arma::pow(1 - alphpow(ualpha, 2 * lag1), -1.0);
     }
     if(ntimes > 2){
      mlag2.diag(0) = -1 * arma::ones(ntimes);
      mlag2.diag(-2) = arma::ones(ntimes - 2);
      diff2 = times.t() * mlag2;
      lag2 = diff2.cols(0, ntimes - 3);
      a2 = arma::pow(1 - alphpow(alpha, 2 * lag2), -1.0);
      if(diffmeth == "numeric"){
       la2 = arma::pow(1 - alphpow(lalpha, 2 * lag2), -1.0);
       ua2 = arma::pow(1 - alphpow(ualpha, 2 * lag2), -1.0);
      }
     } else {
      a2.reset();
      if(diffmeth == "numeric"){
       la2.reset();
       ua2.reset();
      }
     }
     arma::mat diagnum = arma::join_rows(arma::join_rows(arma::ones(1), arma::pow(a2, -1.0)), arma::ones(1));
     arma::mat diagdenom1 = arma::trans(arma::repmat(arma::pow(a1.t(), -1.0), 1, 2));
     diagdenom1.reshape(1, 2 *(ntimes - 1));
     arma::mat diagdenom = arma::join_rows(arma::join_rows(arma::ones(1), diagdenom1), arma::ones(1));
     diagdenom.reshape(2, ntimes);
     icmat.diag(0) = diagnum / arma::prod(diagdenom, 0);
     icmat.diag(1) = - a1 % alphpow(alpha, lag1);
     icmat.diag(-1) = - a1 % alphpow(alpha, lag1);


     if(diffmeth == "analytic"){

      // construct gicmat
      arma::mat gicmatodiag = - arma::pow(a1, 2.0) % lag1 % alphpow(alpha, lag1 - 1) % (1 + alphpow(alpha, 2 * lag1));
      arma::mat gicmatediag = 2 * arma::pow(a1, 2.0) % lag1 % alphpow(alpha, 2 * lag1 - 1);
      arma::mat gicmatcdiag;
      if(ntimes > 2){
       gicmatcdiag = gicmatediag.cols(0, ntimes - 3) % (a1.cols(1, ntimes - 2) / a2)
                       - 2 * a1.cols(0, ntimes - 3) % arma::pow(a1.cols(1, ntimes - 2), 2.0)
                       % alphpow(alpha, 2 * lag1.cols(1, ntimes - 2) + 2 * lag2 - 1)
                       % (lag2 % (alphpow(alpha, - 2 * lag1.cols(1, ntimes - 2)) - 1)
                       - lag1.cols(1, ntimes - 2) % (alphpow(alpha, - 2 * lag2) - 1));
      } else {
       gicmatcdiag.reset();
      }
      gicmat.diag(0) = arma::join_rows(arma::join_rows(gicmatediag.cols(0, 0), gicmatcdiag), 
                                                              gicmatediag.cols(ntimes - 2,ntimes - 2));
      gicmat.diag(1) = gicmatodiag;
      gicmat.diag(-1) = gicmatodiag;

      // construct ggicmat
      arma::mat ggicmatodiag = - arma::pow(a1, 3.0) % lag1 % (alphpow(alpha, lag1 - 2)) % ((lag1 + 1) % alphpow(alpha, 4 * lag1)
                                       + 6 * lag1 % alphpow(alpha, 2 * lag1) + (lag1 - 1));
      arma::mat ggicmatediag = 2 * arma::pow(a1, 3.0) % lag1 % (alphpow(alpha, 2 * lag1 - 2))
                                             % ((2 * lag1 - 1) + (2 * lag1 + 1) % alphpow(alpha, 2 * lag1));
      arma::mat ggicmatcdiag;
      if(ntimes > 2){
       ggicmatcdiag = ggicmatediag.cols(0, ntimes - 3) % (a1.cols(1, ntimes - 2) / a2)
                       - 4 * gicmatediag.cols(0, ntimes - 3) % arma::pow(a1.cols(1, ntimes - 2), 2.0) 
                       % alphpow(alpha, 2 * lag1.cols(1, ntimes - 2) + 2 * lag2 - 1)
                       % (lag2 % (alphpow(alpha, - 2 * lag1.cols(1, ntimes - 2))- 1) 
                       - lag1.cols(1, ntimes - 2) % (alphpow(alpha, -2 * lag2) - 1))
                       - 2 * a1.cols(0, ntimes - 3) % arma::pow(a1.cols(1, ntimes - 2), 3.0) 
                       % alphpow(alpha, 2 * lag1.cols(1, ntimes - 2) + 2 * lag2 - 2)
                       % (lag2 % (2 * lag2 - 1) % (1 + (4 * lag1.cols(1, ntimes - 2) / (2 * lag2 - 1) - 1) 
                       % alphpow(alpha, 2 * lag1.cols(1, ntimes - 2))) % (alphpow(alpha, - 2 * lag1.cols(1, ntimes - 2)) - 1) 
                       - lag1.cols(1, ntimes - 2) % (2 * lag1.cols(1, ntimes - 2) - 1)
                       % (1 + (4 * lag1.cols(1, ntimes - 2) / (2 * lag1.cols(1, ntimes - 2) - 1) - 1) 
                       % alphpow(alpha, 2 * lag1.cols(1, ntimes - 2))) % (alphpow(alpha, -2 * lag2) - 1));

      } else {
       ggicmatcdiag.reset();
      }
      ggicmat.diag(0) = arma::join_rows(arma::join_rows(ggicmatediag.cols(0, 0), ggicmatcdiag), 
                                                              ggicmatediag.cols(ntimes - 2, ntimes - 2));
      ggicmat.diag(1) = ggicmatodiag;
      ggicmat.diag(-1) = ggicmatodiag;

      } else if(diffmeth == "numeric"){

      // upper and lower icmat
      arma::mat ldiagnum = arma::join_rows(arma::join_rows(arma::ones(1), arma::pow(la2, -1.0)), arma::ones(1));
      arma::mat ldiagdenom1 = arma::trans(arma::repmat(arma::pow(la1.t(), -1.0), 1, 2));
      ldiagdenom1.reshape(1, 2 *(ntimes - 1));
      arma::mat ldiagdenom = arma::join_rows(arma::join_rows(arma::ones(1), ldiagdenom1), arma::ones(1));
      ldiagdenom.reshape(2, ntimes);
      arma::mat udiagnum = arma::join_rows(arma::join_rows(arma::ones(1), arma::pow(ua2, -1.0)), arma::ones(1));
      arma::mat udiagdenom1 = arma::trans(arma::repmat(arma::pow(ua1.t(), -1.0), 1, 2));
      udiagdenom1.reshape(1, 2 *(ntimes - 1));
      arma::mat udiagdenom = arma::join_rows(arma::join_rows(arma::ones(1), udiagdenom1), arma::ones(1));
      udiagdenom.reshape(2, ntimes);
      gicmat.diag(0) = ldiagnum / arma::prod(ldiagdenom, 0);
      gicmat.diag(1) = - la1 % alphpow(lalpha, lag1);
      gicmat.diag(-1) = - la1 % alphpow(lalpha, lag1); 
      ggicmat.diag(0) = udiagnum / arma::prod(udiagdenom, 0);
      ggicmat.diag(1) = - ua1 % alphpow(ualpha, lag1);
      ggicmat.diag(-1) = - ua1 % alphpow(ualpha, lag1);

      }

     } else if(corrmod == "uniform") {
      
     // construct icmat
     double lphi = std::log(alpha) - std::log(1 - alpha) - h;
     double uphi = std::log(alpha) - std::log(1 - alpha) + h;
     double lalpha = std::exp(lphi) / (1 + std::exp(lphi));
     double ualpha = std::exp(uphi) / (1 + std::exp(uphi));
     arma::sp_mat oneone(1,1);
     oneone(0,0) = 1;
     arma::sp_mat onemat = arma::repmat(arma::repmat(oneone, ntimes, 1), 1, ntimes);
     arma::sp_mat cmat = alpha * onemat;
     cmat.diag(0) = - (1 + (xntimes - 2) * alpha) * arma::ones(ntimes);
     double analmult = (xntimes - 1) * (alpha - 1) * (alpha + ntimes1);
     icmat = std::pow(analmult, -1.0) * cmat;

     if(diffmeth == "analytic"){

      // construct gicmat
      arma::sp_mat gcmat = - (1 + (xntimes -1) * std::pow(alpha, 2.0)) * onemat;
      gcmat.diag(0) = (alpha * (xntimes - 1) * (2 + (xntimes - 2) * alpha)) * arma::ones(ntimes);
      gicmat = std::pow(analmult, -2.0) * gcmat;

      // construct ggicmat
      arma::sp_mat ggcmat = (2 * alpha * (xntimes - 1) * ((xntimes - 1) * std::pow(alpha, 2.0) + 3) - 2 * (xntimes - 2)) * onemat;
      ggcmat.diag(0) = -(2 * (xntimes - 1) *(std::pow(alpha, 3.0) * (xntimes - 2) * 
                         (xntimes - 1) + 3 * std::pow(alpha, 2.0) * (xntimes - 1) + 1)) * arma::ones(ntimes);
      ggicmat = std::pow(analmult, -3.0) * ggcmat;

     } else if(diffmeth == "numeric"){

     // upper and lower icmat
     double nummultl = (xntimes - 1) * (lalpha - 1) * (lalpha + ntimes1);
     double nummultu = (xntimes - 1) * (ualpha - 1) * (ualpha + ntimes1);
     arma::sp_mat lcmat = lalpha * onemat;
     arma::sp_mat ucmat = ualpha * onemat;
     lcmat.diag(0) = - (1 + (xntimes - 2) * lalpha) * arma::ones(ntimes);
     ucmat.diag(0) = - (1 + (xntimes - 2) * ualpha) * arma::ones(ntimes);
     gicmat = std::pow(nummultl, -1.0) * lcmat;
     ggicmat = std::pow(nummultu, -1.0) * ucmat;

     }

    }

   } else {

    // one time point only

     // construct icmat, gicmat and ggicmat
     icmat(0,0) = 1;
     gicmat(0,0) = 1;
     ggicmat(0,0) = 1;

   }

  if(diffmeth == "analytic"){

  // output data
  return Rcpp::List::create(Rcpp::Named("icmat")=icmat,
                            Rcpp::Named("gicmat")=gicmat,
                            Rcpp::Named("ggicmat")=ggicmat);
  } else {

  // output data
  return Rcpp::List::create(Rcpp::Named("icmat")=icmat,
                            Rcpp::Named("licmat")=gicmat,
                            Rcpp::Named("uicmat")=ggicmat);

  }

}

