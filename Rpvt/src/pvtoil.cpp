// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "pvtgas_header.h"



// ************************************** OIL PVT CORRELATIONS ****************************************

// ************ STANDING CORRELATIONS ************
// ************ STANDING CORRELATIONS ************
// ************ STANDING CORRELATIONS ************

// [[Rcpp::export]]
double PB_STANDING(double t, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  t = t - 459.67;
  double Pb;
  Pb = 18.2 * (std::pow((rsi / spgr),0.83) * std::pow(10.,(0.00091 * t - 0.0125 * api)) - 1.4);
  return(Pb);
}


// [[Rcpp::export]]
arma::vec RSi_STANDING(double t, const double pb, const double api, const double spgr) {

  // t in R
  // p in psia
  arma::vec Rsi(2);
  t = t - 459.67;
  Rsi(0) = spgr * std::pow((pb / 18.2 + 1.4) * std::pow(10.,(0.0125 * api - 0.00091 * t)),(1. / 0.83));
  Rsi(1) = spgr * (1. / 0.83) * (1. / 18.2) * std::pow((pb / 18.2 + 1.4),(1. / 0.83) - 1.) * std::pow(std::pow(10.,0.0125 * api - 0.00091 * t),1. / 0.83);
  return(Rsi);
}


// [[Rcpp::export]]
arma::vec RS_STANDING(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  double Pb;
  arma::vec Rs(2);
  Pb = PB_STANDING(t,api,spgr,rsi);
  t = t - 459.67;
  if (p >= Pb) {
    Rs(0) = rsi;
    Rs(1) = 0.;
  } else {
    Rs(0) = spgr * std::pow((p / 18.2 + 1.4) * std::pow(10.,(0.0125 * api - 0.00091 * t)),(1. / 0.83));
    Rs(1) = spgr * (1. / 0.83) * (1. / 18.2) * std::pow((p / 18.2 + 1.4),(1. / 0.83) - 1.) * std::pow(std::pow(10.,0.0125 * api - 0.00091 * t),1. / 0.83);
  }
  return(Rs);
}


// [[Rcpp::export]]
arma::vec BOB_STANDING(double t, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  t = t - 459.67;
  double spgr_o;
  arma::vec Bo(2);
  spgr_o = 141.5 / (131.5 + api);
  Bo(0) = 0.9759 + 0.00012 * std::pow((rsi * std::pow(spgr / spgr_o,0.5) + 1.25 * t),1.2);
  Bo(1) = 1.2 * 0.00012 * std::pow((rsi * std::pow(spgr / spgr_o,0.5) + 1.25 * t),0.2) * std::pow(spgr / spgr_o,0.5);
  return(Bo);
}


// [[Rcpp::export]]
double CO_UNDERSAT_STANDING(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  double spgr_o; double Pb; double Bob; double RHOob; double Co;
  spgr_o = 141.5 / (131.5 + api);
  Pb = PB_STANDING(t,api,spgr,rsi);
  Bob = BOB_STANDING(t,api,spgr,rsi)(0);
  RHOob = (62.37 * spgr_o + 0.0136 * spgr * rsi) / Bob;
  Co = 1e-6 * exp((RHOob + 0.004347 * (p - Pb) - 79.1) / (7.141e-4 * (p - Pb) - 12.938));
  return(Co);
}


// [[Rcpp::export]]
double CO_UNDERSAT_SPIVEY(double t, const double p, const double pb, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  arma::vec C0n = {3.011, -0.0835, 3.51, 0.327, -1.918, 2.52};
  arma::vec C1n = {-2.6254, -0.259, -0.0289, -0.608, -0.642, -2.73};
  arma::vec C2n = {0.497, 0.382, -0.0584, 0.0911, 0.154, 0.429};
  arma::vec Zn(6);
  arma::vec Xn(6);
  double derivative_1; double derivative_2;double Cofb;double Co;
  t = t - 459.67;
  Xn = {log(api), log(spgr), log(pb), log(p/pb), log(rsi), log(t)};
  double sum_z = 0.;
  for (int i = 0; i < 6; i++) {
    Zn(i) = C0n(i) + C1n(i) * Xn(i) + C2n(i) * Xn(i) *  Xn(i);
    sum_z = sum_z + Zn(i);
  }
  Cofb = exp(2.434 + 0.475 * sum_z + 0.048 * sum_z * sum_z);
  derivative_1 = (-0.608 + 0.1822 * log(p/pb))/p;
  derivative_2 = Cofb * (0.475 + 0.096 * sum_z) * derivative_1;
  Co = (Cofb + (p - pb) * derivative_2) * 1e-6;
  return(Co);
}

// [[Rcpp::export]]
arma::vec BO_STANDING(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  double Pb; double Bob; double Rs; double Co;
  arma::vec Bo(2);
  Pb = PB_STANDING(t,api,spgr,rsi);
  if (p >= Pb) {
    Bob = BOB_STANDING(t,api,spgr,rsi)(0);
    // Co = CO_UNDERSAT_STANDING(t,p,api,spgr,rsi);
    Co = CO_UNDERSAT_SPIVEY(t,p,Pb,api,spgr,rsi);
    Bo(0) = Bob * exp(Co * (Pb - p));
    Bo(1) = 0.;
  } else {
    Rs = RS_STANDING(t,p,api,spgr,rsi)(0);
    Bo(0) = BOB_STANDING(t,api,spgr,Rs)(0);
    Bo(1) = BOB_STANDING(t,api,spgr,Rs)(1);
  }
  return(Bo);
}




// [[Rcpp::export]]
double CO_STANDING(double t, const double p, const double api, const double spgr, const double rsi, const double tsc, const double psca, const double tpc, const double ppc) {

  // t in R
  // p in psia
  double Pb; double Bg; double Rs; double dRs_dP; double Bo; double dBo_dRs; double Co;
  Pb = PB_STANDING(t,api,spgr,rsi);
  if (p >= Pb) {
    // Co = CO_UNDERSAT_STANDING(t,p,api,spgr,rsi);
    Co = CO_UNDERSAT_SPIVEY(t,p,Pb,api,spgr,rsi);
  } else {
    Bg = B_GAS_DAK(t,p,tsc,psca,tpc,ppc)(1);
    Rs = RS_STANDING(t,p,api,spgr,rsi)(0);
    dRs_dP = RS_STANDING(t,p,api,spgr,rsi)(1);
    Bo = BOB_STANDING(t,api,spgr,Rs)(0);
    dBo_dRs = BOB_STANDING(t,api,spgr,Rs)(1);
    Co = (Bg - dBo_dRs) * dRs_dP / Bo;
  }
  return(Co);
}



// [[Rcpp::export]]
double DENSITY_STANDING(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  double spgr_o; double Bo; double Rs; double Density;
  spgr_o = 141.5 / (131.5 + api);
  Rs = RS_STANDING(t,p,api,spgr,rsi)(0);
  Bo = BO_STANDING(t,p,api,spgr,rsi)(0);
  Density = (62.37 * spgr_o + 0.0136 * spgr * Rs) / Bo;
  return(Density);
}




// ************ VASQUEZ AND BEGGS CORRELATIONS ************
// ************ VASQUEZ AND BEGGS CORRELATIONS ************
// ************ VASQUEZ AND BEGGS CORRELATIONS ************

// [[Rcpp::export]]
double PB_VASQUEZ_BEGGS(double t, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  arma::vec c(3);
  if (api <= 30) {
    c= {0.0362,1.0937,25.724};
  } else {
    c = {0.0178,1.187,23.931};
  }
  double Pb;
  Pb = std::pow((rsi / c(0) / spgr / exp(c(2) * (api / t))), (1. / c(1)));
  return(Pb);
}


// [[Rcpp::export]]
arma::vec RSi_VASQUEZ_BEGGS(double t, const double pb, const double api, const double spgr) {

  // Correlations for Fluid Physical Property Prediction", M.E. Vasquez and H.D. Beggs, JPT 968 - 70, June 1980
  // t in R
  // p in psia
  arma::vec Rsi(2);
  arma::vec c(3);
  if (api <= 30) {
    c= {0.0362,1.0937,25.724};
  } else {
    c = {0.0178,1.187,23.931};
  }
  Rsi(0) = c(0) * spgr * std::pow(pb,c(1)) * exp(c(2) * (api / t));
  Rsi(1) = c(0) * c(1) * spgr * std::pow(pb,(c(1)-1)) * exp(c(2) * (api / t));
  return(Rsi);
}


// [[Rcpp::export]]
arma::vec RS_VASQUEZ_BEGGS(double t, const double p, const double api, const double spgr, const double rsi) {

  // Correlations for Fluid Physical Property Prediction", M.E. Vasquez and H.D. Beggs, JPT 968 - 70, June 1980
  // t in R
  // p in psia
  double Pb;
  arma::vec Rs(2);
  arma::vec c(3);
  if (api <= 30) {
    c= {0.0362,1.0937,25.724};
  } else {
    c = {0.0178,1.187,23.931};
  }
  Pb = PB_VASQUEZ_BEGGS(t,api,spgr,rsi);
  if (p >= Pb) {
    Rs(0) = rsi;
    Rs(1) = 0.;
  } else {
    Rs(0) = c(0) * spgr * std::pow(p,c(1)) * exp(c(2) * (api / t));
    Rs(1) = c(0) * c(1) * spgr * std::pow(p,(c(1)-1.)) * exp(c(2) * (api / t));
  }
  return(Rs);
}


// [[Rcpp::export]]
arma::vec BOB_VASQUEZ_BEGGS(double t, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  arma::vec Bo(2);
  arma::vec c(3);
  if (api <= 30) {
    c = {4.677e-4,1.751e-5,-1.811e-8};
  } else {
    c = {4.67e-4,1.1e-5,1.377e-9};
  }
  t = t - 459.67;
  Bo(0) = 1. + c(0) * rsi + c(1) * (api / spgr) * (t - 60.) + c(2) * rsi * (api / spgr) * (t - 60.);
  Bo(1) = c(0) + c(2) * (api / spgr) * (t - 60.);
  return(Bo);
}


// [[Rcpp::export]]
double CO_UNDERSAT_VASQUEZ_BEGGS(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  double Co;
  t = t - 459.67;
  Co = 1e-5 * (-1433 + 5. * rsi + 17.2 * t - 1180 * spgr + 12.61 * api) / p;
  return(Co);
}


// [[Rcpp::export]]
arma::vec BO_VASQUEZ_BEGGS(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  double Pb; double Bob; double Rs; double Co;
  arma::vec Bo(2);
  Pb = PB_VASQUEZ_BEGGS(t,api,spgr,rsi);
  if (p >= Pb) {
    Bob = BOB_VASQUEZ_BEGGS(t,api,spgr,rsi)(0);
    // Co = CO_UNDERSAT_VASQUEZ_BEGGS(t,p,api,spgr,rsi);
    Co = CO_UNDERSAT_SPIVEY(t,p,Pb,api,spgr,rsi);
    Bo(0) = Bob * exp(Co * (Pb - p));
    Bo(1) = 0.;
  } else {
    Rs = RS_VASQUEZ_BEGGS(t,p,api,spgr,rsi)(0);
    Bo(0) = BOB_VASQUEZ_BEGGS(t,api,spgr,Rs)(0);
    Bo(1) = BOB_VASQUEZ_BEGGS(t,api,spgr,Rs)(1);
  }
  return(Bo);
}




// [[Rcpp::export]]
double CO_VASQUEZ_BEGGS(double t, const double p, const double api, const double spgr, const double rsi, const double tsc, const double psca, const double tpc, const double ppc) {

  // t in R
  // p in psia
  double Pb; double Bg; double Rs; double dRs_dP; double Bo; double dBo_dRs; double Co;
  Pb = PB_VASQUEZ_BEGGS(t,api,spgr,rsi);
  if (p >= Pb) {
    // Co = CO_UNDERSAT_VASQUEZ_BEGGS(t,p,api,spgr,rsi);
    Co = CO_UNDERSAT_SPIVEY(t,p,Pb,api,spgr,rsi);
  } else {
    Bg = B_GAS_DAK(t,p,tsc,psca,tpc,ppc)(1);
    Rs = RS_VASQUEZ_BEGGS(t,p,api,spgr,rsi)(0);
    dRs_dP = RS_VASQUEZ_BEGGS(t,p,api,spgr,rsi)(1);
    Bo = BOB_VASQUEZ_BEGGS(t,api,spgr,Rs)(0);
    dBo_dRs = BOB_VASQUEZ_BEGGS(t,api,spgr,Rs)(1);
    Co = (Bg - dBo_dRs) * dRs_dP / Bo;
  }
  return(Co);
}



// [[Rcpp::export]]
double DENSITY_VASQUEZ_BEGGS(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  double spgr_o; double Bo; double Rs; double Density;
  spgr_o = 141.5 / (131.5 + api);
  Rs = RS_VASQUEZ_BEGGS(t,p,api,spgr,rsi)(0);
  Bo = BO_VASQUEZ_BEGGS(t,p,api,spgr,rsi)(0);
  Density = (62.37 * spgr_o + 0.0136 * spgr * Rs) / Bo;
  return(Density);
}











// ************ FARSHAD AND PETROSKY CORRELATIONS ************
// ************ FARSHAD AND PETROSKY CORRELATIONS ************
// ************ FARSHAD AND PETROSKY CORRELATIONS ************

// [[Rcpp::export]]
double PB_FARSHAD_PETROSKY(double t, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  double Pb;
  t = t - 459.67;
  double x = 4.561e-5 * std::pow(t,1.3911) - 7.916e-4 * std::pow(api,1.541);
  Pb = 112.727 * (std::pow(rsi,0.5774) * std::pow(10., x) / std::pow(spgr,0.8439) - 12.34);
  return(Pb);
}


// [[Rcpp::export]]
arma::vec RSi_FARSHAD_PETROSKY(double t, const double pb, const double api, const double spgr) {

  // t in R
  // p in psia
  arma::vec Rsi(2);
  t = t - 459.67;
  double x = -4.561e-5 * std::pow(t,1.3911) + 7.916e-4 * std::pow(api,1.541);
  Rsi(0) = std::pow((pb / 112.727 + 12.34) * std::pow(spgr,0.8439) * std::pow(10.,x),1.73184);
  Rsi(1) = (1.73184/112.727) * std::pow((pb / 112.727 + 12.34),0.73184) * std::pow(std::pow(spgr,0.8439) * std::pow(10.,x),1.73184);
  return(Rsi);
}


// [[Rcpp::export]]
arma::vec RS_FARSHAD_PETROSKY(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  double Pb;
  arma::vec Rs(2);
  Pb = PB_FARSHAD_PETROSKY(t,api,spgr,rsi);
  t = t - 459.67;
  if (p >= Pb) {
    Rs(0) = rsi;
    Rs(1) = 0.;
  } else {
    double x = -4.561e-5 * std::pow(t,1.3911) + 7.916e-4 * std::pow(api,1.541);
    Rs(0) = std::pow((p / 112.727 + 12.34) * std::pow(spgr,0.8439) * std::pow(10.,x),1.73184);
    Rs(1) = (1.73184/112.727) * std::pow((p / 112.727 + 12.34),0.73184) * std::pow(std::pow(spgr,0.8439) * std::pow(10.,x),1.73184);
  }
  return(Rs);
}


// [[Rcpp::export]]
arma::vec BOB_FARSHAD_PETROSKY(double t, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  arma::vec Bo(2);
  double spgr_o = 141.5 / (131.5 + api);
  t = t - 459.67;
  Bo(0) = 1.0113 + 7.2046e-5 * std::pow((std::pow(rsi,0.3738) * (std::pow(spgr,0.2914) / std::pow(spgr_o, 0.6265)) + 0.24626 * std::pow(t,0.5371)),3.0936);
  Bo(1) = 7.2046e-5*3.0936*std::pow((std::pow(rsi,0.3738) * (std::pow(spgr,0.2914) / std::pow(spgr_o, 0.6265)) + 0.24626 * std::pow(t,0.5371)),2.0936)*0.3738*std::pow(rsi,-0.6262)*(std::pow(spgr,0.2914) / std::pow(spgr_o, 0.6265));
  return(Bo);
}


// [[Rcpp::export]]
double CO_UNDERSAT_FARSHAD_PETROSKY(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  double Co;
  t = t - 459.67;
  Co = 1.705e-7 * std::pow(rsi,0.69357) * std::pow(spgr,0.1885) * std::pow(api,0.3272) * std::pow(t,0.6729) * std::pow(p,-0.5906);
  return(Co);
}


// [[Rcpp::export]]
arma::vec BO_FARSHAD_PETROSKY(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  double Pb; double Bob; double Rs; double Co;
  arma::vec Bo(2);
  Pb = PB_FARSHAD_PETROSKY(t,api,spgr,rsi);
  if (p >= Pb) {
    Bob = BOB_FARSHAD_PETROSKY(t,api,spgr,rsi)(0);
    // Co = CO_UNDERSAT_FARSHAD_PETROSKY(t,p,api,spgr,rsi);
    Co = CO_UNDERSAT_SPIVEY(t,p,Pb,api,spgr,rsi);
    Bo(0) = Bob * exp(Co * (Pb - p));
    Bo(1) = 0.;
  } else {
    Rs = RS_FARSHAD_PETROSKY(t,p,api,spgr,rsi)(0);
    Bo(0) = BOB_FARSHAD_PETROSKY(t,api,spgr,Rs)(0);
    Bo(1) = BOB_FARSHAD_PETROSKY(t,api,spgr,Rs)(1);
  }
  return(Bo);
}



// [[Rcpp::export]]
double CO_FARSHAD_PETROSKY(double t, const double p, const double api, const double spgr, const double rsi, const double tsc, const double psca, const double tpc, const double ppc) {

  // t in R
  // p in psia
  double Pb; double Bg; double Rs; double dRs_dP; double Bo; double dBo_dRs; double Co;
  Pb = PB_FARSHAD_PETROSKY(t,api,spgr,rsi);
  if (p >= Pb) {
    // Co = CO_UNDERSAT_FARSHAD_PETROSKY(t,p,api,spgr,rsi);
    Co = CO_UNDERSAT_SPIVEY(t,p,Pb,api,spgr,rsi);
  } else {
    Bg = B_GAS_DAK(t,p,tsc,psca,tpc,ppc)(1);
    Rs = RS_FARSHAD_PETROSKY(t,p,api,spgr,rsi)(0);
    dRs_dP = RS_FARSHAD_PETROSKY(t,p,api,spgr,rsi)(1);
    Bo = BOB_FARSHAD_PETROSKY(t,api,spgr,Rs)(0);
    dBo_dRs = BOB_FARSHAD_PETROSKY(t,api,spgr,Rs)(1);
    Co = (Bg - dBo_dRs) * dRs_dP / Bo;
  }
  return(Co);
}



// [[Rcpp::export]]
double DENSITY_FARSHAD_PETROSKY(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  double spgr_o; double Bo; double Rs; double Density;
  spgr_o = 141.5 / (131.5 + api);
  Rs = RS_FARSHAD_PETROSKY(t,p,api,spgr,rsi)(0);
  Bo = BO_FARSHAD_PETROSKY(t,p,api,spgr,rsi)(0);
  Density = (62.37 * spgr_o + 0.0136 * spgr * Rs) / Bo;
  return(Density);
}













  // ************ AL MARHOUN CORRELATIONS ************
  // ************ AL MARHOUN CORRELATIONS ************
  // ************ AL MARHOUN CORRELATIONS ************

  // [[Rcpp::export]]
  double PB_AL_MARHOUN(double t, const double api, const double spgr, const double rsi) {

    // t in R
    // p in psia
    double Pb;
    double spgr_o = 141.5 / (131.5 + api);
    arma::vec a = {0.00538088,0.715082,-1.87784,3.1437,1.32657};
    Pb = a(0) * std::pow(rsi,a(1)) * std::pow(spgr,a(2)) * std::pow(spgr_o,a(3)) * std::pow(t,a(4));
    return(Pb);
  }


// [[Rcpp::export]]
arma::vec RSi_AL_MARHOUN(double t, const double pb, const double api, const double spgr) {

  // t in R
  // p in psia
  arma::vec Rsi(2);
  arma::vec a = {1490.28,2.62605,1.398441,-4.396279,-1.85513};
  double spgr_o = 141.5 / (131.5 + api);
  Rsi(0) = a(0) * std::pow(spgr,a(1)) * std::pow(pb,a(2)) * std::pow(spgr_o,a(3)) * std::pow(t,a(4));
  Rsi(1) = a(2) * a(0) * std::pow(spgr,a(1)) * std::pow(pb,(a(2)-1.)) * std::pow(spgr_o,a(3)) * std::pow(t,a(4));
  return(Rsi);
}


// [[Rcpp::export]]
arma::vec RS_AL_MARHOUN(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  double Pb;
  arma::vec Rs(2);
  arma::vec a = {1490.28,2.62605,1.398441,-4.396279,-1.85513};
  double spgr_o = 141.5 / (131.5 + api);
  Pb = PB_AL_MARHOUN(t,api,spgr,rsi);
  if (p >= Pb) {
    Rs(0) = rsi;
    Rs(1) = 0.;
  } else {
    Rs(0) = a(0) * std::pow(spgr,a(1)) * std::pow(p,a(2)) * std::pow(spgr_o,a(3)) * std::pow(t,a(4));
    Rs(1) = a(2) * a(0) * std::pow(spgr,a(1)) * std::pow(p,(a(2)-1.)) * std::pow(spgr_o,a(3)) * std::pow(t,a(4));
  }
  return(Rs);
}


// [[Rcpp::export]]
arma::vec BOB_AL_MARHOUN(double t, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  arma::vec Bo(2);
  double spgr_o = 141.5 / (131.5 + api);
  double f = std::pow(rsi,0.74239) * std::pow(spgr,0.323294) * std::pow(spgr_o,-1.20204);
  double df_drs = 0.74239 * std::pow(rsi,-0.25761) * std::pow(spgr,0.323294) * std::pow(spgr_o,-1.20204);
  Bo(0) = 0.497069 + 0.862963e-3 * t + 0.182594e-2 * f + 0.318099e-5 * f * f;
  Bo(1) = 0.182594e-2 * df_drs + 2 *  0.318099e-5 * f * df_drs;
  return(Bo);
}


// [[Rcpp::export]]
double CO_UNDERSAT_AL_MARHOUN(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  double Co;
  arma::vec a = {-14.1042,2.7314,-0.0000560605,-580.8778};
  double Pb = PB_AL_MARHOUN(t,api,spgr,rsi);
  double Bob = BOB_AL_MARHOUN(t,api,spgr,rsi)(0);
  double spgr_o = 141.5 / (131.5 + api);
  double spgr_ob = (spgr_o + 0.000218 * rsi * spgr) / Bob;
  Co = exp(a(0) + a(1) / spgr_ob + a(2) * (p - Pb) / std::pow(spgr_ob,3) + a(3) / t);
  return(Co);
}


// [[Rcpp::export]]
arma::vec BO_AL_MARHOUN(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  double Pb; double Bob; double Rs; double Co;
  arma::vec Bo(2);
  Pb = PB_AL_MARHOUN(t,api,spgr,rsi);
  if (p >= Pb) {
    Bob = BOB_AL_MARHOUN(t,api,spgr,rsi)(0);
    // Co = CO_UNDERSAT_AL_MARHOUN(t,p,api,spgr,rsi);
    Co = CO_UNDERSAT_SPIVEY(t,p,Pb,api,spgr,rsi);
    Bo(0) = Bob * exp(Co * (Pb - p));
    Bo(1) = 0.;
  } else {
    Rs = RS_AL_MARHOUN(t,p,api,spgr,rsi)(0);
    Bo(0) = BOB_AL_MARHOUN(t,api,spgr,Rs)(0);
    Bo(1) = BOB_AL_MARHOUN(t,api,spgr,Rs)(1);
  }
  return(Bo);
}



// [[Rcpp::export]]
double CO_AL_MARHOUN(double t, const double p, const double api, const double spgr, const double rsi, const double tsc, const double psca, const double tpc, const double ppc) {

  // t in R
  // p in psia
  double Pb; double Bg; double Rs; double dRs_dP; double Bo; double dBo_dRs; double Co;
  Pb = PB_AL_MARHOUN(t,api,spgr,rsi);
  if (p >= Pb) {
    // Co = CO_UNDERSAT_AL_MARHOUN(t,p,api,spgr,rsi);
    Co = CO_UNDERSAT_SPIVEY(t,p,Pb,api,spgr,rsi);
  } else {
    Bg = B_GAS_DAK(t,p,tsc,psca,tpc,ppc)(1);
    Rs = RS_AL_MARHOUN(t,p,api,spgr,rsi)(0);
    dRs_dP = RS_AL_MARHOUN(t,p,api,spgr,rsi)(1);
    Bo = BOB_AL_MARHOUN(t,api,spgr,Rs)(0);
    dBo_dRs = BOB_AL_MARHOUN(t,api,spgr,Rs)(1);
    Co = (Bg - dBo_dRs) * dRs_dP / Bo;
  }
  return(Co);
}



// [[Rcpp::export]]
double DENSITY_AL_MARHOUN(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  double spgr_o; double Bo; double Rs; double Density;
  spgr_o = 141.5 / (131.5 + api);
  Rs = RS_AL_MARHOUN(t,p,api,spgr,rsi)(0);
  Bo = BO_AL_MARHOUN(t,p,api,spgr,rsi)(0);
  Density = (62.37 * spgr_o + 0.0136 * spgr * Rs) / Bo;
  return(Density);
}








// ************ GLASO CORRELATIONS ************
// ************ GLASO CORRELATIONS ************
// ************ GLASO CORRELATIONS ************

// [[Rcpp::export]]
double PB_GLASO(double t, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  double Pb;
  t = t - 459.67;
  double x = std::pow((rsi/spgr),0.816) * std::pow(t,0.172) / std::pow(api,0.989);
  double a = -0.30218 * log10(x) * log10(x) + 1.7447 * log10(x) + 1.7669;
  Pb = std::pow(10.,a);
  return(Pb);
}


// [[Rcpp::export]]
arma::vec RSi_GLASO(double t, const double pb, const double api, const double spgr) {

  // t in R
  // p in psia
  t = t - 459.67;
  double x; double z; double a; double b; double c;
  arma::vec Rsi(2);
  a = -0.30218;
  b = 1.7447;
  c = -log10(pb) + 1.7669;
  double u; double dc_dp; double du_dc; double dx_du;
  double delta = b * b - 4 * a * c;
  double y1 = (-b+sqrt(delta))/(2*a);
  double y2 = (-b-sqrt(delta))/(2*a);
  x = std::min(y1,y2);
  z = std::pow(10.,x);
  double sol = z * std::pow(api,0.989)/std::pow(t,0.172);
  Rsi(0) = std::pow(sol,(1./0.816)) * spgr;
  if (y1 <= y2) {
    u = y1;
    dc_dp = -1./pb/log(10.);
    du_dc = -1. * std::pow(b*b-4.*a*c,-0.5);
    dx_du = log(10.)*std::pow(10.,u);
  } else {
    u = y2;
    dc_dp = -1./pb/log(10.);
    du_dc = 1. * std::pow(b*b-4.*a*c,-0.5);
    dx_du = log(10.)*std::pow(10.,u);
  }
  double dRsi_dx = (1./0.816) * std::pow(z,(1./0.816)-1.) * std::pow(std::pow(api,0.989)/std::pow(t,0.172),(1./0.816)) * spgr;
  Rsi(1) = dRsi_dx * dx_du * du_dc * dc_dp;
  return(Rsi);
}


// [[Rcpp::export]]
arma::vec RS_GLASO(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  double Pb = PB_GLASO(t,api,spgr,rsi);
  t = t - 459.67;
  double x; double z; double a; double b; double c;
  arma::vec Rs(2);
  a =  -0.30218;
  b = 1.7447;
  c = -log10(p) + 1.7669;
  double u; double dc_dp; double du_dc; double dx_du;
  double delta = b * b - 4 * a * c;
  double y1 = (-b+sqrt(delta))/(2*a);
  double y2 = (-b-sqrt(delta))/(2*a);
  x = std::min(y1,y2);
  if (p >= Pb) {
    Rs(0) = rsi;
    Rs(1) = 0.;
  } else {
    z = std::pow(10.,x);
    double sol = z * std::pow(api,0.989)/std::pow(t,0.172);
    Rs(0) = std::pow(sol,(1./0.816)) * spgr;
    if (y1 <= y2) {
      u = y1;
      dc_dp = -1./p/log(10.);
      du_dc = -1. * std::pow(b*b-4.*a*c,-0.5);
      dx_du = log(10.)*std::pow(10.,u);

    } else {
      u = y2;
      dc_dp = -1./p/log(10.);
      du_dc = 1. * std::pow(b*b-4.*a*c,-0.5);
      dx_du = log(10.)*std::pow(10.,u);
    }
    double dRs_dx = (1./0.816) * std::pow(z,(1./0.816)-1.) * std::pow(std::pow(api,0.989)/std::pow(t,0.172),(1./0.816)) * spgr;
    Rs(1) = dRs_dx * dx_du * du_dc * dc_dp;
  }
  return(Rs);
}


// [[Rcpp::export]]
arma::vec BOB_GLASO(double t, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  t = t - 459.67;
  arma::vec Bo(2);
  double spgr_o = 141.5 / (131.5 + api);
  double y = rsi * std::pow((spgr/spgr_o),0.526) + 0.968 * t;
  double log10_y = log10(y);
  double dy_dRs = std::pow(spgr/spgr_o,0.526);
  double f = -6.58511+2.91329*log10_y-0.27683*log10_y*log10_y;
  double df_dlog10_y = 2.91329-2*0.27683*log10_y;
  double dlog10_y_dy = 1./y/log(10.);
  double dB_df = log(10.) * std::pow(10.,f) ;
  Bo(0) = std::pow(10.,f) + 1.;
  Bo(1) = dB_df * df_dlog10_y * dlog10_y_dy * dy_dRs;
  return(Bo);
}


// [[Rcpp::export]]
double CO_UNDERSAT_GLASO(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia
  double Co;
  arma::vec a = {-14.1042,2.7314,-0.0000560605,-580.8778};
  double Pb = PB_GLASO(t,api,spgr,rsi);
  double Bob = BOB_GLASO(t,api,spgr,rsi)(0);
  double spgr_o = 141.5 / (131.5 + api);
  double spgr_ob = (spgr_o + 0.000218 * rsi * spgr) / Bob;
  Co = exp(a(0) + a(1) / spgr_ob + a(2) * (p - Pb) / std::pow(spgr_ob,3.) + a(3) / t);
  return(Co);
}


// [[Rcpp::export]]
arma::vec BO_GLASO(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  double Pb; double Bob; double Rs; double Co;
  arma::vec Bo(2);
  Pb = PB_GLASO(t,api,spgr,rsi);
  if (p >= Pb) {
    Bob = BOB_GLASO(t,api,spgr,rsi)(0);
    // Co = CO_UNDERSAT_GLASO(t,p,api,spgr,rsi);
    Co = CO_UNDERSAT_SPIVEY(t,p,Pb,api,spgr,rsi);
    Bo(0) = Bob * exp(Co * (Pb - p));
    Bo(1) = 0.;
  } else {
    Rs = RS_GLASO(t,p,api,spgr,rsi)(0);
    Bo(0) = BOB_GLASO(t,api,spgr,Rs)(0);
    Bo(1) = BOB_GLASO(t,api,spgr,Rs)(1);
  }
  return(Bo);
}



// [[Rcpp::export]]
double CO_GLASO(double t, const double p, const double api, const double spgr, const double rsi, const double tsc, const double psca, const double tpc, const double ppc) {

  // t in R
  // p in psia
  double Pb; double Bg; double Rs; double dRs_dP; double Bo; double dBo_dRs; double Co;
  Pb = PB_GLASO(t,api,spgr,rsi);
  if (p >= Pb) {
    // Co = CO_UNDERSAT_GLASO(t,p,api,spgr,rsi);
    Co = CO_UNDERSAT_SPIVEY(t,p,Pb,api,spgr,rsi);
  } else {
    Bg = B_GAS_DAK(t,p,tsc,psca,tpc,ppc)(1);
    Rs = RS_GLASO(t,p,api,spgr,rsi)(0);
    dRs_dP = RS_GLASO(t,p,api,spgr,rsi)(1);
    Bo = BOB_GLASO(t,api,spgr,Rs)(0);
    dBo_dRs = BOB_GLASO(t,api,spgr,Rs)(1);
    Co = (Bg - dBo_dRs) * dRs_dP / Bo;
  }
  return(Co);
}



// [[Rcpp::export]]
double DENSITY_GLASO(double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  double spgr_o; double Bo; double Rs; double Density;
  spgr_o = 141.5 / (131.5 + api);
  Rs = RS_GLASO(t,p,api,spgr,rsi)(0);
  Bo = BO_GLASO(t,p,api,spgr,rsi)(0);
  Density = (62.37 * spgr_o + 0.0136 * spgr * Rs) / Bo;
  return(Density);
}










// ************ VISCOSITY CORRELATIONS ************
// ************ VISCOSITY CORRELATIONS ************
// ************ VISCOSITY CORRELATIONS ************


// [[Rcpp::export]]
double MUO_BEGGS_ROBINSON(const std::string pvt_model, double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  double Pb = 0.0; double Rs = 0.0; double z; double y; double x;
  double muo_dead; double muo_live; double A; double B;
  double MUO_BEGGS_ROBINSON;
  if (pvt_model == "Standing") {
    Pb = PB_STANDING(t,api,spgr,rsi);
    Rs = RS_STANDING(t,p,api,spgr,rsi)(0);
  }
  if (pvt_model == "Vasquez_Beggs") {
    Pb = PB_VASQUEZ_BEGGS(t,api,spgr,rsi);
    Rs = RS_VASQUEZ_BEGGS(t,p,api,spgr,rsi)(0);
  }
  if (pvt_model == "Farshad_Petrosky") {
    Pb = PB_FARSHAD_PETROSKY(t,api,spgr,rsi);
    Rs = RS_FARSHAD_PETROSKY(t,p,api,spgr,rsi)(0);
  }
  if (pvt_model == "Al_Marhoun") {
    Pb = PB_AL_MARHOUN(t,api,spgr,rsi);
    Rs = RS_AL_MARHOUN(t,p,api,spgr,rsi)(0);
  }
  if (pvt_model == "Glaso") {
    Pb = PB_GLASO(t,api,spgr,rsi);
    Rs = RS_GLASO(t,p,api,spgr,rsi)(0);
  }
  t = t - 459.67;
  if (p >= Pb) {
    z = 3.0324 - 0.02023 * api;
    y = std::pow(10.,z);
    x = y * std::pow(t,-1.163);
    muo_dead = std::pow(10.,x) - 1.;
    A = 10.715 * std::pow(rsi + 100, -0.515);
    B = 5.44 * std::pow(rsi + 150., -0.338);
    muo_live = A * std::pow(muo_dead,B);
    double C1 = 2.6;
    double C2 = 1.187;
    double C3 = -11.513;
    double C4 = -8.98 * std::pow(10.,-5.);
    double W = C1 * std::pow(p,C2) * exp(C3 + C4 * p);
    MUO_BEGGS_ROBINSON = muo_live * std::pow(p / Pb,W);
  } else {
    z = 3.0324 - 0.02023 * api;
    y = std::pow(10.,z);
    x = y * std::pow(t,-1.163);
    muo_dead = std::pow(10.,x) - 1.;
    A = 10.715 * std::pow(Rs + 100., -0.515);
    B = 5.44 * std::pow(Rs + 150., -0.338);
    MUO_BEGGS_ROBINSON = A * std::pow(muo_dead,B);
  }
  return(MUO_BEGGS_ROBINSON);
}



// [[Rcpp::export]]
double MUO_AL_MARHOUN(const std::string pvt_model, double t, const double p, const double api, const double spgr, const double rsi) {

  // t in R
  // p in psia

  double Pb = 0.0; double Rs = 0.0; double Bo; double Bob;
  double MUO_AL_MARHOUN; double spgr_ob = 0.0; double spgr_ub = 0.0;
  double a_ub; double b_ub; double visc_ub;
  arma::vec a = {54.56805426,-7.179530398,-36.447,4.478878992};
  arma::vec b = {10.715,100.0,-0.515,5.44,150.0,-0.338};
  double spgr_o = 141.5 / (131.5 + api);

  double visc_dead = exp(a(0) + a(1) * log(t-459.67) + a(2) * log(log(api)) + a(3) * log(t-459.67) * log(log(api)));
  double a_ob = b(0) * std::pow(rsi + b(1),b(2));
  double b_ob = b(3) * std::pow(rsi + b(4),b(5));
  double visc_ob = a_ob * std::pow(visc_dead,b_ob);

  if (pvt_model == "Standing") {
    Pb = PB_STANDING(t,api,spgr,rsi);
    Rs = RS_STANDING(t,p,api,spgr,rsi)(0);
    Bo = BOB_STANDING(t,api,spgr,Rs)(0);
    spgr_ub = (spgr_o + 0.000218 * Rs * spgr) / Bo;
    Bob = BOB_STANDING(t,api,spgr,rsi)(0);
    spgr_ob = (spgr_o + 0.000218 * rsi * spgr) / Bob;
  }
  if (pvt_model == "Vasquez_Beggs") {
    Pb = PB_VASQUEZ_BEGGS(t,api,spgr,rsi);
    Rs = RS_VASQUEZ_BEGGS(t,p,api,spgr,rsi)(0);
    Bo = BOB_VASQUEZ_BEGGS(t,api,spgr,Rs)(0);
    spgr_ub = (spgr_o + 0.000218 * Rs * spgr) / Bo;
    Bob = BOB_VASQUEZ_BEGGS(t,api,spgr,rsi)(0);
    spgr_ob = (spgr_o + 0.000218 * rsi * spgr) / Bob;
  }
  if (pvt_model == "Farshad_Petrosky") {
    Pb = PB_FARSHAD_PETROSKY(t,api,spgr,rsi);
    Rs = RS_FARSHAD_PETROSKY(t,p,api,spgr,rsi)(0);
    Bo = BOB_FARSHAD_PETROSKY(t,api,spgr,Rs)(0);
    spgr_ub = (spgr_o + 0.000218 * Rs * spgr) / Bo;
    Bob = BOB_FARSHAD_PETROSKY(t,api,spgr,rsi)(0);
    spgr_ob = (spgr_o + 0.000218 * rsi * spgr) / Bob;
  }
  if (pvt_model == "Al_Marhoun") {
    Pb = PB_AL_MARHOUN(t,api,spgr,rsi);
    Rs = RS_AL_MARHOUN(t,p,api,spgr,rsi)(0);
    Bo = BOB_AL_MARHOUN(t,api,spgr,Rs)(0);
    spgr_ub = (spgr_o + 0.000218 * Rs * spgr) / Bo;
    Bob = BOB_AL_MARHOUN(t,api,spgr,rsi)(0);
    spgr_ob = (spgr_o + 0.000218 * rsi * spgr) / Bob;
  }
  if (pvt_model == "Glaso") {
    Pb = PB_GLASO(t,api,spgr,rsi);
    Rs = RS_GLASO(t,p,api,spgr,rsi)(0);
    Bo = BOB_GLASO(t,api,spgr,Rs)(0);
    spgr_ub = (spgr_o + 0.000218 * Rs * spgr) / Bo;
    Bob = BOB_GLASO(t,api,spgr,rsi)(0);
    spgr_ob = (spgr_o + 0.000218 * rsi * spgr) / Bob;
  }
  if (p >= Pb) {
    MUO_AL_MARHOUN = exp(log(visc_ob) + 0.000151292 * std::pow(spgr_ob,2.) * (p - Pb));
  } else {
    a_ub = b(0) * std::pow(Rs + b(1),b(2));
    b_ub = b(3) * std::pow(Rs + b(4),b(5));
    visc_ub = a_ub * std::pow(visc_dead,b_ub);
    MUO_AL_MARHOUN = exp(log(visc_ub) + 0.000151292 * std::pow(spgr_ub,2.) * (p - Pb));
  }
  return(MUO_AL_MARHOUN);
}


// ************ PVT TABLE GENERATION ************
// ************ PVT TABLE GENERATION ************
// ************ PVT TABLE GENERATION ************


// [[Rcpp::export]]
arma::mat PVT_OIL_PROPERTIES_STANDING(double t, arma::vec p, double spgr, double api, double rsi, double tpc, double ppc) {

  double psca = 14.696;              // Psia
  double tsc = (60. + 459.67);        // R
  int lp = p.size();
  arma::mat results_table(lp,4);
  for (int i = 0; i < lp; i++) {
    results_table(i,0) = RS_STANDING(t,p(i),api,spgr,rsi)(0);
    results_table(i,1) = BO_STANDING(t,p(i),api,spgr,rsi)(0);
    results_table(i,2) = DENSITY_STANDING(t,p(i),api,spgr,rsi);
    results_table(i,3) = CO_STANDING(t,p(i),api,spgr,rsi,tsc,psca,tpc,ppc);
  }
  return(results_table);
}




// [[Rcpp::export]]
arma::mat PVT_OIL_PROPERTIES_VASQUEZ_BEGGS(double t, arma::vec p, double spgr, double api, double rsi, double tpc, double ppc) {

  double psca = 14.696;              // Psia
  double tsc = (60. + 459.67);        // R
  int lp = p.size();
  arma::mat results_table(lp,4);
  for (int i = 0; i < lp; i++) {
    results_table(i,0) = RS_VASQUEZ_BEGGS(t,p(i),api,spgr,rsi)(0);
    results_table(i,1) = BO_VASQUEZ_BEGGS(t,p(i),api,spgr,rsi)(0);
    results_table(i,2) = DENSITY_VASQUEZ_BEGGS(t,p(i),api,spgr,rsi);
    results_table(i,3) = CO_VASQUEZ_BEGGS(t,p(i),api,spgr,rsi,tsc,psca,tpc,ppc);
  }
  return(results_table);
}


// [[Rcpp::export]]
arma::mat PVT_OIL_PROPERTIES_FARSHAD_PETROSKY(double t, arma::vec p, double spgr, double api, double rsi, double tpc, double ppc) {

  double psca = 14.696;              // Psia
  double tsc = (60. + 459.67);        // R
  int lp = p.size();
  arma::mat results_table(lp,4);
  for (int i = 0; i < lp; i++) {
    results_table(i,0) = RS_FARSHAD_PETROSKY(t,p(i),api,spgr,rsi)(0);
    results_table(i,1) = BO_FARSHAD_PETROSKY(t,p(i),api,spgr,rsi)(0);
    results_table(i,2) = DENSITY_FARSHAD_PETROSKY(t,p(i),api,spgr,rsi);
    results_table(i,3) = CO_FARSHAD_PETROSKY(t,p(i),api,spgr,rsi,tsc,psca,tpc,ppc);
  }
  return(results_table);
}


// [[Rcpp::export]]
arma::mat PVT_OIL_PROPERTIES_AL_MARHOUN(double t, arma::vec p, double spgr, double api, double rsi, double tpc, double ppc) {

  double psca = 14.696;              // Psia
  double tsc = (60. + 459.67);        // R
  int lp = p.size();
  arma::mat results_table(lp,4);
  for (int i = 0; i < lp; i++) {
    results_table(i,0) = RS_AL_MARHOUN(t,p(i),api,spgr,rsi)(0);
    results_table(i,1) = BO_AL_MARHOUN(t,p(i),api,spgr,rsi)(0);
    results_table(i,2) = DENSITY_AL_MARHOUN(t,p(i),api,spgr,rsi);
    results_table(i,3) = CO_AL_MARHOUN(t,p(i),api,spgr,rsi,tsc,psca,tpc,ppc);
  }
  return(results_table);
}


// [[Rcpp::export]]
arma::mat PVT_OIL_PROPERTIES_GLASO(double t, arma::vec p, double spgr, double api, double rsi, double tpc, double ppc) {

  double psca = 14.696;              // Psia
  double tsc = (60. + 459.67);        // R
  int lp = p.size();
  arma::mat results_table(lp,4);
  for (int i = 0; i < lp; i++) {
    results_table(i,0) = RS_GLASO(t,p(i),api,spgr,rsi)(0);
    results_table(i,1) = BO_GLASO(t,p(i),api,spgr,rsi)(0);
    results_table(i,2) = DENSITY_GLASO(t,p(i),api,spgr,rsi);
    results_table(i,3) = CO_GLASO(t,p(i),api,spgr,rsi,tsc,psca,tpc,ppc);
  }
  return(results_table);
}


// [[Rcpp::export]]
arma::vec VISC_OIL_PROPERTIES_BEGGS_ROBINSON(const std::string pvt_model, double t, arma::vec p, double spgr, double api, double rsi) {

  int lp = p.size();
  arma::vec results_table(lp);
  for (int i = 0; i < lp; i++) {
    results_table(i) = MUO_BEGGS_ROBINSON(pvt_model,t,p(i),api,spgr,rsi);
  }
  return(results_table);
}



// [[Rcpp::export]]
arma::vec VISC_OIL_PROPERTIES_AL_MARHOUN(const std::string pvt_model, double t, arma::vec p, double spgr, double api, double rsi) {

  int lp = p.size();
  arma::vec results_table(lp);
  for (int i = 0; i < lp; i++) {
    results_table(i) = MUO_AL_MARHOUN(pvt_model,t,p(i),api,spgr,rsi);
  }
  return(results_table);
}





















