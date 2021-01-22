// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us


#ifndef PVTOIL_HEADER_H
#define PVTOIL_HEADER_H

#include "RcppArmadillo.h"

double PB_STANDING(double t, const double api, const double spgr, const double rsi);
arma::vec RSi_STANDING(double t, const double pb, const double api, const double spgr);
arma::vec RS_STANDING(double t, const double p, const double api, const double spgr, const double rsi);
arma::vec BOB_STANDING(double t, const double api, const double spgr, const double rsi);
double CO_UNDERSAT_SPIVEY(double t, const double p, const double api, const double spgr, const double rsi);
arma::vec BO_STANDING(double t, const double p, const double api, const double spgr, const double rsi);
double CO_STANDING(double t, const double p, const double api, const double spgr, const double rsi, const double tsc, const double psca, const double tpc, const double ppc);
double DENSITY_STANDING(double t, const double p, const double api, const double spgr, const double rsi);

double PB_VASQUEZ_BEGGS(double t, const double api, const double spgr, const double rsi);
arma::vec RSi_VASQUEZ_BEGGS(double t, const double pb, const double api, const double spgr);
arma::vec RS_VASQUEZ_BEGGS(double t, const double p, const double api, const double spgr, const double rsi);
arma::vec BOB_VASQUEZ_BEGGS(double t, const double api, const double spgr, const double rsi);
arma::vec BO_VASQUEZ_BEGGS(double t, const double p, const double api, const double spgr, const double rsi);
double CO_VASQUEZ_BEGGS(double t, const double p, const double api, const double spgr, const double rsi, const double tsc, const double psca, const double tpc, const double ppc);
double DENSITY_VASQUEZ_BEGGS(double t, const double p, const double api, const double spgr, const double rsi);

double PB_FARSHAD_PETROSKY(double t, const double api, const double spgr, const double rsi);
arma::vec RSi_FARSHAD_PETROSKY(double t, const double pb, const double api, const double spgr);
arma::vec RS_FARSHAD_PETROSKY(double t, const double p, const double api, const double spgr, const double rsi);
arma::vec BOB_FARSHAD_PETROSKY(double t, const double api, const double spgr, const double rsi);
arma::vec BO_FARSHAD_PETROSKY(double t, const double p, const double api, const double spgr, const double rsi);
double CO_FARSHAD_PETROSKY(double t, const double p, const double api, const double spgr, const double rsi, const double tsc, const double psca, const double tpc, const double ppc);
double DENSITY_FARSHAD_PETROSKY(double t, const double p, const double api, const double spgr, const double rsi);

double PB_AL_MARHOUN(double t, const double api, const double spgr, const double rsi);
arma::vec RSi_AL_MARHOUN(double t, const double pb, const double api, const double spgr);
arma::vec RS_AL_MARHOUN(double t, const double p, const double api, const double spgr, const double rsi);
arma::vec BOB_AL_MARHOUN(double t, const double api, const double spgr, const double rsi);
arma::vec BO_AL_MARHOUN(double t, const double p, const double api, const double spgr, const double rsi);
double CO_AL_MARHOUN(double t, const double p, const double api, const double spgr, const double rsi, const double tsc, const double psca, const double tpc, const double ppc);
double DENSITY_AL_MARHOUN(double t, const double p, const double api, const double spgr, const double rsi);

double PB_GLASO(double t, const double api, const double spgr, const double rsi);
arma::vec RSi_GLASO(double t, const double pb, const double api, const double spgr);
arma::vec RS_GLASO(double t, const double p, const double api, const double spgr, const double rsi);
arma::vec BOB_GLASO(double t, const double api, const double spgr, const double rsi);
arma::vec BO_GLASO(double t, const double p, const double api, const double spgr, const double rsi);
double CO_GLASO(double t, const double p, const double api, const double spgr, const double rsi, const double tsc, const double psca, const double tpc, const double ppc);
double DENSITY_GLASO(double t, const double p, const double api, const double spgr, const double rsi);

double MUO_BEGGS_ROBINSON(std::string pvt_model, double t, const double p, const double api, const double spgr, const double rsi);
double MUO_AL_MARHOUN(std::string pvt_model, double t, const double p, const double api, const double spgr, const double rsi);


arma::mat PVT_OIL_PROPERTIES_STANDING(double t, arma::vec p, const double spgr, const double api, const double rsi, const double tpc, const double ppc);
arma::mat PVT_OIL_PROPERTIES_VASQUEZ_BEGGS(double t, arma::vec p, const double spgr, const double api, const double rsi, const double tpc, const double ppc);
arma::mat PVT_OIL_PROPERTIES_FARSHAD_PETROSKY(double t, arma::vec p, const double spgr, const double api, const double rsi, const double tpc, const double ppc);
arma::mat PVT_OIL_PROPERTIES_AL_MARHOUN(double t, arma::vec p, const double spgr, const double api, const double rsi, const double tpc, const double ppc);
arma::mat PVT_OIL_PROPERTIES_GLASO(double t, arma::vec p, const double spgr, const double api, const double rsi, const double tpc, const double ppc);


arma::vec VISC_OIL_PROPERTIES_BEGGS_ROBINSON(const std::string pvt_model, double t, arma::vec p, const double spgr, const double api, const double rsi);
arma::vec VISC_OIL_PROPERTIES_AL_MARHOUN(const std::string pvt_model, double t, arma::vec p, const double spgr, const double api, const double rsi);



#endif
