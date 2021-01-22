// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us


#ifndef PVTGAS_HEADER_H
#define PVTGAS_HEADER_H

#include "RcppArmadillo.h"

double PPC_SUTTON(const double spgr, const arma::mat & nhc_properties, const std::string fluid);
double TPC_SUTTON(const double spgr, const arma::mat & nhc_properties, const std::string fluid);
double Z_FACTOR_DAK(double t, double p, double tpc, double ppc);
double DENSITY_GAS_DAK(double t, double p, double tpc, double ppc, double spgr);
arma::vec B_GAS_DAK(double t, double p, double tsc, double psca, double tpc, double ppc);
double COMPRESSIBILITY_GAS_DAK(double t, double p, double tpc, double ppc);
double VISCOSITY_GAS_SUTTON(double t, double p, double tpc, double ppc, double spgr);
double PSEUDO_PRESSURE_GAS(double t, double p, double tpc, double ppc, double spgr);
arma::mat PVT_GAS_PROPERTIES_DAK_SUTTON(double t, arma::vec p, const double spgr, const double tpc, const double ppc);

#endif
