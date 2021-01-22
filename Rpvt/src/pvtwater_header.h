// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us


#ifndef PVTWATER_HEADER_H
#define PVTWATER_HEADER_H

#include "RcppArmadillo.h"

arma::vec B_WATER_MCCAIN(double t, double p);
double DENSITY_WATER_MCCAIN(double t, double p, double S);
arma::vec Rs_WATER_MCCAIN(const std::string gas_saturated, double t, double p, double S);
double COMPRESSIBILITY_WATER_MCCAIN(const std::string gas_saturated, double t, double p, double S, double Bg);
double VISCOSITY_WATER_MCCAIN(double t, double p, double S);

double B_WATER_MEEHAN(const std::string gas_saturated, double t, double p, double S);
double DENSITY_WATER_MEEHAN(const std::string gas_saturated, double t, double p, double S);
arma::vec Rs_WATER_MEEHAN(const std::string gas_saturated, double t, double p, double S);
double COMPRESSIBILITY_WATER_MEEHAN(const std::string gas_saturated, double t, double p, double S, double Bg);
double VISCOSITY_WATER_MEEHAN(double t, double p, double S);

arma::vec WATER_SPIVEY(const std::string gas_saturated, double t, double p, double S, double z_factor);
double VISCOSITY_WATER_SPIVEY(double t, double p, double S);

arma::mat PVT_WATER_PROPERTIES_MCCAIN(const std::string gas_saturated, double t, arma::vec p, double S, double tpc, double ppc);
arma::vec VISC_WATER_PROPERTIES_MCCAIN(double t, arma::vec p, double S);
arma::mat PVT_WATER_PROPERTIES_MEEHAN(const std::string gas_saturated, double t, arma::vec p, double S, double tpc, double ppc);
arma::vec VISC_WATER_PROPERTIES_MEEHAN(double t, arma::vec p, double S);
arma::mat PVT_WATER_PROPERTIES_SPIVEY(const std::string gas_saturated, double t, arma::vec p, double S, double tpc, double ppc);
arma::vec VISC_WATER_PROPERTIES_SPIVEY(double t, arma::vec p, double S);

#endif
