// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// ************************************** GAS PVT CORRELATIONS ****************************************


double PPC_SUTTON(const double spgr, const arma::mat & nhc_properties, const std::string fluid) {

    int s = nhc_properties.n_rows;
    arma::vec composition(s);
    arma::vec mw(s);
    arma::vec tcrit(s);
    arma::vec pcrit(s);
    const double mw_air = 28.964;
    double y_HC; double spgr_HC;
    double tcrit_HC, pcrit_HC;
    double tpc_pseudo, ppc_pseudo, epsilon, PPC_SUTTON;

    composition = nhc_properties.col(0);
    mw =  nhc_properties.col(1);
    tcrit =  nhc_properties.col(2);
    pcrit =  nhc_properties.col(3);
    y_HC = 1 - sum(composition);
    spgr_HC = (spgr - sum(composition.t() * mw) / mw_air) / y_HC;
    if (fluid == "dry_gas") {
        tcrit_HC = 120.1 + 429 * spgr_HC  - 62.9 * std::pow(spgr_HC, 2.);
        pcrit_HC = 671.1 - 14 * spgr_HC  - 34.3 * std::pow(spgr_HC, 2.);
    }  else {
        // "Gas Condensate, wet gas, rich gas
        tcrit_HC = 164.3 + 357.7 * spgr_HC  - 67.7 *  std::pow(spgr_HC, 2.);
        pcrit_HC = 744 - 125.4 * spgr_HC  + 5.9 *  std::pow(spgr_HC, 2.);
    }
    tpc_pseudo = y_HC * tcrit_HC + sum(composition.t() * tcrit);
    ppc_pseudo = y_HC * pcrit_HC + sum(composition.t() * pcrit);
    epsilon = 120. * (std::pow(sum(composition.rows(1,2)), 0.9) - std::pow(sum(composition.rows(1,2)), 1.6)) +  15 * (std::pow(composition(1), 0.5) - std::pow(composition(1), 4.));
    PPC_SUTTON = (ppc_pseudo * (tpc_pseudo - epsilon)) / (tpc_pseudo + (composition(1) * (1 - composition(1)) * epsilon));
    return PPC_SUTTON;
}



double TPC_SUTTON(const double spgr, const arma::mat & nhc_properties, const std::string fluid) {

    int s = nhc_properties.n_rows;
    arma::vec composition(s);
    arma::vec mw(s);
    arma::vec tcrit(s);
    arma::vec pcrit(s);
    const double mw_air = 28.964;
    double y_HC, spgr_HC;
    double tcrit_HC, tpc_pseudo, epsilon, TPC_SUTTON;

    composition = nhc_properties.col(0);
    mw =  nhc_properties.col(1);
    tcrit =  nhc_properties.col(2);
    y_HC = 1 - sum(composition);
    spgr_HC = (spgr - sum(composition.t() * mw) / mw_air) / y_HC;
    if (fluid == "dry_gas") {
        tcrit_HC = 120.1 + 429. * spgr_HC  - 62.9 * std::pow(spgr_HC, 2.);
    } else {
        // "Gas Condensate, wet gas, rich gas
        tcrit_HC = 164.3 + 357.7 * spgr_HC  - 67.7 *  std::pow(spgr_HC, 2.);
    }
    tpc_pseudo = y_HC * tcrit_HC + sum(composition.t() * tcrit);
    epsilon = 120. * (std::pow(sum(composition.rows(1,2)), 0.9) - std::pow(sum(composition.rows(1,2)), 1.6)) +  15. * (std::pow(composition(1), 0.5) - std::pow(composition(1), 4.));
    TPC_SUTTON = tpc_pseudo - epsilon;
    return TPC_SUTTON;
}



double Z_FACTOR_DAK(double t, double p, double tpc, double ppc) {

    // Calculates gas compressibility factor using DAK method
    // "Fundamental PVT Calculations for Associated and Gas/Condensate Natural-Gas Systems",
    // SPE Reservoir Evaluation & Engineering, June 2007, 270-284, SPE 97099

    // t in R
    // p in psi
    // A values are based on Poettmann-Carpenter data, 5960 data points

    arma::vec A = {0.3265, -1.07, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.721};
    arma::vec coeff(5);
    double tr, pr, rho, error, g;
    double F_rho, DF_rho, Z_FACTOR_DAK;

    tr = t / tpc;
    pr = p / ppc;
    coeff(0) = A(0) + A(1) / tr + A(2) / std::pow(tr,3.) + A(3) / std::pow(tr,4.) + A(4) / std::pow(tr,5.);
    coeff(1) = A(5) + A(6) / tr + A(7) / std::pow(tr,2.);
    coeff(2) = -A(8) * (A(6) / tr + A(7) / std::pow(tr,2.));
    coeff(3) = A(9) / std::pow(tr,3.);
    coeff(4) = A(9) * A(10) / std::pow(tr,3.);
    rho = 0.27 * pr / tr;
    error = 1.;
    g = 0.27 * pr;
    for(int i = 0; i < 100; i++) {
        F_rho = (g / (tr * rho)) - coeff(0) * rho - coeff(1) * std::pow(rho,2.) - coeff(2) * std::pow(rho,5.) - coeff(3) * std::pow(rho,2.) * exp(-A(10) * std::pow(rho,2.)) - coeff(4) * std::pow(rho,4.) * exp(-A(10) * std::pow(rho,2.)) - 1.;
        DF_rho = (-g / (tr * std::pow(rho,2.))) - coeff(0) - 2. * coeff(1) * rho - 5. * coeff(2) * std::pow(rho,4.) - 2. * coeff(3) * rho * exp(-A(10) * std::pow(rho,2.)) + 2. * A(10) * rho * exp(-A(10) * std::pow(rho,2.)) * coeff(3) * std::pow(rho,2.) - 4. * coeff(4) * std::pow(rho,3.)* exp(-A(10) * std::pow(rho,2.)) + 2. * A(10) * rho * exp(-A(10) * std::pow(rho,2.)) * coeff(4) * std::pow(rho,4.);
        error = -F_rho / DF_rho;
        if ((i < 100) & (std::abs(error) < 1.e-13)) {
            break;
        } else if (i == 99) {
            rho = 123456.;
        } else {
            rho = rho + error;
        }
    }
    // while (std::abs(error) > 1e-13) {
    //     F_rho = (g / (tr * rho)) - coeff(0) * rho - coeff(1) * std::pow(rho,2.) - coeff(2) * std::pow(rho,5.) - coeff(3) * std::pow(rho,2.) * exp(-A(10) * std::pow(rho,2.)) - coeff(4) * std::pow(rho,4.) * exp(-A(10) * std::pow(rho,2.)) - 1.;
    //     DF_rho = (-g / (tr * std::pow(rho,2.))) - coeff(0) - 2. * coeff(1) * rho - 5. * coeff(2) * std::pow(rho,4.) - 2. * coeff(3) * rho * exp(-A(10) * std::pow(rho,2.)) + 2. * A(10) * rho * exp(-A(10) * std::pow(rho,2.)) * coeff(3) * std::pow(rho,2.) - 4. * coeff(4) * std::pow(rho,3.)* exp(-A(10) * std::pow(rho,2.)) + 2. * A(10) * rho * exp(-A(10) * std::pow(rho,2.)) * coeff(4) * std::pow(rho,4.);
    //     error = -F_rho / DF_rho;
    //     rho = rho + error;
    // }
    if (rho == 123456.) {
        Z_FACTOR_DAK = 1.;
    } else {
        Z_FACTOR_DAK = g / (rho * tr);
    }
    return(Z_FACTOR_DAK);
}



double DENSITY_GAS_DAK(double t, double p, double tpc, double ppc, double spgr) {

    // Calculates gas density using DAK method
    // "Fundamental PVT Calculations for Associated and Gas/Condensate Natural-Gas Systems",
    // SPE Reservoir Evaluation & Engineering, June 2007, 270-284, SPE 97099

    const double mw_air = 28.964;
    const double R = 10.73159;
    double mw = spgr * mw_air;
    double Z_factor = Z_FACTOR_DAK(t,p,tpc,ppc);
    double DENSITY_GAS_DAK = p * mw / Z_factor / R / t;         // lb/cuft
    return(DENSITY_GAS_DAK);
}



arma::vec B_GAS_DAK(double t, double p, double tsc, double psca, double tpc, double ppc) {

    // Calculates gas formation volume factor using DAK method
    // "Fundamental PVT Calculations for Associated and Gas/Condensate Natural-Gas Systems",
    // SPE Reservoir Evaluation & Engineering, June 2007, 270-284, SPE 97099

    // t in R
    // p in psi
    arma::vec B_GAS_DAK(2);
    double Z_factor = Z_FACTOR_DAK(t,p,tpc,ppc);
    B_GAS_DAK(0) = (Z_factor * t / p) / (tsc / psca);
    B_GAS_DAK(1) = B_GAS_DAK(0) / 5.615;      // rb/scf
    return(B_GAS_DAK);
}



double COMPRESSIBILITY_GAS_DAK(double t, double p, double tpc, double ppc) {

    // Calculates gas compressibility using DAK method
    // "Compressibility of natural gases Shawket",
    // Journal of Petroleum Science and Engineering, 10 ( 1993 ) 157-162

    // t in R
    // p in psi
    // A values are based on Poettmann-Carpenter data, 5960 data points

    arma::vec A = {0.3265, -1.07, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.721};
    arma::vec coeff(5);
    double tr, pr, rho, g;
    double Z_factor, COMPRESSIBILITY_GAS_DAK;

    tr = t / tpc;
    pr = p / ppc;
    g = 0.27 * pr;
    Z_factor = Z_FACTOR_DAK(t,p,tpc,ppc);
    rho = g / tr / Z_factor;
    coeff(0) = A(0) + A(1) / tr + A(2) / std::pow(tr,3.) + A(3) / std::pow(tr,4.) + A(4) / std::pow(tr,5.);
    coeff(1) = 2. * (A(5) + A(6) / tr + A(7) / std::pow(tr,2.)) * rho;
    coeff(2) = -5. * A(8) * (A(6) / tr + A(7) / std::pow(tr,2.)) *  std::pow(rho,4.);
    coeff(3) = (2. * A(9) * rho / std::pow(tr,3.) + 2 * A(9) * A(10) * std::pow(rho,3.) / std::pow(tr,3.) - 2 * A(9) * A(10) * A(10) * std::pow(rho,5.) / std::pow(tr,3.));
    coeff(4) = exp(-A(10) * std::pow(rho,2.));
    double dZ_drho = coeff(0) + coeff(1) + coeff(2) + coeff(3) * coeff(4);
    double Cr = (1. / pr) - (0.27 / std::pow(Z_factor,2.) / tr) * (dZ_drho / (1. + rho * dZ_drho / Z_factor));
    COMPRESSIBILITY_GAS_DAK = (Cr / ppc);                 // 1/Psi
    return(COMPRESSIBILITY_GAS_DAK);
}



double VISCOSITY_GAS_SUTTON(double t, double p, double tpc, double ppc, double spgr) {

    // Calculate pseudo critical pressure for gas condensate
    // sg is the reservoir gas gravity at or above the dew point pressure
    // Reference equations are (12),(14 for Epilison), (15) and (16)
    // "Fundamental PVT Calculations for Associated and Gas/Condensate Natural-Gas Systems",
    // SPE Reservoir Evaluation & Engineering, June 2007, 270-284, SPE 97099

    // t in R
    // p in psi
    const double mw_air = 28.964;
    double tr = t / tpc;
    double mw = spgr * mw_air;
    double rho = DENSITY_GAS_DAK(t,p,tpc,ppc,spgr) * 16.018463374/ 1000;      // gr/cm3
    double epsilon = 0.9490 * std::pow((tpc / std::pow(mw,3.) / std::pow(ppc,4.)),0.1666667);
    double visc_lowp = std::pow(10.,-4.) * (0.807 * std::pow(tr,0.618) - 0.357 * exp(-0.449 * tr) + 0.340 * exp(-4.058 * tr) + 0.018) / epsilon;
    double X = 3.47 + 1588 / t + 0.0009 * mw;
    double Y = 1.66378 - 0.004679 * X;
    double VISCOSITY_GAS_SUTTON = visc_lowp * exp(X * std::pow(rho,Y));
    return(VISCOSITY_GAS_SUTTON);
}



double PSEUDO_PRESSURE_GAS(double t, double p, double tpc, double ppc, double spgr) {

    // Calculates gas pseudo pressure
    // t in R
    // p in psi
    double psc = 0.;                    // Psig
    double psca = psc + 14.696;        // Psia
    double dp = 10.;
    int lp;

    arma::vec p_table_1 = arma::regspace<arma::vec>(psca, dp, p);
    if (p_table_1.max() < p) {
        lp = p_table_1.n_elem + 1;
    } else {
        lp = p_table_1.n_elem;
    }
    arma::vec p_table_imp(lp);
    arma::mat results_table(lp,5);
    if (p_table_1.max() < p) {
        p_table_imp.subvec(0,lp-2) = p_table_1;                  // Psia
        p_table_imp.subvec(lp-1,lp-1) = p;                       // Psia
    } else {
        p_table_imp.subvec(0,lp-1) = p_table_1;                  // Psia
    }
    for (int i = 0; i < lp; i++) {
        results_table(i,0) = p_table_imp(i);
        results_table(i,1) = Z_FACTOR_DAK(t,p_table_imp(i),tpc,ppc);
        results_table(i,2) = VISCOSITY_GAS_SUTTON(t,p_table_imp(i),tpc,ppc,spgr);
        results_table(i,3) = 2 * results_table(i,0) / results_table(i,1) / results_table(i,2);
        if (i == 0) {
            results_table(i,4) = 0.;
        } else {
            results_table(i,4) = results_table(i-1,4) + 0.5 * (results_table(i,3) + results_table(i-1,3)) * (results_table(i,0) - results_table(i-1,0));
        }
    }
    double sol = results_table(lp-1,4);
    return(sol);
}


arma::mat PVT_GAS_PROPERTIES_DAK_SUTTON(double t, arma::vec p, double spgr, double tpc, double ppc) {

    double psca = 14.696;              // Psia
    double tsc = (60. + 459.67);        // R
    double term_old = 0, term_new;
    int lp = p.size();
    arma::mat results_table(lp,6);
    for (int i = 0; i < lp; i++) {
        results_table(i,0) = Z_FACTOR_DAK(t,p(i),tpc,ppc);
        results_table(i,1) = B_GAS_DAK(t,p(i),tsc,psca,tpc,ppc)(1);
        results_table(i,2) = DENSITY_GAS_DAK(t,p(i),tpc,ppc,spgr);
        results_table(i,3) = COMPRESSIBILITY_GAS_DAK(t,p(i),tpc,ppc);
        results_table(i,4) = VISCOSITY_GAS_SUTTON(t,p(i),tpc,ppc,spgr);
        term_new = 2 * p(i) / results_table(i,0) / results_table(i,4);
        if (i == 0) {
            results_table(i,5) = 0.;
        } else {
            results_table(i,5) = results_table(i-1,5) + 0.5 * (term_new + term_old) * (p(i) - p(i-1));
        }
        term_old = term_new;
    }
    return(results_table);
}

