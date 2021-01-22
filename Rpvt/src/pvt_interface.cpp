// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us

#include "RcppArmadillo.h"
#include "pvtgas_header.h"
#include "pvtoil_header.h"
#include "pvtwater_header.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// **************** GAS PVT INTERFACE ******************

// [[Rcpp::export]]
NumericMatrix PVT_GAS_PROPERTIES(Rcpp::List lst) {
    std::string input_unit = as<std::string>(lst["input_unit"]);
    std::string output_unit = as<std::string>(lst["output_unit"]);
    std::string fluid = as<std::string>(lst["fluid"]);
    std::string pvt_model = as<std::string>(lst["pvt_model"]);
    std::string visc_model = as<std::string>(lst["visc_model"]);
    double t = as<double>(lst["t"]);
    double p = as<double>(lst["p"]);
    double spgr = as<double>(lst["gas_spgr"]);
    arma::vec nhc_composition = as<arma::vec>(lst["nhc_composition"]);
    double cgr = as<double>(lst["cgr"]);
    double api = as<double>(lst["cond_api"]);

    arma::mat nhc_properties(3,4);
    double p_imp; double t_imp; double psc;
    double mw_cond; double spgr_cond; double gor;
    double dp = 10.; double p_min = 0.; int lp;
    double kpa_to_psi = 0.14503773800722;
    double m3_to_scf = 35.314667;
    double m3_to_bbl = 6.289814;
    double kg_m3_to_lb_cuft = 0.06242796;
    double scf_stb_to_sm3_sm3 = 0.1781077;

    arma::vec nhc_mw = {28.01, 34.08, 44.01};
    arma::vec nhc_tcrit = {226.98, 672.35, 547.54};
    arma::vec nhc_pcrit = {492.26, 1299.97, 1070.67};
    nhc_properties.col(0) = nhc_composition;
    nhc_properties.col(1) = nhc_mw;
    nhc_properties.col(2) = nhc_tcrit;
    nhc_properties.col(3) = nhc_pcrit;

    if (input_unit == "SI") {
        p_imp = p * kpa_to_psi;     // kPag to Psig
        t_imp = t * 1.8 + 32.;       // C to F
        t_imp = t_imp + 459.67;     // F to R
        psc = 0.;                    // Psig
        cgr = cgr * (scf_stb_to_sm3_sm3) * 1e06;     // STB/MMSCF
    } else {
        p_imp = p;                  // Psig
        t_imp = t + 459.67;         // F to R
        psc = 0.;                    // Psig
    }
    if (cgr == 0.) {
    } else {
        mw_cond = 6084. / (api - 5.9);
        spgr_cond = 141.5 / (131.5 + api);
        gor = 1.e6 / cgr;
        spgr = (gor * spgr + 4580. * spgr_cond) / (gor + 133000. * spgr_cond / mw_cond);
    }
    double tpc = TPC_SUTTON(spgr,nhc_properties,fluid);
    double ppc = PPC_SUTTON(spgr,nhc_properties,fluid);

    arma::vec p_table_1 = regspace< vec>(std::min(p_min, p_imp), dp, p_imp);
    if ((p_table_1.min() > psc) & (p_table_1.max() < p_imp)) {
        lp = p_table_1.n_elem + 2;
    } else if ((p_table_1.min() > psc) & (p_table_1.max() == p_imp)) {
        lp = p_table_1.n_elem + 1;
    } else if ((p_table_1.min() == psc) & (p_table_1.max() < p_imp)) {
        lp = p_table_1.n_elem + 1;
    } else {
        lp = p_table_1.n_elem;
    }
    arma::vec p_table_imp(lp);
    if ((p_table_1.min() > psc) & (p_table_1.max() < p_imp)) {
        p_table_imp.subvec(0,0) = psc + 14.696;                          // Psig to Psia
        p_table_imp.subvec(1,lp-2) = p_table_1 + 14.696;                 // Psig to Psia
        p_table_imp.subvec(lp-1,lp-1) = p_imp + 14.696;                  // Psig to Psia
    } else if ((p_table_1.min() > psc) & (p_table_1.max() == p_imp)) {
        p_table_imp.subvec(0,0) = psc + 14.696;                          // Psig to Psia
        p_table_imp.subvec(1,lp-1) = p_table_1 + 14.696;                 // Psig to Psia
    } else if ((p_table_1.min() == psc) & (p_table_1.max() < p_imp)) {
        p_table_imp.subvec(0,lp-2) = p_table_1 + 14.696;                 // Psig to Psia
        p_table_imp.subvec(lp-1,lp-1) = p_imp + 14.696;                  // Psig to Psia
    } else {
        p_table_imp.subvec(0,lp-1) = p_table_1 + 14.696;                 // Psig to Psia
    }

    arma::mat pvt_gas_result(lp,6);
    if ((pvt_model == "DAK") && (visc_model == "Sutton")) {
       pvt_gas_result = PVT_GAS_PROPERTIES_DAK_SUTTON(t_imp,p_table_imp,spgr,tpc,ppc);
    } else {
       pvt_gas_result = PVT_GAS_PROPERTIES_DAK_SUTTON(t_imp,p_table_imp,spgr,tpc,ppc);
    }

    arma::mat results_table(lp,2);
    NumericMatrix results_table_(lp,8);
    CharacterVector colname(8);
    for (int i = 0; i < lp; i++) {
        results_table(i,0) = t_imp - 459.67;
        results_table(i,1) = p_table_imp(i) - 14.696;
    }
    results_table.insert_cols(results_table.n_cols,pvt_gas_result);
    if (output_unit == "SI") {
        for (int i = 0; i < lp; i++) {
        results_table(i,0) = (results_table(i,0)-32.)/1.8;                               // F to C
        results_table(i,1) = results_table(i,1) / kpa_to_psi;                           // Psig to kPag
        results_table(i,3) = results_table(i,3) / kg_m3_to_lb_cuft;                     // lb/cuft to kg/m3
        results_table(i,4) = results_table(i,4) * m3_to_scf / m3_to_bbl;                // rm3/Sm3
        results_table(i,5) = results_table(i,5) * kpa_to_psi;                           // 1/Psi to 1/kPa
        results_table(i,7) = results_table(i,7) / kpa_to_psi / kpa_to_psi / 1.e06;       // Psi^2/cp to MPa^2/mPa.s
        }
    }
    results_table_ = wrap(results_table);                                               // converting arma::mat to Rcpp::NumericMatrix
    if (output_unit == "SI") {
        colname = {"T_(C)", "P_(kPag)", "Z-Factor", "Bg_(rm3/sm3)", "Density_(kg/m3)", "Cg_(1/kPaa)", "Viscosity_(mPa.s)", "m(p)_(MPaa^2/mPa.s)"};
        colnames(results_table_) = colname;
    } else {
        colname = {"T_(F)", "P_(Psig)", "Z-Factor", "Bg_(rb/scf)", "Density_(lb/cuft)", "Cg_(1/Psia)", "Viscosity_(cp)", "m(p)_(Psia^2/cp)"};
        colnames(results_table_) = colname;
    }
    if (output_unit == "SI") {
        results_table_.attr("gas pseudocritical temperature (C)") = tpc / 1.8 - 273.15;
        results_table_.attr("gas pseudocritical pressure (kPaa)") = ppc / kpa_to_psi;
    }
    if (output_unit == "Field") {
        results_table_.attr("gas pseudocritical temperature (F)") = tpc - 459.67;
        results_table_.attr("gas pseudocritical pressure (Psia)") = ppc;
    }
    return(results_table_);
}




// **************** OIL PVT INTERFACE ******************

// [[Rcpp::export]]
NumericMatrix PVT_OIL_PROPERTIES(Rcpp::List lst) {

    std::string input_unit = as<std::string>(lst["input_unit"]);
    std::string output_unit = as<std::string>(lst["output_unit"]);
    std::string fluid = as<std::string>(lst["fluid"]);
    std::string pvt_model = as<std::string>(lst["pvt_model"]);
    std::string visc_model = as<std::string>(lst["visc_model"]);
    double t = as<double>(lst["t"]);
    double p = as<double>(lst["p"]);
    double api = as<double>(lst["oil_api"]);
    double spgr = as<double>(lst["gas_spgr"]);
    arma::vec nhc_composition = as<arma::vec>(lst["nhc_composition"]);
    double sat_cond = as<double>(lst["sat_cond"]);
    int flag_sat = as<double>(lst["flag_sat"]);

    arma::mat nhc_properties(3,4);
    double p_imp; double t_imp; double psc;
    double dp = 10.; double p_min = 0.; int lp;
    double kpa_to_psi = 0.14503773800722;
    double m3_to_scf = 35.314667;
    double m3_to_bbl = 6.289814;
    double kg_m3_to_lb_cuft = 0.06242796;
    double scf_stb_to_sm3_sm3 = 0.1781077;
    double pb = 0.; double rsi = 0.;
    String gas_type = "dry_gas";

    arma::vec nhc_mw = {28.01, 34.08, 44.01};
    arma::vec nhc_tcrit = {226.98, 672.35, 547.54};
    arma::vec nhc_pcrit = {492.26, 1299.97, 1070.67};
    nhc_properties.col(0) = nhc_composition;
    nhc_properties.col(1) = nhc_mw;
    nhc_properties.col(2) = nhc_tcrit;
    nhc_properties.col(3) = nhc_pcrit;

    if (input_unit == "SI") {
        p_imp = p * kpa_to_psi;     // kPag to Psig
        t_imp = t * 1.8 + 32.;       // C to F
        t_imp = t_imp + 459.67;     // F to R
        psc = 0.;                    // Psig
    } else {
        p_imp = p;                  // Psig
        t_imp = t + 459.67;         // F to R
        psc = 0.;                    // Psig
    }
    double tpc = TPC_SUTTON(spgr,nhc_properties,gas_type);
    double ppc = PPC_SUTTON(spgr,nhc_properties,gas_type);

    arma::vec p_table_1 = regspace< vec>(std::min(p_min, p_imp), dp, p_imp);
    if ((p_table_1.min() > psc) & (p_table_1.max() < p_imp)) {
        lp = p_table_1.n_elem + 2;
    } else if ((p_table_1.min() > psc) & (p_table_1.max() == p_imp)) {
        lp = p_table_1.n_elem + 1;
    } else if ((p_table_1.min() == psc) & (p_table_1.max() < p_imp)) {
        lp = p_table_1.n_elem + 1;
    } else {
        lp = p_table_1.n_elem;
    }
    arma::vec p_table_imp(lp);
    if ((p_table_1.min() > psc) & (p_table_1.max() < p_imp)) {
        p_table_imp.subvec(0,0) = psc + 14.696;                          // Psig to Psia
        p_table_imp.subvec(1,lp-2) = p_table_1 + 14.696;                 // Psig to Psia
        p_table_imp.subvec(lp-1,lp-1) = p_imp + 14.696;                  // Psig to Psia
    } else if ((p_table_1.min() > psc) & (p_table_1.max() == p_imp)) {
        p_table_imp.subvec(0,0) = psc + 14.696;                          // Psig to Psia
        p_table_imp.subvec(1,lp-1) = p_table_1 + 14.696;                 // Psig to Psia
    } else if ((p_table_1.min() == psc) & (p_table_1.max() < p_imp)) {
        p_table_imp.subvec(0,lp-2) = p_table_1 + 14.696;                 // Psig to Psia
        p_table_imp.subvec(lp-1,lp-1) = p_imp + 14.696;                  // Psig to Psia
    } else {
        p_table_imp.subvec(0,lp-1) = p_table_1 + 14.696;                 // Psig to Psia
    }

    if (flag_sat == 1) {
        if (input_unit == "SI") {
            pb = sat_cond * kpa_to_psi;     // kPag to Psig
        } else {
            pb = sat_cond;
        }
        pb = pb + 14.696;
        if (pvt_model == "Standing") {
            rsi = RSi_STANDING(t_imp,pb,api,spgr)(0);
        }
        if (pvt_model == "Vasquez_Beggs") {
            rsi = RSi_VASQUEZ_BEGGS(t_imp,pb,api,spgr)(0);
        }
        if (pvt_model == "Farshad_Petrosky") {
            rsi = RSi_FARSHAD_PETROSKY(t_imp,pb,api,spgr)(0);
        }
        if (pvt_model == "Al_Marhoun") {
            rsi = RSi_AL_MARHOUN(t_imp,pb,api,spgr)(0);
        }
        if (pvt_model == "Glaso") {
            rsi = RSi_GLASO(t_imp,pb,api,spgr)(0);
        }
    } else {
        if (input_unit == "SI") {
            rsi = sat_cond / scf_stb_to_sm3_sm3;     // m3/m3 to to SCF/STB
        } else {
            rsi = sat_cond;
        }
        if (pvt_model == "Standing") {
            pb = PB_STANDING(t_imp,api,spgr,rsi);
        }
        if (pvt_model == "Vasquez_Beggs") {
            pb = PB_VASQUEZ_BEGGS(t_imp,api,spgr,rsi);
        }
        if (pvt_model == "Farshad_Petrosky") {
            pb = PB_FARSHAD_PETROSKY(t_imp,api,spgr,rsi);
        }
        if (pvt_model == "Al_Marhoun") {
            pb = PB_AL_MARHOUN(t_imp,api,spgr,rsi);
        }
        if (pvt_model == "Glaso") {
            pb = PB_GLASO(t_imp,api,spgr,rsi);
        }
    }
    arma::mat pvt_oil_result(lp,4);
    arma::vec visc_oil_result(lp);
    arma::mat pvt_gas_result(lp,6);

    if (pvt_model == "Standing") {
        pvt_oil_result = PVT_OIL_PROPERTIES_STANDING(t_imp,p_table_imp,spgr,api,rsi,tpc,ppc);
    }
    if (pvt_model == "Vasquez_Beggs") {
        pvt_oil_result = PVT_OIL_PROPERTIES_VASQUEZ_BEGGS(t_imp,p_table_imp,spgr,api,rsi,tpc,ppc);
    }
    if (pvt_model == "Farshad_Petrosky") {
        pvt_oil_result = PVT_OIL_PROPERTIES_FARSHAD_PETROSKY(t_imp,p_table_imp,spgr,api,rsi,tpc,ppc);
    }
    if (pvt_model == "Al_Marhoun") {
        pvt_oil_result = PVT_OIL_PROPERTIES_AL_MARHOUN(t_imp,p_table_imp,spgr,api,rsi,tpc,ppc);
    }
    if (pvt_model == "Glaso") {
        pvt_oil_result = PVT_OIL_PROPERTIES_GLASO(t_imp,p_table_imp,spgr,api,rsi,tpc,ppc);
    }

    if(visc_model == "Beggs_Robinson") {
        visc_oil_result = VISC_OIL_PROPERTIES_BEGGS_ROBINSON(pvt_model,t_imp,p_table_imp,spgr,api,rsi);
    }
    if(visc_model == "Al_Marhoun") {
        visc_oil_result = VISC_OIL_PROPERTIES_AL_MARHOUN(pvt_model,t_imp,p_table_imp,spgr,api,rsi);
    }
    pvt_gas_result = PVT_GAS_PROPERTIES_DAK_SUTTON(t_imp,p_table_imp,spgr,tpc,ppc);

    arma::mat results_table(lp,2);
    NumericMatrix results_table_(lp,13);
    CharacterVector colname(13);
    for (int i = 0; i < lp; i++) {
        results_table(i,0) = t_imp - 459.67;
        results_table(i,1) = p_table_imp(i) - 14.696;
    }
    results_table.insert_cols(results_table.n_cols,pvt_oil_result);
    results_table.insert_cols(results_table.n_cols,visc_oil_result);
    results_table.insert_cols(results_table.n_cols,pvt_gas_result);
    if (output_unit == "SI") {
        for (int i = 0; i < lp; i++) {
            results_table(i,0) = (results_table(i,0)-32.)/1.8;                               // F to C
            results_table(i,1) = results_table(i,1) / kpa_to_psi;                           // Psig to kPag
            results_table(i,2) = results_table(i,2) * scf_stb_to_sm3_sm3;                   // scf/stb to sm3/sm3
            results_table(i,4) = results_table(i,4) / kg_m3_to_lb_cuft;                     // lb/cuft to kg/m3
            results_table(i,5) = results_table(i,5) * kpa_to_psi;                           // 1/Psi to 1/kPa
            results_table(i,8) = results_table(i,8) * (m3_to_scf / m3_to_bbl);              // rm3/Sm3
            results_table(i,9) = results_table(i,9) / kg_m3_to_lb_cuft;                     // lb/cuft to kg/m3
            results_table(i,10) = results_table(i,10) * kpa_to_psi;                         // 1/Psi to 1/kPa
            results_table(i,12) = results_table(i,12) / kpa_to_psi / kpa_to_psi / 1e6;      // MPsi^2/cp to MPa^2/mPa.s
        }
    }
    results_table_ = wrap(results_table);                                               // converting arma::mat to Rcpp::NumericMatrix
    if (output_unit == "SI") {
        colname = {"T_(C)", "P_(kPag)", "Rso_(rm3/sm3)", "Bo_(rm3/sm3)", "Oil_Density_(kg/m3)", "Co_(1/kPaa)", "Oil_Viscosity_(mPa.s)", "Z-Factor", "Bg_(rm3/sm3)", "Gas_Density_(kg/m3)", "Cg_(1/kPaa)", "Gas_Viscosity_(mPa.s)", "m(p)_(MPaa^2/mPa.s)"};
        colnames(results_table_) = colname;
    } else {
        colname = {"T_(F)", "P_(Psig)", "Rso_(scf/stb)", "Bo_(rb/stb)", "Oil_Density_(lb/ft3)", "Co_(1/Psia)", "Oil_Viscosity_(cp)", "Z-Factor", "Bg_(rb/scf)", "Gas_Density_(lb/ft3)", "Cg_(1/Psia)", "Gas_Viscosity_(cp)", "m(p)_(Psia^2/cp)"};
        colnames(results_table_) = colname;
    }
    if (output_unit == "SI") {
        results_table_.attr("gas pseudocritical temperature (C)") = tpc / 1.8 - 273.15;
        results_table_.attr("gas pseudocritical pressure (kPaa)") = ppc / kpa_to_psi;
    }
    if (output_unit == "Field") {
        results_table_.attr("gas pseudocritical temperature (F)") = tpc - 459.67;
        results_table_.attr("gas pseudocritical pressure (Psia)") = ppc;
    }
    return(results_table_);
}



// **************** WATER PVT INTERFACE ******************


// [[Rcpp::export]]
NumericMatrix PVT_WATER_PROPERTIES(Rcpp::List lst) {

    std::string input_unit = as<std::string>(lst["input_unit"]);
    std::string output_unit = as<std::string>(lst["output_unit"]);
    std::string fluid = as<std::string>(lst["fluid"]);
    std::string pvt_model = as<std::string>(lst["pvt_model"]);
    std::string visc_model = as<std::string>(lst["visc_model"]);
    std::string gas_saturated = as<std::string>(lst["gas_saturated"]);
    double t = as<double>(lst["t"]);
    double p = as<double>(lst["p"]);
    double S = as<double>(lst["salinity"]);

    arma::mat nhc_properties(3,4);
    double p_imp; double t_imp; double psc;
    double dp = 10.; double p_min = 0.; int lp;
    double kpa_to_psi = 0.14503773800722;
    double m3_to_scf = 35.314667;
    double m3_to_bbl = 6.289814;
    double kg_m3_to_lb_cuft = 0.06242796;
    double spgr = 0.63;
    String gas_type = "dry_gas";

    arma::vec nhc_comp = {0.0, 0.0, 0.0};
    arma::vec nhc_mw = {28.01, 34.08, 44.01};
    arma::vec nhc_tcrit = {226.98, 672.35, 547.54};
    arma::vec nhc_pcrit = {492.26, 1299.97, 1070.67};
    nhc_properties.col(0) = nhc_comp;
    nhc_properties.col(1) = nhc_mw;
    nhc_properties.col(2) = nhc_tcrit;
    nhc_properties.col(3) = nhc_pcrit;

    if (input_unit == "SI") {
        p_imp = p * kpa_to_psi;     // kPag to Psig
        t_imp = t * 1.8 + 32.;       // C to F
        t_imp = t_imp + 459.67;     // F to R
        psc = 0.;                    // Psig
    } else {
        p_imp = p;                  // Psig
        t_imp = t + 459.67;         // F to R
        psc = 0.;                    // Psig
    }

    double tpc = TPC_SUTTON(spgr,nhc_properties,gas_type);
    double ppc = PPC_SUTTON(spgr,nhc_properties,gas_type);

    arma::vec p_table_1 = regspace< vec>(std::min(p_min, p_imp), dp, p_imp);
    if ((p_table_1.min() > psc) & (p_table_1.max() < p_imp)) {
        lp = p_table_1.n_elem + 2;
    } else if ((p_table_1.min() > psc) & (p_table_1.max() == p_imp)) {
        lp = p_table_1.n_elem + 1;
    } else if ((p_table_1.min() == psc) & (p_table_1.max() < p_imp)) {
        lp = p_table_1.n_elem + 1;
    } else {
        lp = p_table_1.n_elem;
    }
    arma::vec p_table_imp(lp);
    if ((p_table_1.min() > psc) & (p_table_1.max() < p_imp)) {
        p_table_imp.subvec(0,0) = psc + 14.696;                          // Psig to Psia
        p_table_imp.subvec(1,lp-2) = p_table_1 + 14.696;                 // Psig to Psia
        p_table_imp.subvec(lp-1,lp-1) = p_imp + 14.696;                  // Psig to Psia
    } else if ((p_table_1.min() > psc) & (p_table_1.max() == p_imp)) {
        p_table_imp.subvec(0,0) = psc + 14.696;                          // Psig to Psia
        p_table_imp.subvec(1,lp-1) = p_table_1 + 14.696;                 // Psig to Psia
    } else if ((p_table_1.min() == psc) & (p_table_1.max() < p_imp)) {
        p_table_imp.subvec(0,lp-2) = p_table_1 + 14.696;                 // Psig to Psia
        p_table_imp.subvec(lp-1,lp-1) = p_imp + 14.696;                  // Psig to Psia
    } else {
        p_table_imp.subvec(0,lp-1) = p_table_1 + 14.696;                 // Psig to Psia
    }

    arma::mat pvt_water_result(lp,4);
    arma::vec visc_water_result(lp);
    if (pvt_model == "Meehan") {
        pvt_water_result = PVT_WATER_PROPERTIES_MEEHAN(gas_saturated,t_imp,p_table_imp,S,tpc,ppc);
    }
    if (pvt_model == "McCain") {
        pvt_water_result = PVT_WATER_PROPERTIES_MCCAIN(gas_saturated,t_imp,p_table_imp,S,tpc,ppc);
    }
    if (pvt_model == "Spivey") {
        pvt_water_result = PVT_WATER_PROPERTIES_SPIVEY(gas_saturated,t_imp,p_table_imp,S,tpc,ppc);
    }
    if (visc_model == "Meehan") {
        visc_water_result = VISC_WATER_PROPERTIES_MEEHAN(t_imp,p_table_imp,S);
    }
    if (visc_model == "McCain") {
        visc_water_result = VISC_WATER_PROPERTIES_MCCAIN(t_imp,p_table_imp,S);
    }
    if (visc_model == "Spivey") {
        visc_water_result = VISC_WATER_PROPERTIES_SPIVEY(t_imp,p_table_imp,S);
    }


    arma::mat results_table(lp,2);
    NumericMatrix results_table_(lp,7);
    CharacterVector colname(7);
    for (int i = 0; i < lp; i++) {
        results_table(i,0) = t_imp - 459.67;
        results_table(i,1) = p_table_imp(i) - 14.696;
    }
    results_table.insert_cols(results_table.n_cols,pvt_water_result);
    results_table.insert_cols(results_table.n_cols,visc_water_result);
    if (output_unit == "SI") {
        for (int i = 0; i < lp; i++) {
            results_table(i,0) = (results_table(i,0)-32.)/1.8;                               // F to C
            results_table(i,1) = results_table(i,1) / kpa_to_psi;                           // Psig to kPag
            results_table(i,2) = results_table(i,2) / (m3_to_scf / m3_to_bbl);              // scf/STB to sm3/sm3
            results_table(i,4) = results_table(i,4) / kg_m3_to_lb_cuft;                     // lb/cuft to kg/m3
            results_table(i,5) = results_table(i,5) * kpa_to_psi;                           // 1/Psi to 1/kPa
        }
    }
    results_table_ = wrap(results_table);                        // converting arma::mat to Rcpp::NumericMatrix
    if (output_unit == "SI") {
        colname = {"T_(C)", "P_(kPag)", "Rsw_(rm3/sm3)", "Bw_(rm3/sm3)", "Density_(kg/m3)", "Cw_(1/kPaa)", "Viscosity_(mPa.s)"};
        colnames(results_table_) = colname;
    } else {
        colname = {"T_(F)", "P_(Psig)", "Rsw_(scf/stb)", "Bw_(rb/stb)", "Density_(lb/ft3)", "Cw_(1/Psia)", "Viscosity_(cp)"};
        colnames(results_table_) = colname;
    }
    return(results_table_);
}
