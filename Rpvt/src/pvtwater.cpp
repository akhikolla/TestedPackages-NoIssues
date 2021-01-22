// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "pvtgas_header.h"



// ************************************** WATER PVT CORRELATIONS ****************************************

// ************ MCCAIN CORRELATIONS ************

arma::vec B_WATER_MCCAIN(double t, double p) {

    // "Reservoir Fluid Property Correlations-State of the Art",
    // SPE Reservoir Engineering, May 1991 , 266-272

    // t in R
    // p in psia
    t = t - 459.67;     // R to F
    double dVp;
    double dVpp;
    double dVt;
    arma::vec B_WATER(2);
    dVp = -1.95301e-9 * p * t - 1.72834e-13 * p * p * t - 3.58922e-7 * p - 2.25341e-10 * p * p;
    dVpp = -1.95301e-9 * t - 2 * 1.72834e-13 * p * t - 3.58922e-7 - 2 * 2.25341e-10 * p;
    dVt = -1.0001e-2 + 1.33391e-4 * t + 5.50654e-7 * t * t;
    B_WATER(0) = (1+dVp) * (1+dVt);                     // rb/STB
    B_WATER(1) = dVpp * (1+dVt);                       // rb/STB/psia
    return(B_WATER);
}


double DENSITY_WATER_MCCAIN(double t, double p, double S) {

    // "Reservoir Fluid Property Correlations-State of the Art",
    // SPE Reservoir Engineering, May 1991 , 266-272

    // t in R
    // p in psia
    // S salinity in weight percent TDS
    // 1 wt% = 10000 ppm
    // 1 ppm = 1 mg/lit
    double Bw = B_WATER_MCCAIN(t,p)(0);                         // rb/STB
    double density_std;
    double DENSITY_WATER;
    density_std = 62.368 + 0.438603 * S + 1.60074e-3 * S * S;   // lbm/ft3
    DENSITY_WATER = density_std / Bw;             // lbm/ft3
    return(DENSITY_WATER);
}


arma::vec Rs_WATER_MCCAIN(const std::string gas_saturated, double t, double p, double S) {

    // "Reservoir Fluid Property Correlations-State of the Art",
    // SPE Reservoir Engineering, May 1991 , 266-272

    // t in R
    // p in psia
    // S salinity in weight percent TDS
    // 1 wt% = 10000 ppm
    // 1 ppm = 1 mg/lit
    t = t - 459.67;     // R to F
    double A;
    double B;
    double C;
    arma::vec Rs_WATER(2);
    if (gas_saturated == "yes") {
        A = 8.15839 - 6.12265e-2 * t + 1.91663e-4 * t * t - 2.1654e-7 * t * t * t;
        B = 1.01021e-2 - 7.44241e-5 * t + 3.05553e-7 * t * t - 2.94883e-10 * t * t * t;
        C = -1e-7 * (9.02505 - 0.130237 * t + 8.53425e-4 * t * t - 2.34122e-6 * t * t * t + 2.37049e-9 * t * t * t * t);
        double Rs_water_pure = A + B * p + C * p * p;
        double dRs_water_pure_dp = B + 2 * C * p;
        Rs_WATER(0) = Rs_water_pure * pow(10,-0.0840655 * S * pow(t,-0.285854));      // scf/STB
        Rs_WATER(1) = dRs_water_pure_dp * pow(10,-0.0840655 * S * pow(t,-0.285854));      // scf/STB/psia
    } else {
        Rs_WATER(0) = 0;      // scf/STB
        Rs_WATER(1) = 0;      // scf/STB/psia
    }

    return(Rs_WATER);
}


double COMPRESSIBILITY_WATER_MCCAIN(const std::string gas_saturated, double t, double p, double S, double Bg) {

    // "Reservoir Fluid Property Correlations-State of the Art",
    // SPE Reservoir Engineering, May 1991 , 266-272

    // t in R
    // p in psia
    // S salinity in weight percent TDS
    // 1 wt% = 10000 ppm
    // 1 ppm * rho = 1 mg/lit
    double tsc = 60;
    double psc = 14.696;
    double kg_m3_to_lb_cuft = 0.06242796;
    double rho_w = DENSITY_WATER_MCCAIN(tsc,psc,S) / kg_m3_to_lb_cuft / 1000;      // lbm/ft3 to g/cm3
    double Bw = B_WATER_MCCAIN(t,p)(0);                         // rb/STB
    // double dBw_dp = B_WATER(t,p)(1);                     // rb/STB/psia
    double dRsw_dp = Rs_WATER_MCCAIN(gas_saturated,t,p,S)(1);                 // scf/STB/psia
    t = t - 459.67;                      // R to F
    double S1 = S * rho_w * 10000;      // wt% to mg/lit
    double Cw1;
    double Cw2;
    double COMPRESSIBILITY_WATER;
    if (gas_saturated == "yes") {
        Cw1 = 1 / (7.033 * p + 0.5415 * S1 - 537 * t + 403300);
        // Cw1 = -1 / Bw * dBw_dp;
        Cw2 = Bg / Bw * dRsw_dp;
    } else {
        Cw1 = 1 / (7.033 * p + 0.5415 * S1 - 537 * t + 403300);
        // Cw1 = -1 / Bw * dBw_dp;
        Cw2 = 0;
    }
    COMPRESSIBILITY_WATER = Cw1 + Cw2;
    return(COMPRESSIBILITY_WATER);
}


double VISCOSITY_WATER_MCCAIN(double t, double p, double S) {

    // "Reservoir Fluid Property Correlations-State of the Art",
    // SPE Reservoir Engineering, May 1991 , 266-272

    // t in R
    // p in psia
    // S salinity in weight percent TDS
    t = t - 459.67;     // R to F
    double MUw;
    double MUw1;
    double A;
    double B;
    double VISCOSITY_WATER;

    A = 109.574 - 8.40564 * S + 0.313314 * S * S + 8.72213e-3 * S * S * S;
    B = 1.12166 - 2.63951e-2 * S + 6.79461e-4 * S * S + 5.47119e-5 * S * S * S - 1.55586e-6 * S * S * S * S;
    MUw1 = A * pow(t,-B);
    MUw = MUw1 * (0.9994 + 4.0295e-5 * p + 3.1062e-9 * p * p);
    VISCOSITY_WATER = MUw;                                  // cp
    return(VISCOSITY_WATER);
}


arma::mat PVT_WATER_PROPERTIES_MCCAIN(const std::string gas_saturated, double t, arma::vec p, double S, double tpc, double ppc) {

    double psca = 14.696;              // Psia
    double tsc = (60 + 459.67);        // R
    double Bg;
    int lp = p.size();
    arma::mat results_table(lp,4);

    for (int i = 0; i < lp; i++) {
        results_table(i,0) = Rs_WATER_MCCAIN(gas_saturated,t,p(i),S)(0);                          // rb/STB
        results_table(i,1) = B_WATER_MCCAIN(t,p(i))(0);
        results_table(i,2) = DENSITY_WATER_MCCAIN(t,p(i),S);
        Bg = B_GAS_DAK(t,p(i),tsc,psca,tpc,ppc)(1);                                                 // rb/STB
        results_table(i,3) = COMPRESSIBILITY_WATER_MCCAIN(gas_saturated,t,p(i),S,Bg);
    }
    return(results_table);
}


arma::vec VISC_WATER_PROPERTIES_MCCAIN(double t, arma::vec p, double S) {

    int lp = p.size();
    arma::vec results_table(lp);
    for (int i = 0; i < lp; i++) {
        results_table(i) = VISCOSITY_WATER_MCCAIN(t,p(i),S);
    }
    return(results_table);
}



// ************ MEEHAN CORRELATIONS ************


double B_WATER_MEEHAN(const std::string gas_saturated, double t, double p, double S) {

    // t in R
    // p in psia
    // S salinity in weight percent TDS
    // 1 wt% = 10000 ppm
    t = t - 459.67;     // R to F
    double a; double b; double c;
    double Sc;
    double Bw_Meehan;
    Sc = 1 + S * (5.1e-8 * p + (5.47e-6 - 1.96e-18 * p) * (t - 60) + (-3.23e-8 + 8.5e-13 * p) * (t - 60) * (t - 60));
    if (gas_saturated == "yes") {
        a = 0.9911 + 6.35e-6 * t + 8.5e-7 * t * t;
        b = -1.093e-6 - 3.497e-9 * t + 4.57e-12 * t * t;
        c = -5e-11 + 6.429e-13 * t - 1.43e-15 * t * t;
    } else {
        a = 0.9947 + 5.8e-6 * t + 1.02e-6 * t * t;
        b = -4.228e-6 + 1.8376e-8 * t - 6.77e-11 * t * t;
        c = 1.3e-10 - 1.3855e-12 * t + 4.285e-15 * t * t;
    }
    Bw_Meehan = (a + b * p + c * p * p) * Sc;                     // rb/STB
    return(Bw_Meehan);
}


double DENSITY_WATER_MEEHAN(const std::string gas_saturated, double t, double p, double S) {

    // t in R
    // p in psia
    // S salinity in weight percent TDS
    // 1 wt% = 10000 ppm
    // 1 ppm = 1 mg/lit
    double Bw = B_WATER_MEEHAN(gas_saturated,t,p,S) ;                // rb/STB
    double density_std;
    double DENSITY_WATER_MEEHAN;
    density_std = 62.368 + 0.438603 * S + 1.60074e-3 * S * S;              // lbm/ft3
    DENSITY_WATER_MEEHAN = density_std / Bw;                           // lbm/ft3
    return(DENSITY_WATER_MEEHAN);
}


arma::vec Rs_WATER_MEEHAN(const std::string gas_saturated, double t, double p, double S) {

    // t in R
    // p in psia
    // S salinity in weight percent TDS
    // 1 wt% = 10000 ppm
    // 1 ppm = 1 mg/lit
    t = t - 459.67;     // R to F
    double A;
    double B;
    double C;
    arma::vec Rs_WATER(2);
    if (gas_saturated == "yes") {
        A = 8.15839 - 6.12265e-2 * t + 1.91663e-4 * t * t - 2.1654e-7 * t * t * t;
        B = 1.01021e-2 - 7.44241e-5 * t + 3.05553e-7 * t * t - 2.94883e-10 * t * t * t;
        C = -1e-7 * (9.02505 - 0.130237 * t + 8.53425e-4 * t * t - 2.34122e-6 * t * t * t + 2.37049e-9 * t * t * t * t);
        double Rs_water_pure = A + B * p + C * p * p;
        double dRs_water_pure_dp = B + 2 * C * p;
        Rs_WATER(0) = Rs_water_pure * pow(10,-0.0840655 * S * pow(t,-0.285854));          // scf/STB
        Rs_WATER(1) = dRs_water_pure_dp * pow(10,-0.0840655 * S * pow(t,-0.285854));      // scf/STB/psia
    } else {
        Rs_WATER(0) = 0;      // scf/STB
        Rs_WATER(1) = 0;      // scf/STB/psia
    }

    return(Rs_WATER);
}


double COMPRESSIBILITY_WATER_MEEHAN(const std::string gas_saturated, double t, double p, double S, double Bg) {

    // t in R
    // p in psia
    // S salinity in weight percent TDS
    double Rsw = Rs_WATER_MEEHAN(gas_saturated,t,p,S)(0);
    double dRsw_dp = Rs_WATER_MEEHAN(gas_saturated,t,p,S)(1);                 // scf/STB/psia
    double Bw = B_WATER_MEEHAN(gas_saturated,t,p,S);                          // rb/STB
    // double dBw_dp = B_WATER(t,p)(1);                     // rb/STB/psia
    t = t - 459.67;                      // R to F
    double C0; double C1; double C2;
    double Cwf; double Cw; double Sc; double Cw1; double Cw2;
    double COMPRESSIBILITY_WATER_MEEHAN;
    C0 = 3.8546 - 0.000134 * p;
    C1 = -0.01052 + 4.77e-7 * p;
    C2 = 3.9267e-5 - 8.8e-10 * p;
    Cwf = 1e-6 * (C0 + C1 * t + C2 * t * t);
    Cw = Cwf * (1 + 8.9e-3 * Rsw);
    Sc = 1 + (-0.052 + 2.7e-4 * t - 1.14e-6 * t * t + 1.121e-9 * t * t * t) * pow(S,0.7);
    Cw1 = Cw * Sc;
    if (gas_saturated == "yes") {
        Cw2 = Bg / Bw * dRsw_dp;
    } else {
        Cw2 = 0;
    }
    COMPRESSIBILITY_WATER_MEEHAN = Cw1 + Cw2;
    return(COMPRESSIBILITY_WATER_MEEHAN);
}


double VISCOSITY_WATER_MEEHAN(double t, double p, double S) {

    // t in R
    // p in psia
    // S salinity in weight percent TDS
    t = t - 459.67;     // R to F
    double MUw;
    double MUw1;
    double A;
    double B;
    double f;
    double VISCOSITY_WATER_MEEHAN;

    A = -0.04518 + 0.009313 * S - 0.000393 * S * S;
    B = 70.634 + 0.09576 * S * S;
    f = 1 + 3.5e-12 * p * p * (t - 40);
    MUw1 = A + B / t;
    MUw = MUw1 * f;
    VISCOSITY_WATER_MEEHAN = MUw;                                  // cp
    return(VISCOSITY_WATER_MEEHAN);
}


arma::mat PVT_WATER_PROPERTIES_MEEHAN(const std::string gas_saturated, double t, arma::vec p, double S, double tpc, double ppc) {

    double psca = 14.696;              // Psia
    double tsc = (60 + 459.67);        // R
    double Bg;
    int lp = p.size();
    arma::mat results_table(lp,4);

    for (int i = 0; i < lp; i++) {
        results_table(i,0) = Rs_WATER_MEEHAN(gas_saturated,t,p(i),S)(0);                          // rb/STB
        results_table(i,1) = B_WATER_MEEHAN(gas_saturated,t,p(i),S);
        results_table(i,2) = DENSITY_WATER_MEEHAN(gas_saturated,t,p(i),S);
        Bg = B_GAS_DAK(t,p(i),tsc,psca,tpc,ppc)(1);                                                 // rb/STB
        results_table(i,3) = COMPRESSIBILITY_WATER_MEEHAN(gas_saturated,t,p(i),S,Bg);
    }
    return(results_table);
}


arma::vec VISC_WATER_PROPERTIES_MEEHAN(double t, arma::vec p, double S) {

    int lp = p.size();
    arma::vec results_table(lp);
    for (int i = 0; i < lp; i++) {
        results_table(i) = VISCOSITY_WATER_MEEHAN(t,p(i),S);
    }
    return(results_table);
}



// ************ SPIVEY CORRELATIONS ************
// [[Rcpp::export]]
arma::vec WATER_SPIVEY(const std::string gas_saturated, double t, double p, double S, double z_factor) {

    // "Estimating Density, Formation Volume Factor, Compressibility, Methane Solubility,
    // and Viscosity for Oilfield Brines at Temperatures From 0 to 275˚ C,
    // Pressures to 200 MPa, and Salinities to 5.7 mole/kg"
    // J.P. SPIVEY, W.D. MCCAIN, JR., Texas A&M University
    // JCPT, July 2004, Volume 43, No. 7, 52-61

    // t in R
    // p in psia
    // S salinity in weight percent TDS
    // 1 wt% = 10000 ppm
    // 1 ppm = 1 mg/lit

    arma::vec results(4);
    const double kpa_to_psi = 0.14503773800722;
    const double Mw_NaCl = 58.4428; // g/gmole
    const double Mw_CH4 = 16.043; // g/gmole

    double T = t / 1.8 - 273.15;         // C
    double T_sc = (60 - 32) / 1.8;       // C
    double P = p / kpa_to_psi / 1000;   // MPa
    double P_sc = 14.696 / kpa_to_psi / 1000;         // MPa
    double m = 1000. * (S / 100.) / (Mw_NaCl * (1 - S / 100.));  // gmole/Kg
    const double Tc = 647.096;
    const double Pc = 22.064;
    const double R = 8.314467;
    double T_K = T + 273.15;
    double tau = 1 - T_K / Tc;

    arma::vec DENSw = {-0.127213, 0.645486, 1.03265, -0.070291, 0.639589};
    arma::vec Ew = {4.221, -3.478, 6.221, 0.5182, -0.4405};
    arma::vec Fw = {-11.403, 29.932, 27.952, 0.20684, 0.3768};
    double DENSITYw_70 = (DENSw(0) * pow((T / 100.), 2.0) + DENSw(1) * (T / 100.) + DENSw(2)) / (DENSw(3) * pow((T / 100.), 2.0) + DENSw(4) * (T / 100.) + 1);
    double Ew_T = (Ew(0) * pow((T / 100.), 2.0) + Ew(1) * (T / 100.) + Ew(2)) / (Ew(3) * pow((T / 100.), 2.0) + Ew(4) * (T / 100.) + 1);
    double Fw_T = (Fw(0) * pow((T / 100.), 2.0) + Fw(1) * (T / 100.) + Fw(2)) / (Fw(3) * pow((T / 100.), 2.0) + Fw(4) * (T / 100.) + 1);
    // double I_W_70 = 1 / Ew_T * std::log(std::abs(Ew_T + Fw_T));
    // double I_W_P = 1 / Ew_T * std::log(std::abs(Ew_T * (P / 70.) + Fw_T));
    // double DENSITY_WATER = DENSITYw_70 * std::exp(I_W_P - I_W_70);
    // double COMPRESSIBILITY_WATER = (1 / 70.) / (Ew_T * (P / 70.) + Fw_T);
    double DENSITYw_70_sc = (DENSw(0) * pow((T_sc / 100.), 2.0) + DENSw(1) * (T_sc / 100.) + DENSw(2)) / (DENSw(3) * pow((T_sc / 100.), 2.0) + DENSw(4) * (T_sc / 100.) + 1);
    double Ew_T_sc = (Ew(0) * pow((T_sc / 100.), 2.0) + Ew(1) * (T_sc / 100.) + Ew(2)) / (Ew(3) * pow((T_sc / 100.), 2.0) + Ew(4) * (T_sc / 100.) + 1);
    double Fw_T_sc = (Fw(0) * pow((T_sc / 100.), 2.0) + Fw(1) * (T_sc / 100.) + Fw(2)) / (Fw(3) * pow((T_sc / 100.), 2.0) + Fw(4) * (T_sc / 100.) + 1);
    // double I_W_70_sc = 1 / Ew_T_sc * std::log(std::abs(Ew_T_sc + Fw_T_sc));
    // double I_W_P_sc = 1 / Ew_T_sc * std::log(std::abs(Ew_T_sc * (P_sc / 70.) + Fw_T_sc));
    // double DENSITY_WATER_sc = DENSITYw_70_sc * std::exp(I_W_P_sc - I_W_70_sc);

    arma::vec Dm_2 = {-1.1149e-4, 1.7105e-4, -4.3766e-4, 0.0, 0.0};
    arma::vec Dm_3_2 = {-8.878e-4, -1.388e-4, -2.96318e-3, 0.0, 0.51103};
    arma::vec Dm_1 = {2.1466e-3, 1.2427e-2, 4.2648e-2, -8.1009e-2, 0.525417};
    arma::vec Dm_1_2 = {2.356e-4, -3.636e-4, -2.278e-4, 0.0, 0.0};
    double Dm_2_T = (Dm_2(0) * pow((T / 100.), 2.0) + Dm_2(1) * (T / 100.) + Dm_2(2)) / (Dm_2(3) * pow((T / 100.), 2.0) + Dm_2(4) * (T / 100.) + 1);
    double Dm_3_2_T = (Dm_3_2(0) * pow((T / 100.), 2.0) + Dm_3_2(1) * (T / 100.) + Dm_3_2(2)) / (Dm_3_2(3) * pow((T / 100.), 2.0) + Dm_3_2(4) * (T / 100.) + 1);
    double Dm_1_T = (Dm_1(0) * pow((T / 100.), 2.0) + Dm_1(1) * (T / 100.) + Dm_1(2)) / (Dm_1(3) * pow((T / 100.), 2.0) + Dm_1(4) * (T / 100.) + 1);
    double Dm_1_2_T = (Dm_1_2(0) * pow((T / 100.), 2.0) + Dm_1_2(1) * (T / 100.) + Dm_1_2(2)) / (Dm_1_2(3) * pow((T / 100.), 2.0) + Dm_1_2(4) * (T / 100.) + 1);
    double DENSITYb_70 = DENSITYw_70 + Dm_2_T * pow(m, 2.0) + Dm_3_2_T * pow(m, 1.5) + Dm_1_T * pow(m, 1.0) + Dm_1_2_T * pow(m, 0.5);
    double Dm_2_T_sc = (Dm_2(0) * pow((T_sc / 100.), 2.0) + Dm_2(1) * (T_sc / 100.) + Dm_2(2)) / (Dm_2(3) * pow((T_sc / 100.), 2.0) + Dm_2(4) * (T_sc / 100.) + 1);
    double Dm_3_2_T_sc = (Dm_3_2(0) * pow((T_sc / 100.), 2.0) + Dm_3_2(1) * (T_sc / 100.) + Dm_3_2(2)) / (Dm_3_2(3) * pow((T_sc / 100.), 2.0) + Dm_3_2(4) * (T_sc / 100.) + 1);
    double Dm_1_T_sc = (Dm_1(0) * pow((T_sc / 100.), 2.0) + Dm_1(1) * (T_sc / 100.) + Dm_1(2)) / (Dm_1(3) * pow((T_sc / 100.), 2.0) + Dm_1(4) * (T_sc / 100.) + 1);
    double Dm_1_2_T_sc = (Dm_1_2(0) * pow((T_sc / 100.), 2.0) + Dm_1_2(1) * (T_sc / 100.) + Dm_1_2(2)) / (Dm_1_2(3) * pow((T_sc / 100.), 2.0) + Dm_1_2(4) * (T_sc / 100.) + 1);
    double DENSITYb_70_sc = DENSITYw_70_sc + Dm_2_T_sc * pow(m, 2.0) + Dm_3_2_T_sc * pow(m, 1.5) + Dm_1_T_sc * pow(m, 1.0) + Dm_1_2_T_sc * pow(m, 0.5);

    arma::vec Em = {0.0, 0.0, 0.1249, 0.0, 0.0};
    arma::vec Fm_3_2 = {-0.617, -0.747, -0.4339, 0.0, 10.26};
    arma::vec Fm_1 = {0.0, 9.917, 5.1128, 0.0, 3.892};
    arma::vec Fm_1_2 = {0.0365, -0.0369, 0.0, 0.0, 0.0};

    double Em_T = (Em(0) * pow((T / 100.), 2.0) + Em(1) * (T / 100.) + Em(2)) / (Em(3) * pow((T / 100.), 2.0) + Em(4) * (T / 100.) + 1);
    double Fm_3_2_T = (Fm_3_2(0) * pow((T / 100.), 2.0) + Fm_3_2(1) * (T / 100.) + Fm_3_2(2)) / (Fm_3_2(3) * pow((T / 100.), 2.0) + Fm_3_2(4) * (T / 100.) + 1);
    double Fm_1_T = (Fm_1(0) * pow((T / 100.), 2.0) + Fm_1(1) * (T / 100.) + Fm_1(2)) / (Fm_1(3) * pow((T / 100.), 2.0) + Fm_1(4) * (T / 100.) + 1);
    double Fm_1_2_T = (Fm_1_2(0) * pow((T / 100.), 2.0) + Fm_1_2(1) * (T / 100.) + Fm_1_2(2)) / (Fm_1_2(3) * pow((T / 100.), 2.0) + Fm_1_2(4) * (T / 100.) + 1);
    double Ebm = Ew_T + Em_T * m;
    double Fbm = Fw_T + Fm_3_2_T * pow(m, 1.5) + Fm_1_T * pow(m, 1.0) + Fm_1_2_T * pow(m, 0.5);
    double I_B_70 = 1 / Ebm * std::log(std::abs(Ebm + Fbm));
    double I_B_P = 1 / Ebm * std::log(std::abs(Ebm * (P / 70.) + Fbm));
    double DENSITY_BRINE = DENSITYb_70 * std::exp(I_B_P - I_B_70);
    double V_BRINE = 1 / DENSITY_BRINE;
    double COMPRESSIBILITY_BRINE = (1 / 70.) / (Ebm * (P / 70.) + Fbm);
    double GAS_SOLUBILITY;

    double Em_T_sc = (Em(0) * pow((T_sc / 100.), 2.0) + Em(1) * (T_sc / 100.) + Em(2)) / (Em(3) * pow((T_sc / 100.), 2.0) + Em(4) * (T_sc / 100.) + 1);
    double Fm_3_2_T_sc = (Fm_3_2(0) * pow((T_sc / 100.), 2.0) + Fm_3_2(1) * (T_sc / 100.) + Fm_3_2(2)) / (Fm_3_2(3) * pow((T_sc / 100.), 2.0) + Fm_3_2(4) * (T_sc / 100.) + 1);
    double Fm_1_T_sc = (Fm_1(0) * pow((T_sc / 100.), 2.0) + Fm_1(1) * (T_sc / 100.) + Fm_1(2)) / (Fm_1(3) * pow((T_sc / 100.), 2.0) + Fm_1(4) * (T_sc / 100.) + 1);
    double Fm_1_2_T_sc = (Fm_1_2(0) * pow((T_sc / 100.), 2.0) + Fm_1_2(1) * (T_sc / 100.) + Fm_1_2(2)) / (Fm_1_2(3) * pow((T_sc / 100.), 2.0) + Fm_1_2(4) * (T_sc / 100.) + 1);
    double Ebm_sc = Ew_T_sc + Em_T_sc * m;
    double Fbm_sc = Fw_T_sc + Fm_3_2_T_sc * pow(m, 1.5) + Fm_1_T_sc * pow(m, 1.0) + Fm_1_2_T_sc * pow(m, 0.5);
    double I_B_70_sc = 1 / Ebm_sc * std::log(std::abs(Ebm_sc + Fbm_sc));
    double I_B_P_sc = 1 / Ebm_sc * std::log(std::abs(Ebm_sc * (P_sc / 70.) + Fbm_sc));
    double DENSITY_BRINE_sc = DENSITYb_70_sc * std::exp(I_B_P_sc - I_B_70_sc);
    double V_BRINE_sc = 1 / DENSITY_BRINE_sc;
    double B_BRINE = V_BRINE / V_BRINE_sc;

    GAS_SOLUBILITY = 0.0;
    results(0) = GAS_SOLUBILITY;
    results(1) = B_BRINE;
    results(2) = DENSITY_BRINE;
    results(3) = COMPRESSIBILITY_BRINE;

    if (gas_saturated == "yes") {
        arma::vec a = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502};
        double P_SAT = Pc * std::exp((Tc / T_K) * (a(0) * tau + a(1) * pow(tau, 1.5) + a(2) * pow(tau, 3.0) + a(3) * pow(tau, 3.5) + a(4) * pow(tau, 4.0) + a(5) * pow(tau, 7.5)));

        arma::vec A = {0.0, -0.004462, -0.06763, 0.0, 0.0};
        arma::vec B = {-0.03602, 0.18917, 0.97242, 0.0, 0.0};
        arma::vec C = {0.6855, -3.1992, -3.7968, 0.07711, 0.2229};
        double A_T = (A(0) * pow((T / 100.), 2.0) + A(1) * (T / 100.) + A(2)) / (A(3) * pow((T / 100.), 2.0) + A(4) * (T / 100.) + 1);
        double B_T = (B(0) * pow((T / 100.), 2.0) + B(1) * (T / 100.) + B(2)) / (B(3) * pow((T / 100.), 2.0) + B(4) * (T / 100.) + 1);
        double C_T = (C(0) * pow((T / 100.), 2.0) + C(1) * (T / 100.) + C(2)) / (C(3) * pow((T / 100.), 2.0) + C(4) * (T / 100.) + 1);
        double m_C4_w = std::exp(A_T * pow(std::log(P - P_SAT), 2.0) + B_T * std::log(P - P_SAT) + C_T);

        arma::vec c_lambda = {-0.80898, 1.0827e-3, 183.85, 3.924e-4, -1.97e-6};
        double lambda_C4b_TP = c_lambda(0) + c_lambda(1) * T_K + c_lambda(2) / T_K + c_lambda(3) * P + c_lambda(4) * pow(P, 2.0);
        double kici_C4b = -3.89e-3;
        double m_C4_b = m_C4_w *std::exp(-2. * lambda_C4b_TP * m - kici_C4b * pow(m, 2.0));

        arma::vec c_mu = {8.3143711, -7.2772168e-4, 2.1489858e3, -1.4019672e-5, -6.6743449e5, 7.698589e-2, -5.0253331e-5, -30.092013, 4.8468502e3, 0.0};
        double V_m_C4_b = R * T_K * (c_mu(5) + c_mu(6) * T_K + c_mu(7) / T_K + c_mu(8) / pow(T_K, 2.0) + 2 * m * (c_lambda(3) + 2 * c_lambda(4) * P) + m * m * 0);

        double DENSITY_BRINE_SAT = (1000. + m * Mw_NaCl + m_C4_b * Mw_CH4) / ((1000. + m * Mw_NaCl) * V_BRINE + m_C4_b * V_m_C4_b);
        // double V_BRINE_SAT = 1 / DENSITY_BRINE_SAT;

        double V_m_C4_g = z_factor * R * T_K / P;
        double part_I = (1000. + m * Mw_NaCl) * (-V_BRINE * COMPRESSIBILITY_BRINE);
        double part_II = m_C4_b * (R * T_K * (0.0 + 2 * m * 2.0 * c_lambda(4) + m * m * 0.0));
        double part_III = m_C4_b * ((2 * A_T * std::log(P - P_SAT) + B_T) / (P - P_SAT) - 2 * (c_lambda(3) + 2 * c_lambda(4) * P) * m) * (V_m_C4_b - V_m_C4_g);
        double COMPRESSIBILITY_BRINE_SAT = -(part_I + part_II + part_III) / ((1000. + m * Mw_NaCl) * V_BRINE + m_C4_b * V_m_C4_b);
        double B_BRINE_SAT = ((1000. + m * Mw_NaCl) * V_BRINE + m_C4_b * V_m_C4_b) / ((1000. + m * Mw_NaCl) * V_BRINE_sc);
        GAS_SOLUBILITY = m_C4_b * (1 * R * (T_sc + 273.15) / P_sc) / (1000.0 + m * Mw_NaCl) / V_BRINE_sc;
        results(0) = GAS_SOLUBILITY;
        results(1) = B_BRINE_SAT;
        results(2) = DENSITY_BRINE_SAT;
        results(3) = COMPRESSIBILITY_BRINE_SAT;
    }
    return results;
}


double VISCOSITY_WATER_SPIVEY(double t, double p, double S) {

    // t in R
    // p in psia
    // S salinity in weight percent TDS
    // 1 wt% = 10000 ppm
    // 1 ppm = 1 mg/lit

    double VISCOSITY_WATER_SPIVEY;
    const double kpa_to_psi = 0.14503773800722;
    const double Mw_NaCl = 58.4428; // g/gmole

    double T = t / 1.8 - 273.15;         // C
    double P = p / kpa_to_psi / 1000;   // MPa
    double m = 1000. * (S / 100.) / (Mw_NaCl * (1 - S / 100.));  // gmole/Kg
    double T_K = T + 273.15;


    arma::vec DENSw = {-0.127213, 0.645486, 1.03265, -0.070291, 0.639589};
    arma::vec Ew = {4.221, -3.478, 6.221, 0.5182, -0.4405};
    arma::vec Fw = {-11.403, 29.932, 27.952, 0.20684, 0.3768};
    double DENSITYw_70 = (DENSw(0) * pow((T / 100.), 2.0) + DENSw(1) * (T / 100.) + DENSw(2)) / (DENSw(3) * pow((T / 100.), 2.0) + DENSw(4) * (T / 100.) + 1);
    double Ew_T = (Ew(0) * pow((T / 100.), 2.0) + Ew(1) * (T / 100.) + Ew(2)) / (Ew(3) * pow((T / 100.), 2.0) + Ew(4) * (T / 100.) + 1);
    double Fw_T = (Fw(0) * pow((T / 100.), 2.0) + Fw(1) * (T / 100.) + Fw(2)) / (Fw(3) * pow((T / 100.), 2.0) + Fw(4) * (T / 100.) + 1);
    double I_W_70 = 1 / Ew_T * std::log(std::abs(Ew_T + Fw_T));
    double I_W_P = 1 / Ew_T * std::log(std::abs(Ew_T * (P / 70.) + Fw_T));
    double DENSITY_WATER = DENSITYw_70 * std::exp(I_W_P - I_W_70);

    arma::vec d = {0.28853170e7, -0.11072577e5, -0.90834095e1, 0.30925651e-1, -0.27407100e-4, -0.19283851e7, 0.56216046e4, 0.13827250e2, -0.47609523e-1, 0.35545041e-4};
    arma::vec aa = {-0.21319213, 0.13651589e-2, -0.12191756e-5};
    arma::vec bb = {0.69161945e-1, -0.27292263e-3, 0.20852448e-6};
    arma::vec cc = {-0.25988855e-2, 0.77989227e-5};

    double sum_mu_1 = 0.0;
    for (int i = 0; i < 5; i++) {
        sum_mu_1 += d(i) * std::pow(T_K, i - 2);
    }
    double sum_mu_2 = 0.0;
    for (int i = 5; i < 10; i++) {
        sum_mu_2 += d(i) * std::pow(T_K, i - 7);
    }
    double MU_WATER = std::exp(sum_mu_1 + DENSITY_WATER * sum_mu_2);
    double AA = aa(0) + aa(1) * T_K + aa(2) * std::pow(T_K, 2.0);
    double BB = bb(0) + bb(1) * T_K + bb(2) * std::pow(T_K, 2.0);
    double CC = cc(0) + cc(1) * T_K;
    double MU_RATIO = std::exp(AA * m + BB * m * m + CC * m * m * m);
    double MU_BRINE = MU_WATER * MU_RATIO;

    VISCOSITY_WATER_SPIVEY = MU_BRINE * 1000;               // cp
    return(VISCOSITY_WATER_SPIVEY);
}


arma::mat PVT_WATER_PROPERTIES_SPIVEY(const std::string gas_saturated, double t, arma::vec p, double S, double tpc, double ppc) {

    double kpa_to_psi = 0.14503773800722;
    double m3_to_scf = 35.314667;
    double m3_to_bbl = 6.289814;
    double kg_m3_to_lb_cuft = 0.06242796;
    double z_factor;
    arma::vec temp(5);
    int lp = p.size();
    arma::mat results_table(lp,4);

    for (int i = 0; i < lp; i++) {
        z_factor = Z_FACTOR_DAK(t, p(i), tpc,ppc);                   // rb/STB
        temp = WATER_SPIVEY(gas_saturated, t, p(i), S, z_factor);
        results_table(i,0) = temp(0) * (m3_to_scf / m3_to_bbl);      // scf/STB
        results_table(i,1) = temp(1);
        results_table(i,2) = temp(2) * 1000 * kg_m3_to_lb_cuft;      // lb/cuft
        results_table(i,3) = temp(3) / 1000 / kpa_to_psi;            // 1/psi
    }
    return(results_table);
}

arma::vec VISC_WATER_PROPERTIES_SPIVEY(double t, arma::vec p, double S) {

    int lp = p.size();
    arma::vec results_table(lp);
    for (int i = 0; i < lp; i++) {
        results_table(i) = VISCOSITY_WATER_SPIVEY(t,p(i),S);
    }
    return(results_table);
}

