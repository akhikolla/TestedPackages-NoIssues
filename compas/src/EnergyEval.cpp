#include <Rcpp.h>
using namespace Rcpp;

#include <stdio.h>
#include <iostream>

// double edis(Eigen::VectorXd x, Eigen::VectorXd y){
//   const Eigen::VectorXd diff = x - y;
//   return sqrt(diff.dot(diff));
// }

// [[Rcpp::export]]
double dfeval(DataFrame chain, NumericVector atypes, NumericMatrix& etable){
    int nrows = chain.nrows();
    double energy = 0;
    //StringVector resids = chain["resid"];
    //StringVector atoms = chain["elety"];
    NumericVector resno = chain["resno"];
    NumericVector x = chain["x"];
    NumericVector y = chain["y"];
    NumericVector z = chain["z"];
    for (int i=0; i<nrows; i++){
        //std::string resid1 = as<std::string>(resids[i]); 
        //std::string atom1 = as<std::string>(atoms[i]);
        //Eigen::VectorXd coord1(3);
        //coord1 << x[i], y[i], z[i];
        for (int j=i+1; j<nrows; j++){
          if (resno[i] == resno[j])
            continue;
          //std::string resid2 = as<std::string>(resids[j]); 
          //std::string atom2 = as<std::string>(atoms[j]); 
          //Eigen::VectorXd coord2(3);
          //coord2 << x[j], y[j], z[j];
          int adis = sqrt((x[i]-x[j])*(x[i]-x[j]) +(y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]));
          int bin = floor(2.0 * adis);
          //bin = bin>29 ? 29: bin;
          //energy += lookup(chart, resid1, atom1, resid2, atom2, bin);
          if (bin < 30)
            energy += etable(atypes[i]*167 + atypes[j],bin);
        }
    }
    return energy;
}

