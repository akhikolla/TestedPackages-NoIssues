// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List upalpha(Rcpp::List hgmat, double alpha, Rcpp::String diffmeth, double h){

   
  // option for updating alpha
  
  if(diffmeth == "analytic"){

    // inputs
    arma::sp_mat hmat = Rcpp::as<arma::sp_mat>(hgmat["hmat"]);
    arma::sp_mat ghmat = Rcpp::as<arma::sp_mat>(hgmat["ghmat"]);
    arma::sp_mat gghmat = Rcpp::as<arma::sp_mat>(hgmat["gghmat"]);
    arma::sp_mat gmat = Rcpp::as<arma::sp_mat>(hgmat["gmat"]);
    arma::sp_mat ggmat = Rcpp::as<arma::sp_mat>(hgmat["ggmat"]);
    arma::sp_mat gggmat = Rcpp::as<arma::sp_mat>(hgmat["gggmat"]);

    // calculate new variables
    // could leave as sparse and use ARPACK
    arma::mat dhmat(hmat), dgmat(gmat);
    arma::mat dghmat(ghmat), dggmat(ggmat);
    unsigned int ncoeff = dhmat.n_cols;


    // block diagonal locations
    arma::mat clocblock = arma::linspace(0, ncoeff * ncoeff -1, ncoeff * ncoeff);
    arma::mat vrow = arma::repmat(clocblock,1,ncoeff);
    clocblock.reshape(ncoeff, ncoeff);
    arma::mat vcol = arma::repmat(clocblock.t(),1,ncoeff);
    arma::vec loccol = arma::vectorise(vcol.t());
    arma::vec locrow = arma::vectorise(vrow.t());
    arma::mat locats = arma::join_rows(locrow,loccol);
    arma::umat ulocats = arma::conv_to<arma::umat>::from(locats.t());


    // eigen analysis
    arma::vec xheigval, xgeigval;
    arma::mat xheigvec, xgeigvec;
    arma::eig_sym(xheigval, xheigvec, dhmat);
    arma::eig_sym(xgeigval, xgeigvec, dgmat);
    arma::vec heigval = arma::flipud(xheigval);
    arma::mat heigvec = arma::fliplr(xheigvec);
    arma::vec geigval = arma::flipud(xgeigval);
    arma::mat geigvec = arma::fliplr(xgeigvec);


    // construct sparse matrix
    arma::mat cprodh = heigvec.t() * dghmat * heigvec;
    double gldethmat = arma::sum(cprodh.diag() / heigval);
    arma::mat cprodg = geigvec.t() * dggmat * geigvec;
    double gldetgmat = arma::sum(cprodg.diag() / geigval);
    // hmat
    arma::mat heigdenom = arma::repmat(heigval,1,ncoeff) - arma::repmat(heigval.t(),ncoeff,1);
    arma::mat ex1heigvec = arma::repmat(heigvec, ncoeff, 1);
    ex1heigvec.reshape(ncoeff, ncoeff * ncoeff);
    arma::mat exheigvec = heigvec;
    exheigvec.reshape(ncoeff * ncoeff,1);
    arma::mat ex2heigvec = arma::repmat(exheigvec, 1, ncoeff);
    arma::mat pheigvec = ex1heigvec.t() % ex2heigvec;
    arma::vec vfhnum = arma::vectorise(pheigvec.t());
    arma::sp_mat fhnum(ulocats, vfhnum);
    // gmat
    arma::mat geigdenom = arma::repmat(geigval,1,ncoeff) - arma::repmat(geigval.t(),ncoeff,1);
    arma::mat ex1geigvec = arma::repmat(geigvec, ncoeff, 1);
    ex1geigvec.reshape(ncoeff, ncoeff * ncoeff);
    arma::mat exgeigvec = geigvec;
    exgeigvec.reshape(ncoeff * ncoeff,1);
    arma::mat ex2geigvec = arma::repmat(exgeigvec, 1, ncoeff);
    arma::mat pgeigvec = ex1geigvec.t() % ex2geigvec;
    arma::vec vfgnum = arma::vectorise(pgeigvec.t());
    arma::sp_mat fgnum(ulocats, vfgnum);


    // structures to store intermediate results
    arma::sp_mat intmat = arma::repmat(arma::speye<arma::sp_mat>(ncoeff,ncoeff),1,ncoeff);
    // hmat
    arma::mat kheigdenom(ncoeff, ncoeff);
    arma::uvec qh;
    arma::vec vkheigdenom(ncoeff * ncoeff), ddheig(ncoeff);
    arma::sp_mat fhmat(ncoeff,ncoeff), dheig(ncoeff,ncoeff), bhmat(ncoeff*ncoeff, ncoeff*ncoeff); 
    // gmat
    arma::mat kgeigdenom(ncoeff, ncoeff);
    arma::uvec qg;
    arma::vec vkgeigdenom(ncoeff * ncoeff), ddgeig(ncoeff);
    arma::sp_mat fgmat(ncoeff,ncoeff), dgeig(ncoeff,ncoeff), bgmat(ncoeff*ncoeff, ncoeff*ncoeff);
   


    for(unsigned int k=0; k<ncoeff; k++){

     // hmat
     kheigdenom = arma::repmat(heigdenom.col(k), 1, ncoeff);
     vkheigdenom = arma::vectorise(kheigdenom.t());
     qh = arma::find(vkheigdenom);
     vkheigdenom.elem(qh) = 1 / vkheigdenom.elem(qh);
     bhmat.diag() = vkheigdenom;
     fhmat = intmat * (bhmat * fhnum) * intmat.t();
     dheig = gghmat - 2 * ghmat * fhmat * ghmat;
     ddheig(k) =  arma::as_scalar(heigvec.col(k).t() * dheig * heigvec.col(k));

     // gmat
     kgeigdenom = arma::repmat(geigdenom.col(k), 1, ncoeff);
     vkgeigdenom = arma::vectorise(kgeigdenom.t());
     qg = arma::find(vkgeigdenom);
     vkgeigdenom.elem(qg) = 1 / vkgeigdenom.elem(qg);
     bgmat.diag() = vkgeigdenom;
     fgmat = intmat * (bgmat * fgnum) * intmat.t();
     dgeig = gggmat - 2 * ggmat * fgmat * ggmat;
     ddgeig(k) =  arma::as_scalar(geigvec.col(k).t() * dgeig * geigvec.col(k));


    }

    // calculate gvb and ggvb
    double ggldethmat = arma::as_scalar(sum(ddheig / heigval)) -
                         arma::as_scalar(arma::sum(arma::pow(cprodh.diag() / heigval, 2.0)));
    double ggldetgmat = arma::as_scalar(arma::sum(ddgeig / geigval)) -
                         arma::as_scalar(arma::sum(arma::pow(cprodg.diag() / geigval, 2.0)));
    double gvb = gldetgmat - 2 * gldethmat;
    double ggvb = ggldetgmat - 2 * ggldethmat;

    // update alpha
    double phi = std::log(alpha) - std::log(1 - alpha);
    double gphi = std::exp(phi) / std::pow(1 + std::exp(phi), 2.0);
    double ggphi = (std::exp(phi)*(1 - std::exp(phi))) / std::pow(1 + std::exp(phi), 3.0);
    double gvbphi = gphi * gvb;
    double ggvbphi = ggphi * gvb + std::pow(gphi, 2.0) * ggvb;
    double phiinc = gvbphi / std::abs(ggvbphi);
    double nphi, nalpha;
    if (std::abs(phiinc) > 1){
          if (phiinc > 1){nphi = phi - 1;} else {nphi = phi + 1;}
       }
         else {
       nphi = phi - phiinc;
     }
    nalpha = std::exp(nphi) / (1 + std::exp(nphi));
    if (nalpha > 0.95){nalpha = 0.95;}
    if (nalpha < 0.05){nalpha = 0.05;}

    // output
    return Rcpp::List::create(Rcpp::Named("gvb")=gvbphi,
                              Rcpp::Named("ggvb")=ggvbphi,
                              Rcpp::Named("alpha")=nalpha);


   } else if(diffmeth == "numeric"){

    // inputs 
    arma::sp_mat hmat = Rcpp::as<arma::sp_mat>(hgmat["hmat"]);
    arma::sp_mat hmatl = Rcpp::as<arma::sp_mat>(hgmat["lhmat"]);
    arma::sp_mat hmatu = Rcpp::as<arma::sp_mat>(hgmat["uhmat"]);
    arma::sp_mat gmat = Rcpp::as<arma::sp_mat>(hgmat["gmat"]);
    arma::sp_mat gmatl = Rcpp::as<arma::sp_mat>(hgmat["lgmat"]);
    arma::sp_mat gmatu = Rcpp::as<arma::sp_mat>(hgmat["ugmat"]);


    // variance matrices
    arma::mat xgmat(gmat), xgmatl(gmatl), xgmatu(gmatu), 
              xhmat(hmat), xhmatl(hmatl), xhmatu(hmatu);
    double ldetvcov = std::log(arma::det(xgmat)) - 2 * std::log(arma::det(xhmat));
    double ldetvcovl = std::log(arma::det(xgmatl)) - 2 * std::log(arma::det(xhmatl));
    double ldetvcovu = std::log(arma::det(xgmatu)) - 2 * std::log(arma::det(xhmatu));

    // gvb, ggvb and alpha
    double gvbphi = (ldetvcovu - ldetvcovl)/(2 * arma::as_scalar(h));
    double ggvbphi = (ldetvcovu - 2 * ldetvcov + ldetvcovl)/(std::pow(arma::as_scalar(h), 2.0));
    double phi = std::log(alpha) - std::log(1 - alpha);
    double phiinc = gvbphi / std::abs(ggvbphi);
    double nphi, nalpha;
    if (std::abs(phiinc) > 1){
          if (phiinc > 1){nphi = phi - 1;} else {nphi = phi + 1;}
       }
         else {
       nphi = phi - phiinc;
     }
    nalpha = std::exp(nphi) / (1 + std::exp(nphi));
    if (nalpha > 0.95){nalpha = 0.95;}
    if (nalpha < 0.05){nalpha = 0.05;}

    // output
    return Rcpp::List::create(Rcpp::Named("gvb")=gvbphi,
                              Rcpp::Named("ggvb")=ggvbphi,
                              Rcpp::Named("alpha")=nalpha);

   } else {

   // output
   return 0;

   }

}


