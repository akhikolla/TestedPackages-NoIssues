// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List hgmat(List mod, List smat, List cmat, Rcpp::List X, Rcpp::String modtype,
                            Rcpp::String diffmeth){

    // inputs
    arma::sp_mat Xmat = Rcpp::as<arma::sp_mat>(X["design"]);
    arma::rowvec linpred = Rcpp::as<arma::rowvec>(mod["linear.predictors"]);
    arma::rowvec fitted = Rcpp::as<arma::rowvec>(mod["fitted.values"]);
    arma::rowvec y = Rcpp::as<arma::rowvec>(mod["y"]);
    arma::rowvec residuals = y - fitted;
    arma::sp_mat icmat = Rcpp::as<arma::sp_mat>(cmat["icmat"]);
    arma::sp_mat gicmat, ggicmat;
    if(diffmeth == "analytic"){
     gicmat = Rcpp::as<arma::sp_mat>(cmat["gicmat"]);
     ggicmat = Rcpp::as<arma::sp_mat>(cmat["ggicmat"]);
    } else if(diffmeth == "numeric"){
     gicmat = Rcpp::as<arma::sp_mat>(cmat["licmat"]);
     ggicmat = Rcpp::as<arma::sp_mat>(cmat["uicmat"]);
    }
    arma::sp_mat ismat = Rcpp::as<arma::sp_mat>(smat["ismat"]);
    arma::sp_mat ssmat = Rcpp::as<arma::sp_mat>(smat["smat"]);
    arma::rowvec id;
    if (modtype == "glm"){
      Rcpp::List datalist = Rcpp::as<List>(mod["data"]);
      id = Rcpp::as<arma::rowvec>(datalist["subjects"]);
      } else if (modtype == "gee") {
      id = Rcpp::as<arma::rowvec>(mod["id"]);
     }  

    // calculate new variables
    unsigned int maxid = arma::max(id);
    unsigned int ntimes = icmat.n_rows;
    unsigned int cats = ismat.n_rows + 1;
    unsigned int nid = ntimes * (cats - 1);
    arma::mat mismat(ismat);
    arma::mat micmat(icmat), mgicmat(gicmat), mggicmat(ggicmat);
    arma::rowvec varmat = arma::sqrt(fitted % (1 - fitted));
    arma::sp_mat dmat = arma::zeros<arma::sp_mat>(nid * maxid, nid * maxid);
    arma::sp_mat vmat = arma::zeros<arma::sp_mat>(nid * maxid, nid * maxid);
    vmat.diag() = 1 / varmat;
    dmat.diag() = arma::exp(linpred) / arma::pow(1 + arma::exp(linpred), 2.0);
    arma::sp_mat ddmat = Xmat.t() * dmat;
     
    // block diagonal locations
    arma::mat clocblock = arma::linspace(0, nid * maxid-1, nid * maxid);
    arma::mat vrow = arma::repmat(clocblock,1,nid);
    clocblock.reshape(nid, maxid);
    arma::mat vcol = arma::repmat(clocblock.t(), 1, nid);
    arma::vec loccol = arma::vectorise(vcol.t());
    arma::vec locrow = arma::vectorise(vrow.t());
    arma::mat locats = arma::join_rows(locrow,loccol);
    arma::umat ulocats = arma::conv_to<arma::umat>::from(locats.t());
   
    // icormat
    arma::mat icor = arma::kron(micmat, mismat);
    arma::mat ficor = arma::repmat(icor,1, maxid);
    arma::vec vicor = arma::vectorise(ficor);
    arma::sp_mat icormat(ulocats, vicor);

    // gicormat and ggicormat
    arma::mat gicor = arma::kron(mgicmat, mismat);
    arma::mat fgicor = arma::repmat(gicor, 1, maxid);
    arma::vec vgicor = arma::vectorise(fgicor);
    arma::sp_mat gicormat(ulocats, vgicor);
    arma::mat ggicor = arma::kron(mggicmat, mismat);
    arma::mat fggicor = arma::repmat(ggicor, 1, maxid);
    arma::vec vggicor = arma::vectorise(fggicor);
    arma::sp_mat ggicormat(ulocats, vggicor);
  
    // raw residual matrix
    arma::mat qresmat(1,nid * maxid);
    qresmat.submat(0, 0, 0, nid * maxid - 1) = residuals;
    arma::mat resnid = arma::repmat(qresmat.t(), 1, nid);
    qresmat.reshape(nid, maxid);
    arma::mat nqresmat = arma::repmat(qresmat, nid, 1);
    nqresmat.reshape(nid, nid*maxid);
    arma::mat qresprod = resnid.t() % nqresmat;
    arma::vec vresprod = arma::vectorise(qresprod);
    arma::sp_mat qrawres(ulocats, vresprod);

    // between category covariance blocks
    arma::mat qssmat(ssmat);
    qssmat.reshape((cats - 1) * (cats - 1) ,1);
    arma::mat fssmat = arma::repmat(qssmat, ntimes * maxid, 1);
    arma::vec vssmat = arma::vectorise(fssmat);
    arma::mat fvarmat1 = arma::repmat(varmat, (cats - 1), 1);
    arma::vec varmat1 = arma::vectorise(fvarmat1);
    arma::mat nvarmat = arma::conv_to<arma::mat>::from(varmat);   
    nvarmat.reshape((cats - 1), ntimes * maxid);
    arma::mat fvarmat2 = arma::repmat(nvarmat, (cats - 1), 1);
    arma::vec varmat2 = arma::vectorise(fvarmat2);
    arma::vec insvcov = varmat1 % vssmat % varmat2;

    // between category covariance block locations
    arma::mat vcovblock = arma::linspace(0, nid * maxid - 1, nid * maxid);
    arma::mat rvcovblock = arma::repmat(vcovblock.t(),(cats - 1), 1);
    vcovblock.reshape((cats-1), ntimes * maxid);
    arma::mat cvcovblock = arma::repmat(vcovblock.t(), 1, (cats - 1));
    arma::vec rowvcov = arma::vectorise(rvcovblock);
    arma::vec colvcov = arma::vectorise(cvcovblock.t());
    arma::mat blocats = arma::join_rows(rowvcov,colvcov);
    arma::umat ublocats = arma::conv_to<arma::umat>::from(blocats.t());
    
    // replace covariance blocks
    arma::sp_mat selblock(ublocats, arma::ones<arma::vec>(insvcov.n_elem));
    arma::sp_mat repvcovblock = selblock % qrawres;
    arma::vec resrep = arma::nonzeros(repvcovblock);
    arma::umat alllocats = arma::join_rows(ulocats, arma::join_rows(ublocats, ublocats));
    arma::vec allvalues = arma::join_cols(vresprod, arma::join_cols(-resrep, insvcov));
    arma::sp_mat qfullresmat(true, alllocats, allvalues, maxid *nid, maxid * nid, true, true);
    
    // hmat, ghmat and gghmat and gmat, ggmat and gggmat
     arma::sp_mat gmat(Xmat.n_cols, Xmat.n_cols), ggmat(Xmat.n_cols, Xmat.n_cols),
                 gggmat(Xmat.n_cols, Xmat.n_cols), hmat(Xmat.n_cols, Xmat.n_cols),
                 ghmat(Xmat.n_cols, Xmat.n_cols), gghmat(Xmat.n_cols, Xmat.n_cols);

    if(diffmeth == "analytic"){

     // calculate hmat, ghmat and gghmat
     arma::sp_mat vcovmat = vmat * icormat * vmat;
     arma::sp_mat gvcovmat = vmat * gicormat * vmat;
     arma::sp_mat ggvcovmat = vmat * ggicormat * vmat;  
     hmat = ddmat * vcovmat * ddmat.t();
     ghmat = ddmat * gvcovmat * ddmat.t();
     gghmat = ddmat * ggvcovmat * ddmat.t();

     // calculate gmat, ggmat and gggmat
     gmat = ddmat * vcovmat * qfullresmat * vcovmat * ddmat.t();
     ggmat = ddmat * gvcovmat * qfullresmat * vcovmat * ddmat.t() + 
                             ddmat * vcovmat * qfullresmat * gvcovmat * ddmat.t();
     gggmat = ddmat * ggvcovmat * qfullresmat * vcovmat * ddmat.t() +
              2 * ddmat * gvcovmat * qfullresmat * gvcovmat * ddmat.t() +
              ddmat * vcovmat * qfullresmat * ggvcovmat * ddmat.t();
   
    } else if(diffmeth == "numeric"){

     // calculate matrices for alpha updating
     arma::sp_mat vcovmat = vmat * icormat * vmat;
     arma::sp_mat gvcovmat = vmat * gicormat * vmat;
     arma::sp_mat ggvcovmat = vmat * ggicormat * vmat;  
     hmat = ddmat * vcovmat * ddmat.t();
     ghmat = ddmat * gvcovmat * ddmat.t();
     gghmat = ddmat * ggvcovmat * ddmat.t();
     gmat = ddmat * vcovmat * qfullresmat * vcovmat * ddmat.t();
     ggmat = ddmat * gvcovmat * qfullresmat * gvcovmat * ddmat.t();
     gggmat = ddmat * ggvcovmat * qfullresmat * ggvcovmat * ddmat.t();

    }


    if(diffmeth == "analytic"){

     // output
     return Rcpp::List::create(Rcpp::Named("icormat")=icormat,
                               Rcpp::Named("gmat")=gmat,
                               Rcpp::Named("ggmat")=ggmat,
                               Rcpp::Named("gggmat")=gggmat,
                               Rcpp::Named("hmat")=hmat,
                               Rcpp::Named("ghmat")=ghmat,
                               Rcpp::Named("gghmat")=gghmat);

     } else if(diffmeth == "numeric"){

     // output
     return Rcpp::List::create(Rcpp::Named("icormat")=icormat,
                               Rcpp::Named("gmat")=gmat,
                               Rcpp::Named("lgmat")=ggmat,
                               Rcpp::Named("ugmat")=gggmat,
                               Rcpp::Named("hmat")=hmat,
                               Rcpp::Named("lhmat")=ghmat,
                               Rcpp::Named("uhmat")=gghmat);

     } else {

     // output
     return 0;

    }

}
