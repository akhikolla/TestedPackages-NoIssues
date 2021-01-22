#include <RcppArmadillo.h>

// [[Rcpp::export(name = ".ACtestCpp")]]
SEXP ACtestCpp(arma::mat z_e, arma::mat z, arma::mat e, int h, bool univariate,
               bool LM, bool HC0, bool HC1, bool HC2, bool HC3){
  // This function computes the test statistics for ACtest.
  
  // z = the (N - p) x (K * p + numberOf(const, trend, exogen)) matrix of regressors from the estimated VAR(p)
  // e = the (N - p) x K matrix of residuals from the estimated VAR(p)
  // z_e = 'z' concatenated with the h lagged residuals
  // h = lag length of the test
  
  int K = e.n_cols, N = e.n_rows;
  
  // defines matrices:
  arma::mat W0(K*z_e.n_cols, K*z_e.n_cols, arma::fill::zeros),
  W2(K*z_e.n_cols, K*z_e.n_cols, arma::fill::zeros),
  W3(K*z_e.n_cols, K*z_e.n_cols, arma::fill::zeros);
  arma::mat V(K*z_e.n_cols, K*z_e.n_cols, arma::fill::zeros);
  arma::mat tempMat(z_e.n_cols, z_e.n_cols, arma::fill::zeros),
  ee(K, K, arma::fill::zeros);
  
  // holders for the return values:
  Rcpp::NumericMatrix multiReturnValue(1, 5), uniReturnValue(K, 5);
  std::fill(multiReturnValue.begin(), multiReturnValue.end(), NA_REAL);
  std::fill(uniReturnValue.begin(), uniReturnValue.end(), NA_REAL);
  Rcpp::colnames(multiReturnValue) = Rcpp::CharacterVector::create("LM", "HC0", "HC1", "HC2", "HC3");
  Rcpp::colnames(uniReturnValue) = Rcpp::CharacterVector::create("LM", "HC0", "HC1", "HC2", "HC3");
  
  // computes the non-HCCME version of eq. (6):
  if(LM){
    arma::mat e_ = z_e(arma::span::all, arma::span(z.n_cols, z_e.n_cols - 1));
    arma::mat Sc = e_.t() * e_ - e_.t() * z * arma::inv(z.t() * z) * z.t() * e_;
    Sc = kron(Sc, e.t() * e); // as in eq. (9)
    arma::vec cv = vectorise(e.t() * e_); // as in eq. (8)
    tempMat = N * cv.t() * arma::inv(Sc) * cv;
    multiReturnValue(0, 0) = tempMat(0,0);
  }
  
  // continues with the HCCME versions of eq. (6):
  
  // scalings for HC1 - HC3:
  double ht1 = (double)N /(N - K - 1); 
  arma::vec ht2 = 1 - arma::diagvec((z * arma::inv(z.t() * z).t()) * z.t());
  
  if(HC0 || HC1 || HC2 || HC3){
    
    // computes the V and W matricies from eq. (11) in Ahlgren & Catani (2016):
    for(size_t row = 0; row < z_e.n_rows; ++row){
      
      tempMat = z_e.row(row).t() * z_e.row(row);
      ee = e.row(row).t() * e.row(row);
      
      if(HC0 || HC1) W0 += kron(tempMat, ee);
      if(HC2) W2 += kron(tempMat, ee / ht2(row));
      if(HC3) W3 += kron(tempMat, ee / pow(ht2(row), 2));
      V += kron(tempMat, arma::eye( K, K ));
    }
    
    // Regresses 'e' on ' 'z_e':
    arma::mat olsCoef = arma::inv(z_e.t() * z_e) * z_e.t() * e;
    // takes the relevant coefficients and vectorises them to the form of psi in eq. (5):
    arma::vec olsCoefPsi = arma::vectorise(olsCoef(arma::span(z_e.n_cols - K * h, z_e.n_cols - 1), arma::span::all).t());
    
    // inverts 'V':
    arma::mat Vinv = arma::inv(V);
    
    // HC0:
    if(HC0 || HC1){
      // computes the matrix in eq. (11):
      arma::mat VWV0 = Vinv * W0 * Vinv;
      
      // extracts the relevant block and computes HC0 as in eq. (6):
      tempMat = olsCoefPsi.t() *
        arma::inv(VWV0(arma::span(V.n_cols - K * K * h, V.n_cols - 1), arma::span(V.n_cols - K * K * h, V.n_cols - 1))) *
        olsCoefPsi;
      multiReturnValue(0, 1) = tempMat(0,0);
    }
    
    // HC1:
    if(HC1) multiReturnValue(0, 2) = multiReturnValue(0, 1) / ht1;
    
    // HC2:
    if(HC2){
      // computes the matrix in eq. (11):
      arma::mat VWV2 = Vinv * W2 * Vinv;
      
      // extracts the relevant block and computes HC0 as in eq. (6):
      tempMat = olsCoefPsi.t() *
        arma::inv(VWV2(arma::span(V.n_cols - K * K * h, V.n_cols - 1), arma::span(V.n_cols - K * K * h, V.n_cols - 1))) *
        olsCoefPsi;
      multiReturnValue(0, 3) = tempMat(0,0);
    }
    
    // HC3:
    if(HC3){
      // computes the matrix in eq. (11):
      arma::mat VWV3 = Vinv * W3 * Vinv;
      
      // extracts the relevant block and computes HC0 as in eq. (6):
      tempMat = olsCoefPsi.t() *
        arma::inv(VWV3(arma::span(V.n_cols - K * K * h, V.n_cols - 1), arma::span(V.n_cols - K * K * h, V.n_cols - 1))) *
        olsCoefPsi;
      multiReturnValue(0, 4) = tempMat(0,0);
    }
  }
  
  //defines univariate variables:
  arma::uvec uniIndices;
  arma::mat W0uni, W2uni, W3uni;
  arma::mat olsCoefUni, z_eUni;
  arma::vec olsCoefPsiUni;
  arma::mat VUni, VWVuni;
  arma::mat e_uni, ScUni; 
  arma::vec cvUni;
  
  // computes the univariate versions of eq. (6):
  if(univariate){
    for(int j = 0; j < K; ++j){
      
      // computes the non-HCCME version LM of eq. (6):
      if(LM){
        e_uni = z_e(arma::regspace<arma::uvec>(0, z_e.n_rows - 1), arma::regspace<arma::uvec>(z.n_cols + j, K, z_e.n_cols - 1));
        ScUni = e_uni.t() * e_uni - e_uni.t() * z * arma::inv(z.t() * z) * z.t() * e_uni;
        ScUni = kron(ScUni, e.col(j).t() * e.col(j));
        cvUni = vectorise(e.col(j).t() * e_uni);
        tempMat = N * cvUni.t() * arma::inv(ScUni) * cvUni;
        uniReturnValue(j, 0) = tempMat(0,0);
      }
      
      // continues with the HCCME versions of eq. (6):
      if(HC0 || HC1 || HC2 || HC3){
        // computes the indexes of the W and V matricies that corresponds to the univariate versions of them:
        uniIndices = arma::regspace<arma::uvec>(0 + j, K, z.n_cols * K - 1);
        uniIndices = arma::join_cols(uniIndices, arma::regspace<arma::uvec>(z.n_cols * K + j * (K + 1), K * K, z_e.n_cols * K - 1));
        
        // creates the univariate version of 'z_e':
        z_eUni = arma::join_rows(z, z_e.cols(arma::regspace<arma::uvec>(z.n_cols + j, K, z_e.n_cols - 1)));
        
        // Regresses 'e.col(j)' on ' 'z_eUni' and takes the univariate psi as in eq. (5):
        olsCoefUni = arma::inv(z_eUni.t() * z_eUni) * z_eUni.t() * e.col(j);
        olsCoefPsiUni = olsCoefUni(arma::span(z_eUni.n_cols - h, z_eUni.n_cols - 1), 0);
        
        VUni = arma::inv(V(uniIndices, uniIndices));
        
        // HC0:
        if(HC0 || HC1){
          // computes the matrix in eq. (11):
          VWVuni = VUni * W0(uniIndices, uniIndices) * VUni;
          
          // extracts the relevant block and computes HC0 as in eq. (6):
          arma::mat tempV = olsCoefPsiUni.t() * 
            arma::inv(VWVuni(arma::span(VUni.n_cols - h, VUni.n_cols - 1), arma::span(VUni.n_cols - h, VUni.n_cols - 1))) * 
            olsCoefPsiUni;
          uniReturnValue(j, 1) = tempV(0,0);
        }
        
        // HC1:
        if(HC1) uniReturnValue(j, 2) = uniReturnValue(j, 1) / ht1;
        
        // HC2:
        if(HC2){
          // computes the matrix in eq. (11):
          VWVuni = VUni * W2(uniIndices, uniIndices) * VUni;
          
          // extracts the relevant block and computes HC0 as in eq. (6):
          arma::mat tempV = olsCoefPsiUni.t() *
            arma::inv(VWVuni(arma::span(VUni.n_cols - h, VUni.n_cols - 1), arma::span(VUni.n_cols - h, VUni.n_cols - 1))) *
            olsCoefPsiUni;
          uniReturnValue(j, 3) = tempV(0,0);
        }
        
        // HC3:
        if(HC3){
          // computes the matrix in eq. (11):
          VWVuni = VUni * W3(uniIndices, uniIndices) * VUni;
          
          // extracts the relevant block and computes HC0 as in eq. (6):
          arma::mat tempV = olsCoefPsiUni.t() *
            arma::inv(VWVuni(arma::span(VUni.n_cols - h, VUni.n_cols - 1), arma::span(VUni.n_cols - h, VUni.n_cols - 1))) *
            olsCoefPsiUni;
          uniReturnValue(j, 4) = tempV(0,0);
        }
        
      }
    }
  }
  
  Rcpp::List returnList;  
  returnList["multiReturnValue"] = multiReturnValue;
  returnList["uniReturnValue"] = uniReturnValue;
  return returnList;
}
