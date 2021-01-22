#include <RcppArmadillo.h>

// [[Rcpp::export(name = ".makeVar")]]
arma::mat makeVar(arma::mat Ystart, 
                  arma::mat e, 
                  arma::mat param,
                  int p, 
                  bool constFlag, 
                  bool trendFlag, 
                  bool exogenFlag, 
                  arma::mat exogen){
  // this function simulates a VAR(p)
  // Ystart is the first p start observations
  // e is the (N-p) x K matrix of errors, where N is the length of the output simulated series
  // param is the (p * K + numberOf(constFlag, trendFlag, exogen.columns)) x K matrix of parameters 
  //  (ordered as in the 'VARfit' function)
  // p is the lag order of VAR(p)
  // exogen is the matrix of exogenous variables of length (N-p)
  
  int K = e.n_cols, N = e.n_rows + p, exoCols = 0;
  double value = 0;
  
  if(exogenFlag) exoCols = exogen.n_cols;
  if((int)param.n_rows != (K * p + constFlag + trendFlag + exoCols) ||
     param.n_cols != e.n_cols) Rcpp::stop("Wrong size of the 'param' matrix");
  
  arma::mat Y(N, K);
  
  // fills the first rows with 'Ystart'
  Y(arma::span(0, p - 1), arma::span(0, K - 1)) = Ystart;
  
  // adds the errors
  Y(arma::span(p, N - 1), arma::span(0, K - 1)) = e;
  
  // adds constants
  if(constFlag == true){
    for(int i = p; i < N; ++i) Y.row(i) = Y.row(i) + param.row(0);
  } 
  
  // adds trend
  if(trendFlag == true){
    
    arma::vec oneToEnd(N - p);
    for(int i = 0; i < N - p; ++i) oneToEnd(i) = i + 1;
    
    Y(arma::span(p, N - 1), arma::span(0, K - 1)) = 
      Y(arma::span(p, N - 1), arma::span(0, K - 1))
      + oneToEnd * param.row(constFlag);
  }
  
  // adds exogen
  if(exogenFlag == true){
    int exogenRowIndex = constFlag + trendFlag;
    
    for(int i = 0; i < exoCols; ++i){
      
      Y(arma::span(p, N - 1), arma::span(0, K - 1)) = 
        Y(arma::span(p, N - 1), arma::span(0, K - 1))
      + exogen.col(i) * param.row(exogenRowIndex);
    }
    
  }
  
  // adds the endogenous part:
  int endoParaIndex = constFlag + trendFlag;
  if(exogenFlag == true) endoParaIndex += exoCols;
  for(int row = p; row < N; ++row){
    
    for(int j = 0; j < K; ++j){
      
      value = 0;
      for(int k = 0; k < K; ++k){
        for(int l = 0; l < p; ++l){
          
          value += param(endoParaIndex + k + K * l, j) * Y(row - l - 1, k);
          
        }
      }
      
      Y(row, j) += value;
    }
  }
  
  return(Y);
}
