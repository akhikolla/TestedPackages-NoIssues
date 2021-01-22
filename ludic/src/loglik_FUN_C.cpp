#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]

//'C++ implementation of the pseudo-likelihood computation 
//'
//'@param Bmat \code{K x nB} matrix of the observations to be matched.
//'
//'@param Amat \code{nA x K} matrix the database into which a match is looked for.
//'
//'@param eps_p a vector of length \code{K} giving the prior discrepancy rate 
//'expected from A to B for the positives, for each variable.
//'
//'@param eps_n a vector of length \code{K} giving the prior discrepancy rate 
//'expected from A to B for the negatives, for each variable.
//'
//'@param piA a vector of length \code{K} giving the prior probabilities of 
//'observing each variable in A.
//'
//'@param piB a vector of length \code{K} giving the prior probabilities of 
//'observing each variable in B.
//'
//'@rdname loglikC_bin
//'
//'@description \code{loglikC_bin} implements an even faster C++ implementation of the pseudo-likelihood computation for binary
//'variables
//'
//'@export
// [[Rcpp::export]]
NumericMatrix loglikC_bin(arma::mat Bmat, 
                          arma::mat Amat,
                          NumericVector eps_p,
                          NumericVector eps_n,
                          NumericVector piA,
                          NumericVector piB){
  
  mat B_t = Bmat.t();
  mat one_minus_B_t = 1-B_t;
  colvec e_p = as<colvec>(eps_p);
  colvec e_n = as<colvec>(eps_n);
  colvec pi_A = as<colvec>(piA);
  colvec pi_B = as<colvec>(piB);
  
  //   int K = A.n_rows;
  //   int n = B_t.n_cols;
  
  mat oneone = B_t.each_col() % log((1-e_n)/pi_B);
  mat zerozero = one_minus_B_t.each_col() % log((1-e_p)/(1-pi_B));
  mat onezero = one_minus_B_t.each_col() % log(e_n/(1-pi_B));
  mat zeroone = B_t.each_col() % log(e_p/pi_B);
  
  mat res = Amat*oneone + (1-Amat)*zerozero + Amat*onezero + (1-Amat)*zeroone;
  
  return(wrap(res));
}





//'Splitting a character string in C++
//' 
//'@param s a character string to be split
//'
//'@param sep a character that delimits the splits
//' 
//'@keywords internal
//' 
//'@examples
//'strsplitC(c(";aaa;bb;cccc;ee;"), sep=";")
//'@export
// [[Rcpp::export]]
std::vector<std::string> strsplitC(std::string s, char sep){
  
  int len = s.size();
  int nb_split = 0;
  
  for(int i=0; i < len; i++ ){
    if(s[i] == sep){
      nb_split += 1;
    }
  }
  
  
  std::vector<std::string> s_splitted = std::vector<std::string>(nb_split+1);
  int which_split = 0;
  IntegerVector split_index = IntegerVector(nb_split+1);
  split_index[0] = -1;
  int nb_empty =0;
  
  for(int i=0; i < len; i++ ){
    if(s[i] == sep){
      which_split += 1;
      split_index[which_split] = i;
      s_splitted[which_split-1] = (s.substr(split_index[which_split-1]+1, i-(split_index[which_split-1]+1)));
      if(s_splitted[which_split-1] == ""){
        nb_empty += 1;
      }
    }
  }
  
  
  if(split_index[which_split] != len-1){
    s_splitted[which_split] = s.substr(split_index[which_split]+1, len);
  }
  if(s_splitted[which_split] == ""){
    nb_empty += 1;
  }
  
  if(nb_empty>0){
    int which_notempty = 0;
    std::vector<std::string> out = std::vector<std::string>(nb_split + 1 - nb_empty);
    for(int i=0; i<nb_split + 1 ; i++){
      if(s_splitted[i] != ""){
        which_notempty += 1;
        out[which_notempty-1] = s_splitted[i];
      }
    }
    return(out);
  }
  else{
    return(s_splitted);
  }
}





//'@rdname loglikC_bin
//'
//'@description \code{loglikC_bin_wDates} implements a C++ implementation of the pseudo-likelihood computation for binary
//'variables with dates
//'
//'@param Bdates \code{nB x K} matrix of the dates for each observations to be matched.
//'
//'@param Adates \code{nA x K} matrix of the dates for database into which a match is looked for.
//'
//'@export
// [[Rcpp::export]]
NumericMatrix loglikC_bin_wDates(arma::mat Bmat, 
                                 arma::mat Amat,
                                 StringMatrix Bdates, 
                                 StringMatrix Adates,
                                 NumericVector eps_p,
                                 NumericVector eps_n,
                                 NumericVector piA,
                                 NumericVector piB){
  mat B_t = Bmat.t();
  mat one_minus_B_t = 1-B_t;
  colvec e_p = as<colvec>(eps_p);
  colvec e_n = as<colvec>(eps_n);
  colvec pi_A = as<colvec>(piA);
  colvec pi_B = as<colvec>(piB);
  
  int nA = Amat.n_rows;
  int nB = Bmat.n_rows;
  //int K = Amat.n_cols;
  
  NumericMatrix res = NumericMatrix(nA,nB);
  mat restemp;
  
  char sepB = ',';
  char sepA = ';';
  
  for(int j=0; j<nB; j++){
    Rcout << j << nB << std::endl;
    
    colvec Btemp = B_t.col(j);
    colvec oneone = Btemp%log(1/pi_B);
    colvec zerozero = (1-Btemp)%log((1-e_p)/(1-pi_B));
    colvec onezero = (1-Btemp)%log(e_n/(1-pi_B));
    colvec zeroone = Btemp%log(e_p/pi_B);
    
    for(int i=0; i<nA; i++){
      //Rcout << "i: " << i << " ; j: " << j << std::endl;
      rowvec Atemp = Amat.row(i);
      
      uvec index11 = find(Btemp%Atemp.t()==1);
      int nb11 = index11.size();
      
      if(nb11>0){
        
        for(int l=0; l<nb11; l++){
          
          bool PossibleMatch = false ;
          
          std::string datesA_temp_string = Rcpp::as<std::string>(Adates(i, index11(l)));
          std::vector<std::string> datesA_stringvect = strsplitC(datesA_temp_string, sepA);
          int nbDates_A = datesA_stringvect.size();
          
          if(nbDates_A>0){
            
            std::string datesB_temp_string = Rcpp::as<std::string>(Bdates(j, index11(l)));
            std::vector<std::string> datesB_stringvect = strsplitC(datesB_temp_string, sepB);
            int nbDates_B = datesB_stringvect.size();
            if(nbDates_B>0){
              //Rcout << "i: " << i << " ; j: " << j << " ; l: " << l << std::endl;
              int da=0;
              while(da<nbDates_A && !PossibleMatch){
                std::vector<std::string> date_range = strsplitC(datesA_stringvect[da], ':');
                Date datea_debut = Date(date_range[0], "%d%b%Y");
                Date datea_fin = Date(date_range[1], "%d%b%Y");
                int db=0; 
                while(db<nbDates_B && !PossibleMatch){
                  Date dateb_test = Date(datesB_stringvect[db], " %b %d %Y");
                  //Rcout << "Strings: " <<  datesB_stringvect[db] << " in " << date_range[0] << ":" << date_range[1] << std::endl;
                  //Rcout << "Dates: " <<  wrap(dateb_test) << " in " << wrap(datea_debut) << ":" << wrap(datea_fin) << std::endl;
                  //Rcout << "test 1: " << (dateb_test >= datea_debut) << std::endl ;
                  //Rcout << "test 2: " << (dateb_test <= datea_fin) << std::endl ;
                  if((dateb_test >= datea_debut) && (dateb_test <= datea_fin)){
                    //if((Date(date_range[0], "%d%b%Y") - Date(datesB_stringvect[db], "%b%d%Y"))<=360){
                    PossibleMatch = true ;
                    //Rcout << "Possible match: " << PossibleMatch << std::endl ;
                  }
                  db += 1;
                }
                da += 1;
              }
            }
            if(!PossibleMatch){
              oneone(index11(l)) = 0;
              //Rcout << "No Match" << std::endl;
            }
          }
        }
      }
      vec lik_temp = Atemp*oneone + (1-Atemp)*zerozero + Atemp*onezero + (1-Atemp)*zeroone;
      res(i,j) = lik_temp(0);
    }
  }
  
  return(res);
}















//'Fast C++ implementation of the log-likelihood ratio computation for 
//'differentiating variables
//'
//'@param d_max a numeric vector of length \code{K} giving the minimum difference 
//'from which it is considered a discrepancy.
//'
//'@param cost a numeric vector of length \code{K} giving the arbitrary cost of discrepancy.
//'
//'@rdname loglikC_bin
//'
//'@export
// [[Rcpp::export]]
NumericMatrix loglikratioC_diff_arbitrary(arma::mat Bmat, 
                                          arma::mat Amat,
                                          NumericVector d_max,
                                          NumericVector cost){
  
  colvec cc = as<colvec>(cost);
  colvec dmax = as<colvec>(d_max);
  
  int nA = Amat.n_cols;
  int nB = Bmat.n_cols;
  
  mat res = mat(nA,nB);
  
  for(int j=0; j<nB; j++){
    mat A_minus_Btemp = abs(Amat.each_col() - Bmat.col(j));
    umat A_minus_B_inf_u = A_minus_Btemp < repmat(dmax, 1, nA);
    mat A_minus_B_inf = conv_to<mat>::from(A_minus_B_inf_u);
    mat A_minus_B_sup = (1.0-A_minus_B_inf);
    mat inf_mat = (A_minus_B_sup.each_col() % cc);
    rowvec lik_temp = sum(inf_mat, 0);
    res.col(j) = lik_temp.t();
  }
  
  return(wrap(res));
}
