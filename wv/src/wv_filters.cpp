
#include <RcppArmadillo.h>
#include "wv_filters.h"
// We use reverse_vec
#include "armadillo_manipulations.h"


// Creates an index to call functions by memory address.
struct A{
  
  // @title Construct Filter Selection Map
  static std::map<std::string, arma::field<arma::vec> (*)()> create_map()
  {
    // Makes a map of function addresses indexed by a call string
    std::map<std::string, arma::field<arma::vec> (*)()> filterMap; 
    
    
    filterMap["haar"] = &haar_filter;

    return filterMap;
  }
  
  static const std::map<std::string, arma::field<arma::vec> (*)()> filterMap;
};

const std::map<std::string, arma::field<arma::vec> (*)()> A::filterMap =  A::create_map();

//' @title Quadrature Mirror Filter
//' @description Calculate the series quadrature mirror filter (QMF). Requires a series of an even length.
//' @usage qmf(g, inverse)
//' @param g A \code{vector} that contains the filter constants.
//' @param inverse A \code{bool} that indicates whether the inverse quadrature mirror filter is computed. 
//' By default, the inverse quadrature mirror is computed.
//' @return A \code{vector} that contains either the forward QMF (evalute in order) or the inverse QMF (reverse order). 
//' @author James Balamuta
//' @keywords internal
arma::vec qmf(arma::vec g, bool inverse = true) {
  
  unsigned int L = g.n_elem;
  
  arma::vec rev_g = reverse_vec(g);
    
  for(unsigned int i = 0; i < L; i++){
  
    if( (i+!inverse) % 2 != 0){
      rev_g(i) = rev_g(i)*-1;
    }
    
  }
  
  return rev_g;
}

//' @title Haar filter construction
//' @description Creates the haar filter
//' @usage haar_filter()
//' @return A \code{field<vec>} that contains:
//' \itemize{
//'  \item{"L"}{A \code{integer} specifying the length of the filter}
//'  \item{"h"}{A \code{vector} containing the coefficients for the wavelet filter}
//'  \item{"g"}{A \code{vector} containing the coefficients for the scaling filter}
//' }
//' @details
//' This template can be used to increase the amount of filters available for selection.
//' @author James Balamuta
//' @keywords internal
arma::field<arma::vec> haar_filter() {
  
    arma::vec L(1);
    L(0) = 2.0;
    
    arma::vec g(2);
    g.fill(0.7071067811865475);
    
    arma::vec h = qmf(g);
    
    arma::field<arma::vec> out(3);
    
    out(0)=L;
    out(1)=h;
    out(2)=g;
    
    return out;
}

//' @title Select the Wavelet Filter
//' @description Constructs the wavelet filter to be used.
//' @usage select_filter(filter_name)
//' @param filter_name A \code{String} that must receive: \code{"haar"}.
//' @return info A \code{field<vec>} that contains:
//' \itemize{
//'  \item{"L"}{A \code{integer} specifying the length of the filter}
//'  \item{"h"}{A \code{vector} containing the coefficients for the wavelet filter}
//'  \item{"g"}{A \code{vector} containing the coefficients for the scaling filter}
//' }
//' @details 
//' The package is oriented toward using only the haar filter. If the package extends at a later time, then the supporting infrastructure is there.
//' @author James Balamuta
//' @keywords internal
arma::field<arma::vec> select_filter(std::string filter_name = "haar")
{
  
  arma::field<arma::vec> info(3);
  
  std::map<std::string,arma::field<arma::vec> (*)()>::const_iterator it = A::filterMap.find(filter_name);
  if(it != A::filterMap.end())
  {
    //element found;
    info = (*(it->second))();
  }else{
    Rcpp::stop("Wave Filter is not supported! See ?select_filter for supported types."); 
  }
  
  return info;
}
