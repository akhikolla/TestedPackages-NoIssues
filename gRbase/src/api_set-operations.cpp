/* ****************************************************

  Set operations for gRbase and related packages
  
  Author: Søren Højsgaard

  ***************************************************** */

#include <RcppArmadillo.h>
#include "R_like.h"
//[[Rcpp::interfaces(r,cpp)]]
//[[Rcpp::depends(RcppEigen,RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

typedef Rcpp::NumericVector   numVec;
typedef Rcpp::IntegerVector   intVec;
typedef Rcpp::CharacterVector chrVec;
typedef Rcpp::LogicalVector   logVec;


//[[Rcpp::export]]
IntegerVector which_ (SEXP x){
  arma::vec u = as<arma::vec>(x);
  arma::uvec r = find(abs(u) > 1e-6);
  //Rf_PrintValue(wrap(r));
  intVec out = IntegerVector(r.begin(), r.end());
  //Rf_PrintValue(out);
  return out;
}






/* *************************************************

   Implements: 
   get_superset, 
   get_subset
   
   ************************************************* */

// FIXME: REPLACES get_superset_
//[[Rcpp::export]]
IntegerVector get_superset_ (CharacterVector x, List setlist, bool all=false){
  IntegerVector vec(setlist.length());
  int k=0;
  
  for (int i=0; i<setlist.length(); ++i){
    CharacterVector set=setlist[i];
    bool not_contained = (any(is_na(match(x, set))));
    if (!not_contained){
      vec[k++] = i+1;
      if (!all) break;
    }
  }
  
  IntegerVector out = IntegerVector(k);
  if (k > 0) for (int i=0; i<k; ++i) out[i]=vec[i];
  
  return out;
}



// FIXME: REPLACES get_subset_
//[[Rcpp::export]]
IntegerVector get_subset_ (CharacterVector x, List setlist, bool all=false){
  IntegerVector vec(setlist.length());
  int k=0;
  
  for (int i=0; i<setlist.length(); ++i){
    CharacterVector set=setlist[i];
    bool not_contained = (any(is_na(match(set, x))));
    if (!not_contained){
      vec[k++] = i+1;
      if (!all) break;
    }
  }
  
  IntegerVector out = IntegerVector(k);
  if (k > 0) for (int i=0; i<k; ++i) out[i]=vec[i];
  return out;
}






/* *************************************************

   Implements: 
   
   is_subsetof
   
   ************************************************* */

// is_subsetof_ is used in gRain

template <int RTYPE>
bool is_subsetof_impl(Vector<RTYPE> set, Vector<RTYPE> set2) {
  if (set.length() > set2.length()) return false;
  else {
    IntegerVector m = match(set, set2);
    //Rf_PrintValue(m);
    bool out = any(is_na(m));
    return !out;
  }
}

// [[Rcpp::export]]
bool is_subsetof_(SEXP set, SEXP set2) {
    switch(TYPEOF(set)) {
    case INTSXP:  return is_subsetof_impl<INTSXP>(set, set2);
    case REALSXP: return is_subsetof_impl<REALSXP>(set, set2);
    case STRSXP:  return is_subsetof_impl<STRSXP>(set, set2);
    default: stop("Unsupported type.");
    }
}


/* **************************************************************

   Implements:

   allSubsets_cpp(x)  - 'any' vector x
   allSubsets0_cpp(x) - integer vector x
      
   *************************************************************** */

// Works for integer input vector
//[[Rcpp::export]]
List allSubsets0_(const IntegerVector& x)
{
  int nx = x.length(), nout=pow(2., nx), i, k, ny=1;
  double z;
  List out( nout );
  out[0] = -1;
  
  for (i=0; i<nx; ++i){
    z = x[i];
    for (k=0; k<ny; ++k){
      IntegerVector tmp = out[k];
      tmp.push_back( z );
      out[ny + k] = tmp;
    }
    ny = 2 * ny;
  }
  
  for (i=1; i<nout; ++i){
    IntegerVector aa=out[i];
    int M = aa.length();
    IntegerVector vv = no_init( M-1 );
    for (k=1; k<M; ++k){
      vv[k-1] = aa[k];
    }
    out[i-1] = vv;
  }
  out.erase(out.end() - 1, out.end());
  return out;
}

template <int RTYPE>
List do_allSubsets (Vector<RTYPE> vn)
{
  IntegerVector sq = seq_len( vn.size() );
  List lst = allSubsets0_( sq );
  int N=lst.size(), i;
  for (i=0; i<N; ++i){
    lst[i] = vn[ IntegerVector( lst[i] )-1 ];
  }
  return lst;
}

// Works for any type of input vector
// [[Rcpp::export]]
SEXP allSubsets_( SEXP& x){
  int type = TYPEOF(x) ; //Rprintf("type=%i\n", type);
  switch( type ){
  case INTSXP  : return allSubsets0_( x ) ;
  case REALSXP : return do_allSubsets<REALSXP>( x ) ;
  case STRSXP  : return do_allSubsets<STRSXP> ( x ) ;
  default: stop("Unsupported type.");
  }
}


/*** R


allSubsets(1:5)
allSubsets(letters[1:5])

library(microbenchmark)
x <- 1:4
microbenchmark(allSubsets0_R( x ), allSubsets1_R( x ), allSubsets(x))

x <- 1:10
microbenchmark(allSubsets0_R( x ), allSubsets1_R( x ), allSubsets(x))

x <- 1:15
microbenchmark(allSubsets0_R( x ), allSubsets(x),times=10)

*/



/* *************************************************

   Implements:

   all_pairs

   ************************************************* */

CharacterMatrix do_names2pairs_(CharacterVector x, CharacterVector y){
  int lenx = x.length(), leny = y.length();
  int i,j,k;
    
  if (leny == 0){
    if (lenx == 1){
      CharacterMatrix out(lenx * leny, 2);
      return out;
    } else {
      CharacterMatrix out(lenx * (lenx - 1) / 2, 2);
      k = 0;
      for (i=0; i<lenx; ++i){
	for (j=i+1; j<lenx; ++j){
	  out(k,   0) = x[i];
	  out(k++, 1) = x[j];
	}
      }
      return out;
    }
  } else {
    CharacterMatrix out(lenx*leny, 2);
    k = 0;
    for (i=0; i<lenx; ++i){
      for (j=0; j<leny; ++j){
	out(k,   0) = x[i];
	out(k++, 1) = y[j];
      }
    }
    return out;
  }
}

CharacterMatrix sortmat_(CharacterMatrix X){
  CharacterMatrix X2 = clone(X);
  CharacterVector tmp(1);
  for (int i=0; i < X2.rows(); ++i){
    if (X2(i, 0) > X2(i, 1)){
      tmp[0]   = X2(i, 0);
      X2(i, 0) = X2(i, 1);
      X2(i, 1) = tmp[0];
    }
  }
  return X2;
}

//[[Rcpp::export]]
SEXP all_pairs__(CharacterVector x,
		 CharacterVector y=CharacterVector(0),
		 bool sort=false,
		 std::string result="matrix"){
  CharacterMatrix out = do_names2pairs_(x, y);
  
  if (sort) out = sortmat_(out);
  
  if (result == "list"){
    List outlist(out.nrow());
    for (int i=0; i < out.nrow(); ++i)
      outlist[i] = out(i, _);
    return (outlist);
  } else { 
    return(out);
  } 
}






inline int get_length(CharacterVector x){return x.length();}

//[[Rcpp::export]]
SEXP max_set_(const List L, bool index=false){
  int nset = L.length();
  intVec len = sapply(L, get_length);

  intVec perm = order2_(len, true);
  //Rcout << "perm: "; print(perm);
  
  List L2 = L[perm - 1];
  //print(L2);
  
  intVec keepvec2 = rep(1, nset);
  //keepvec2[0] = 0;
  for (int i=0; i<nset-1; i++){
    if (keepvec2[i] == 1){
      for (int j=i+1; j<nset; j++){
	if (keepvec2[j] == 1){
	  //Rf_PrintValue(L2[j]);	  Rf_PrintValue(L2[i]);
	  bool omit = is_subsetof_(L2[j], L2[i]);
	  if (omit) keepvec2[j]=0;
	}
      }
    }
  }

  // Rcout << "keepvec2 : "; print(keepvec2);

  if (index){
    intVec one2n = seq(1, nset); 
    intVec perm_inv = match(one2n, perm);
    // Rcout << "perm_inv : "; print(perm_inv);
    
    intVec idx2 = which_(keepvec2);
    // Rcout << "idx2 : "; print(idx2); The good ones in L2
    
    intVec idx1 = match(idx2 + 1, perm_inv) - 1;
    //  Rcout << "idx1 : "; print(idx1);
    
    intVec out = rep(0, nset);
    out[idx1] = 1;
    //Rcout << "out : "; print(out);
    return out;
  }
  else {  
    L2 = L2[(LogicalVector) keepvec2];
    return L2;
  }
}



//[[Rcpp::export]]
SEXP min_set_(const List L, bool index=false){
  int nset = L.length();
  intVec len = sapply(L, get_length);
  intVec perm = order2_(len, false);
  List L2 = L[perm - 1];  
  intVec keepvec2 = rep(1, nset);
  
  for (int i=0; i<nset-1; i++){
    if (keepvec2[i] == 1){
      for (int j=i+1; j<nset; j++){
	if (keepvec2[j] == 1){
	  bool omit = is_subsetof_(L2[i], L2[j]);
	  if (omit) keepvec2[j]=0;
	}
      }
    }
  }
  if (index){
    intVec one2n = seq(1, nset); 
    intVec perm_inv = match(one2n, perm);
    // Rcout << "perm_inv : "; print(perm_inv);
    
    intVec idx2 = which_(keepvec2);
    // Rcout << "idx2 : "; print(idx2); The good ones in L2
    
    intVec idx1 = match(idx2 + 1, perm_inv) - 1;
    //  Rcout << "idx1 : "; print(idx1);
    
    intVec out = rep(0, nset);
    out[idx1] = 1;
    //Rcout << "out : "; print(out);
    return out;
  }
  else {  
    L2 = L2[(LogicalVector) keepvec2];
    return L2;
  }

}

//[[Rcpp::export]]
SEXP isin_(List L, SEXP set, bool index=false){
  int nset = L.length();

  if (index){
    // Rcout << "here.." << endl;
    intVec out = rep(0, nset);
    for (int i=0; i<nset; i++){
      if (is_subsetof_(set, L[i])) out[i] = 1;
    }
    return out;
  } else {
    // Rcout << "there.." << endl;
    for (int i=0; i<nset; i++){
      if (is_subsetof_(set, L[i])) return wrap(true);
    }
    return wrap(false);
  }
}





















/*** R
library(microbenchmark)
x <- letters[4:1]
y <- letters[4:6]
library(gRbase)
names2pairs(x, y, result="matrix")
names2pairsM(x, y)

names2pairs("x", NULL, result="matrix")
names2pairsM("x", character(0))

names2pairs(x, NULL, result="matrix")
names2pairsM(x, character(0))

m<-matrix(c("a","b","b","a"), ncol=2, byrow=T)


names2pairs2 <- function(x, y=NULL, sort=FALSE, result="list"){
  if (is.null(y))
    y <- character(0)
  out <- names2pairsM(x, y, sort, result)
  out
}

sortmat <- function(out){
     idx <- out[, 1] > out[, 2]
     out[idx, 1:2] <- out[idx, 2:1]
out
}

microbenchmark(sortmat(m), sortmat_(m))

microbenchmark(
  names2pairs(x, y, sort=FALSE, result="matrix"),
  names2pairsM(x, y, sort=FALSE, result="matrix"),
  names2pairs(x, y, sort=TRUE, result="matrix"),
  names2pairsM(x, y, sort=TRUE, result="matrix"),
	names2pairs("x", sort=FALSE, result="matrix"),
	names2pairsM("x", sort=FALSE, result="matrix"),
	names2pairs("x", sort=TRUE, result="matrix"),
	names2pairsM("x", sort=TRUE, result="matrix"),
	names2pairs(x, sort=FALSE, result="matrix"),
	names2pairsM(x, sort=FALSE, result="matrix"),
	names2pairs(x, sort=TRUE, result="matrix"),
	names2pairsM(x, sort=TRUE, result="matrix")
	)

*/







/*** R
x1 <- c("b","a")
x2 <- c("a","k")
set <- letters[1:4]
setlist <- list(x1, x2, c("a","b","k"))
str(setlist)

is_subsetof_(x1, set)
is_subsetof_(x2, set)

is_subsetof_(set, x1)
is_subsetof_(set, x2)

get_superset_(x1, setlist)
get_superset_(c("a","r"), setlist)

get_superset_(x1, setlist, all=T)
get_superset_(c("a","r"), setlist, all=T)

setlist <- list(x1, x2, c("a","b","k"), c("a","r","k"))
str(setlist)

get_subset_(x1, setlist)
get_subset_(c("a", "b", "k"), setlist, all=F)
get_subset_(c("a", "b", "k"), setlist, all=T)
get_subset_(x1, setlist, all=T)
get_subset_(c("a","r"), setlist, all=T)

*/

















































