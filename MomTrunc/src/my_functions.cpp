// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(mvtnorm)]]
#include <mvtnormAPI.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <R_ext/Rdynload.h>  /* required by R */

#include <Rmath.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat TT_GS_sp(arma::uword n, arma::mat R, double nu, arma::vec x, arma::vec lower, arma::vec upper) {
  arma::mat Rinv=R.i();
  //  vec x=(lower+upper)/2;
  uword p=R.n_cols;
  double nup=nu+p;
  mat X(n,p,fill::zeros);
  umat minusj(p-1,p,fill::zeros);
  for(uword j=0;j<p;j++){
    uword k=0;
    for(uword l=0;l<p;l++){
      if(l!=j){
        minusj(k,j)=l;
        k++;
      }
    }
  }
  
  for(uword i=0;i<n;i++){
    double delta=as_scalar(x.t()*Rinv*x);
    double y=arma::randu<double>()*exp(-nup*log(1+delta/nu)/2);
    double kap = nu * (pow(y,-2/nup) - 1);
    for(uword j=0;j<p;j++){
      uvec pj=minusj.col(j);
      vec xj=x(pj);
      rowvec a1=xj.t()*Rinv.rows(pj);
      double ss=as_scalar(a1.cols(pj)*xj);
      double mj=-a1(j)/Rinv(j,j);
      double tj=sqrt(mj*mj+(kap-ss)/Rinv(j,j));
      double lv=(lower(j)<mj-tj)?(mj-tj):(lower(j));
      double rv=(upper(j)<mj+tj)?(upper(j)):(mj+tj);
      double xij=lv+(rv-lv)*arma::randu<double>();
      X(i,j)=xij;
      x(j)=xij;
    }
  }
  
  return X;
}


//[[Rcpp::export]]
arma::vec triangl(const arma::mat& X){
  int n = X.n_cols;
  arma::vec res(n * (n-1) / 2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      res(j + i * (i-1) / 2) = X(i, j);
    }
  }
  return res;
}

// [[Rcpp::export]]
double Fpmvt_cpp(arma::vec& upper,
                 arma::mat& R, int nu = 0,
                 double abseps = 1e-3){
  
  arma::vec lowertrivec = triangl(R);
  
  int n = upper.n_elem;
  int maxpts = 25000;     // default in mvtnorm: 25000
  double releps = 0;      // default in mvtnorm: 0
  int rnd = 1;            // Get/PutRNGstate
  
  double* upper_ = upper.memptr();
  double* correlationMatrix = lowertrivec.memptr();
  double* lower = new double[n];
  int* infin = new int[n];
  double* delta = new double[n];
  
  for (int i = 0; i < n; ++i) {
    infin[i] = 0; // (-inf, bound]
    lower[i] = 0.0;
    delta[i] = 0.0;
  }
  
  // return values
  double error;
  double value;
  int inform;
  
  mvtnorm_C_mvtdst(&n, &nu, lower, upper_,
                   infin, correlationMatrix, delta,
                   &maxpts, &abseps, &releps,
                   &error, &value, &inform, &rnd);
                   delete[] (lower);
                   delete[] (infin);
                   delete[] (delta);
                   
                   return value;
}


// [[Rcpp::export]]
double pmvt_cpp(arma::vec& lower,arma::vec& upper,
                arma::mat& R, int nu = 0,
                double abseps = 1e-3){
  
  arma::vec lowertrivec = triangl(R);
  
  int n = lower.n_elem;
  int maxpts = 25000;     // default in mvtnorm: 25000
  double releps = 0;      // default in mvtnorm: 0
  int rnd = 1;            // Get/PutRNGstate
  
  double* lower_ = lower.memptr();
  double* upper_ = upper.memptr();
  double* correlationMatrix = lowertrivec.memptr();
  int* infin = new int[n];
  double* delta = new double[n];
  
  int i = 0;
  for (i = 0; i < n; ++i){
    delta[i] = 0.0;
    if(std::isinf(lower_[i])){
      if(std::isinf(upper_[i])){
        infin[i] = -1;
      }else{
        infin[i] = 0;
      }
    }else{
      if(std::isinf(upper_[i])){
        infin[i] = 1;
      }else{
        infin[i] = 2;
      }
    }
  }
  
  // return values
  double error;
  double value;
  int inform;
  
  mvtnorm_C_mvtdst(&n, &nu, lower_, upper_,
                   infin, correlationMatrix, delta,
                   &maxpts, &abseps, &releps,
                   &error, &value, &inform, &rnd);
                   delete[] (infin);
                   delete[] (delta);
                   
                   return value;
}

//*************************************************************************
//**************************** aux functions ******************************
//*************************************************************************

//FACTORIAL FUNCTION IN C++
double fact(int x){
  if(x==0)
    return 1;
  double r = x;
  for(int i = r-1;i>0;i--)
  {
    r = r*i;
  }
  return r;
}

//*************************************************************************
//BINOMIAL COEFFICIENT IN C++
double binom(int n,int k){
  double fn = fact(n);
  double fk = fact(k);
  double fnk = fact(n-k);
  double result = fn/(fk*fnk);
  return result;
}

//*************************************************************************
// REPETITIONS GIVEN EACH
arma::colvec rep_each(const arma::colvec& u, int each){
  int n = u.n_elem;
  arma::colvec result(n*each);
  result = trans(vectorise(repmat(u, 1, each), 1)); // this seems inefficient
  return(result);
}

//*************************************************************************
// REPETITIONS GIVEN TIMES
arma::colvec rep_times(const arma::colvec& u, int times){
  int n = u.n_elem;
  arma::colvec result(n*times);
  result = repmat(u, times, 1); // this seems inefficient
  return(result);
}


//*************************************************************************
//VECTORIZED EXPONENCIALIZATION

struct pow_wrapper {
  public: double operator()(double a, double b) {
    return ::pow(a, b);
  }
};

NumericVector vecpow(const NumericVector base, const NumericVector exp){
  Rcpp::NumericVector out(base.size());
  std::transform(base.begin(), base.end(),
                 exp.begin(), out.begin(), pow_wrapper());
  return out;
}

//*************************************************************************
//SIGNS VECTOR GENERATOR
arma::rowvec DecToSigns(int number, int n){
  arma::rowvec out(n,arma::fill::ones);
  out = -out;
  double c = 0;
  do
  {
    if((number&1)!= 0)
      out(c)=1;
    
    number >>= 1;
    c += 1;
  }while(number);
  return arma::reverse(out);
}

//*************************************************************************
// GRID GENERATOR
// [[Rcpp::export]]
arma::mat mygrid(const arma::vec kappa){
  int n = kappa.n_elem;
  arma::vec dims = kappa + 1;
  arma::mat C(prod(dims),n);
  for(int j=1; j<=n; ++j){
    arma::vec seqj = arma::regspace<arma::vec>(0,1,kappa(j-1));
    arma::uvec index1 = arma::regspace<arma::uvec>(1,1,j-1);
    int each = prod(dims(index1-1));
    arma::uvec index2 = arma::regspace<arma::uvec>(j+1,1,n);
    int times = prod(dims(index2-1));
    arma::vec out1 = rep_each(seqj,each);
    arma::vec out2 = rep_times(out1,times);
    C.col(j-1) = out2;
  }
  return C;
}

//*************************************************************************
//FIND LOCATION: AUX FOR INTEGER COMPOSITION OF K INTEGERS SUMMING N
Rcpp::List genColex(int k1, int r, int count, arma::mat X, arma::rowvec c,int n){
  count += 1;
  X.row(count-1) = c;
  if(k1>0){
    X(count-1,n-1) += k1;
    for(int i = n-1; i>=r; --i){
      c(i-1) += 1;
      Rcpp::List W = genColex(k1-1,i,count,X,c,n);
      count = W["count"];
      X = Rcpp::as<arma::mat>(W["X"]);
      c = Rcpp::as<arma::rowvec>(W["c"]);
      c(i-1) -= 1;
    }
  }
  Rcpp::List out = Rcpp::List::create(_["count"] = count, _["X"] = X, _["c"] = c);
  return out;
}


// [[Rcpp::export]]
arma::mat colex(int k, int n){
  int l;
  int count;
  l = binom(n+k-1,n-1);
  arma::mat X(l,n,arma::fill::zeros);
  arma::rowvec c(n,arma::fill::zeros);
  count = 0;
  Rcpp::List W = genColex(k,1,count,X,c,n);
  X = Rcpp::as<arma::mat>(W["X"]);
  return X;
}

//*************************************************************************
//LOCATION FINDER FOR KAPPA
arma::uvec colexind2(arma::mat kk, arma::mat C){
  int mm = kk.n_rows;
  int nn = kk.n_cols;
  arma::colvec yy(mm,1,arma::fill::ones);
  arma::uvec jjj(1);
  arma::uvec ind = arma::linspace<arma::uvec>(1,nn,nn);
  arma::uvec retro = arma::reverse(ind);
  for(int jj=1; jj<=mm; ++jj){
    jjj = jj;
    arma::rowvec ktemp = kk(jjj-1,retro-1);
    arma::rowvec s = arma::cumsum(ktemp);
    for(int r=1; r<=nn-1; ++r){
      if(kk(jj-1,r-1)>0){
        int n1 = nn-r+1;
        yy(jj-1) += C(s(n1-1)+n1-1,n1-1)-C(s(n1-2)+n1-1,n1-1);
      }
    }
  }
  arma::ucolvec yy2 = arma::conv_to<arma::ucolvec>::from(yy);
  return yy2;
}

//*************************************************************************
//Q-FUNCTION FOR THE MEAN-VAR FOR TRUNCATED NORMAL CASE
// [[Rcpp::export]]
Rcpp::List Rcppqfun(arma::vec a,arma::vec b,arma::mat S){
  Rcpp::List out;
  int n = a.n_elem;
  arma::colvec s = sqrt(arma::diagvec(S));
  if(n==1){
    double s1 = s[0];
    double aa = a[0]/s1;
    double bb = b[0]/s1;
    double qa = R::dnorm(aa,0,1,0)/s1;
    double qb = R::dnorm(bb,0,1,0)/s1;
    out["qa"] = qa;
    out["qb"] = qb;
    return out;
  }else{
    arma::vec qa(n,arma::fill::zeros);
    arma::vec qb(n,arma::fill::zeros);
    arma::uvec seqq = arma::linspace<arma::uvec>(1,n,n);
    double da, db;
    arma::uvec ii(1);
    for(int i=1; i<=n; ++i){
      ii = i;
      arma::uvec ind2 = seqq.elem(find(seqq != i));
      arma::mat SSi = S(ind2-1,ind2-1)-S(ind2-1,ii-1)*S(ii-1,ind2-1)/arma::as_scalar(S(ii-1,ii-1));
      arma::colvec si = sqrt(arma::diagvec(SSi));
      arma::mat Ri = SSi%(1/(si * si.t()));
      if(std::isfinite(a(i-1))){
        da = R::dnorm(a(i-1),0,s(i-1),0);
        arma::vec mai = S(ind2-1,ii-1)/S(i-1,i-1)*a(i-1);
        arma::vec haia = (a(ind2-1)-mai)/si;
        arma::vec hbia = (b(ind2-1)-mai)/si;
        qa(i-1) = da*pmvt_cpp(haia,hbia,Ri);
        //qa(i-1) = da;
      } 
      if(std::isfinite(b(i-1))){
        db = R::dnorm(b(i-1),0,s(i-1),0);
        arma::vec mbi = S(ind2-1,ii-1)/S(i-1,i-1)*b(i-1);
        arma::vec haib = (a(ind2-1)-mbi)/si;
        arma::vec hbib = (b(ind2-1)-mbi)/si;
        qb(i-1) = db*pmvt_cpp(haib,hbib,Ri);
        //qb(i-1) = db;
      }
    }
    out["qa"] = qa;
    out["qb"] = qb;
    return out;
  }
}

//*************************************************************************
//Q-FUNCTION FOR THE MEAN-VAR FOR TRUNCATED NORMAL CASE
// [[Rcpp::export]]
Rcpp::List Rcppqfun_ab(arma::vec a,arma::vec b,arma::mat S){
  Rcpp::List out;
  int n = a.n_elem;
  arma::colvec s = sqrt(arma::diagvec(S));
  if(n==1){
    double s1 = s[0];
    double aa = a[0]/s1;
    double bb = b[0]/s1;
    double qa = R::dnorm(aa,0,1,0)/s1;
    double qb = R::dnorm(bb,0,1,0)/s1;
    out["qa"] = qa;
    out["qb"] = qb;
    return out;
  }else{
    arma::vec qa(n,arma::fill::zeros);
    arma::vec qb(n,arma::fill::zeros);
    arma::uvec seqq = arma::linspace<arma::uvec>(1,n,n);
    double da, db;
    arma::uvec ii(1);
    for(int i=1; i<=n; ++i){
      ii = i;
      arma::uvec ind2 = seqq.elem(find(seqq != i));
      arma::mat SSi = S(ind2-1,ind2-1)-S(ind2-1,ii-1)*S(ii-1,ind2-1)/arma::as_scalar(S(ii-1,ii-1));
      arma::colvec si = sqrt(arma::diagvec(SSi));
      arma::mat Ri = SSi%(1/(si * si.t()));

        da = R::dnorm(a(i-1),0,s(i-1),0);
        arma::vec mai = S(ind2-1,ii-1)/S(i-1,i-1)*a(i-1);
        arma::vec haia = (a(ind2-1)-mai)/si;
        arma::vec hbia = (b(ind2-1)-mai)/si;
        qa(i-1) = da*pmvt_cpp(haia,hbia,Ri);
        
        db = R::dnorm(b(i-1),0,s(i-1),0);
        arma::vec mbi = S(ind2-1,ii-1)/S(i-1,i-1)*b(i-1);
        arma::vec haib = (a(ind2-1)-mbi)/si;
        arma::vec hbib = (b(ind2-1)-mbi)/si;
        qb(i-1) = db*pmvt_cpp(haib,hbib,Ri);
        //qb(i-1) = db;
    }
    out["qa"] = qa;
    out["qb"] = qb;
    return out;
  }
}

//*************************************************************************
//Q-FUNCTION FOR THE MEAN-VAR FOR TRUNCATED NORMAL CASE
// [[Rcpp::export]]
arma::vec Rcppqfun_b(arma::vec b,arma::mat S){
  int n = b.n_elem;
  arma::vec qb(n,arma::fill::zeros);
  arma::colvec s = sqrt(arma::diagvec(S));
  if(n==1){
    double s1 = s[0];
    double bb = b[0]/s1;
    double qbn = R::dnorm(bb,0,1,0)/s1;
    qb[0] = qbn;
    return qb;
  }else{
    arma::uvec seqq = arma::linspace<arma::uvec>(1,n,n);
    double db;
    arma::uvec ii(1);
    for(int i=1; i<=n; ++i){
      ii = i;
      arma::uvec ind2 = seqq.elem(find(seqq != i));
      arma::mat SSi = S(ind2-1,ind2-1)-S(ind2-1,ii-1)*S(ii-1,ind2-1)/arma::as_scalar(S(ii-1,ii-1));
      arma::colvec si = sqrt(arma::diagvec(SSi));
      arma::mat Ri = SSi%(1/(si * si.t()));
      
      if(std::isfinite(b(i-1))){
        db = R::dnorm(b(i-1),0,s(i-1),0);
        arma::vec mbi = S(ind2-1,ii-1)/S(i-1,i-1)*b(i-1);
        arma::vec hbib = (b(ind2-1)-mbi)/si;
        qb(i-1) = db*Fpmvt_cpp(hbib,Ri);
        //qb(i-1) = db;
      } 
    }
    return qb;
  }
}

//*************************************************************************
//F-FUNCTION NORMAL CASE INTERNAL (NO GRID)
arma::vec recintab0(arma::vec kappa,arma::vec a,arma::vec b, arma::vec mu, arma::mat S){
  int n = kappa.n_elem;
  arma::vec M;
  if(n==1){
    M.set_size(kappa[0]+1);
    M.zeros();
    double a0 = arma::as_scalar(a);
    double b0 = arma::as_scalar(b);
    double mu0 = arma::as_scalar(mu);
    double S0 = arma::as_scalar(S);
    double s1 = sqrt(S0);
    double aa = (a0-mu0)/s1;
    double bb = (b0-mu0)/s1;
    M[0] = R::pnorm(bb,0,1,1,0) - R::pnorm(aa,0,1,1,0);
    if(kappa[0]>0){
      double pdfa = s1*R::dnorm(aa,0,1,0);
      double pdfb = s1*R::dnorm(bb,0,1,0);
      M(1) = mu0*M(0) + pdfa - pdfb;
      if(std::isinf(a0)){a0 = 0;} 
      if(std::isinf(b0)){b0 = 0;}
      for(int i=2; i<= kappa[0]; ++i){
        pdfa = pdfa*a0;
        pdfb = pdfb*b0;
        M(i) = mu0*M(i-1) + (i-1)*S0*M(i-2) + pdfa - pdfb;
      }
    }
    return M;
  }else{
    arma::uvec seqq = arma::linspace<arma::uvec>(1,n,n);
    int pk = prod(kappa+1);
    M.set_size(pk);
    M.zeros();
    // #
    // #  We create two long vectors G and H to store the two different sets
    // #  of integrals with dimension n-1.
    // #
    arma::vec nn = round(pk/(kappa+1));
    arma::vec begind = arma::cumsum(nn);
    begind = join_cols(arma::zeros<arma::vec>(1),begind);
    int pk1 = begind(n);
    // # number of (n-1)-dimensional integrals
    // #  Each row of cp corresponds to a vector that allows us to map the subscripts
    // #  to the index in the long vectors G and H
    arma::mat cp(n,n,arma::fill::zeros);
    for(int i=1; i<=n; ++i){
      arma::vec kk = kappa;
      kk(i-1) = 0;
      arma::vec cump = cumprod(kk(arma::span(0,n-2))+1);
      cp.row(i-1) = join_horiz(arma::ones<arma::rowvec>(1),cump.t());
    }
    arma::vec G(pk1,arma::fill::zeros);
    arma::vec H(pk1,arma::fill::zeros);
    arma::colvec s = sqrt(arma::diagvec(S));
    Rcpp::NumericVector ha = Rcpp::wrap((a-mu)/s);
    Rcpp::NumericVector hb = Rcpp::wrap((b-mu)/s);
    arma::mat R = S%(1/(s * s.t()));
    arma::colvec pdfa = Rcpp::dnorm(ha);
    pdfa = pdfa/s;
    arma::colvec pdfb = Rcpp::dnorm(hb);
    pdfb = pdfb/s;
    arma::uvec ii(1);
    for(int i=1; i<= n; ++i){
      ii = i;
      arma::uvec ind2 = seqq.elem(find(seqq != i));
      arma::vec kappai = kappa(ind2-1);
      arma::vec ai = a(ind2-1);
      arma::vec bi = b(ind2-1);
      arma::vec mui = mu(ind2-1);
      arma::mat Si = S(ind2-1,ii-1);
      arma::mat SSi = S(ind2-1,ind2-1)-Si*Si.t()/arma::as_scalar(S(ii-1,ii-1));
      arma::uvec ind = arma::regspace<arma::uvec>(begind(i-1)+1,1,begind(i));
      if(std::isfinite(a(i-1))){
        arma::vec mai = mui + Si/S(i-1,i-1)*(a(i-1)-mu(i-1));
        G(ind-1) = pdfa(i-1)*recintab0(kappai,ai,bi,mai,SSi);
      } 
      if(std::isfinite(b(i-1))){
        arma::vec mbi = mui + Si/S(i-1,i-1)*(b(i-1)-mu(i-1));
        H(ind-1) = pdfb(i-1)*recintab0(kappai,ai,bi,mbi,SSi);
      }
    }
    arma::vec a1 = a - mu;
    arma::vec b1 = b - mu;
    M(0) = pmvt_cpp(a1,b1,R);
    //M(0) = 1;
    a.elem(find_nonfinite(a)).zeros();
    b.elem(find_nonfinite(b)).zeros();
    arma::rowvec cp1 = cp.row(n-1);
    arma::mat grid = mygrid(kappa) + 1;
    for(int i=2; i<=pk; ++i){
      arma::rowvec kk = grid.row(i-1);
      //int xx = sum((kk-1)%cp1)+1;
      int i1 = min(seqq.elem(find(kk>1)));
      arma::rowvec kk1 = kk;
      kk1(i1-1)--;
      int ind3 = i - cp1(i1-1);
      M(i-1) = mu(i1-1)*M(ind3-1);
      for(int j=1; j<=n; ++j){
        int kk2 = kk1(j-1)-1;
        if(kk2 > 0){
          M(i-1) += S(i1-1,j-1)*kk2*M(ind3-cp1(j-1)-1);
        }
        int ind4 = begind(j-1) + sum(cp.row(j-1)%(kk1-1)) - cp(j-1,j-1)*kk2 + 1;
        M(i-1) += S(i1-1,j-1)*(pow(a(j-1),kk2)*G(ind4-1)-pow(b(j-1),kk2)*H(ind4-1));
      }
    }
    return M;
  }
}

//*************************************************************************
//F-FUNCTION NORMAL CASE
Rcpp::List recintab1(arma::vec kappa,arma::vec a,arma::vec b, arma::vec mu, arma::mat S){
  Rcpp::List out;
  int n = kappa.n_elem;
  arma::vec M;
  arma::mat index;
  if(n==1){
    M.set_size(kappa[0]+1);
    M.zeros();
    index = arma::linspace<arma::mat>(0,kappa[0],kappa[0]+1);
    double a0 = arma::as_scalar(a);
    double b0 = arma::as_scalar(b);
    double mu0 = arma::as_scalar(mu);
    double S0 = arma::as_scalar(S);
    double s1 = sqrt(S0);
    double aa = (a0-mu0)/s1;
    double bb = (b0-mu0)/s1;
    M[0] = R::pnorm(bb,0,1,1,0) - R::pnorm(aa,0,1,1,0);
    if(kappa[0]>0){
      double pdfa = s1*R::dnorm(aa,0,1,0);
      double pdfb = s1*R::dnorm(bb,0,1,0);
      M(1) = mu0*M(0) + pdfa - pdfb;
      if(std::isinf(a0)){a0 = 0;} 
      if(std::isinf(b0)){b0 = 0;}
      for(int i=2; i<= kappa[0]; ++i){
        pdfa = pdfa*a0;
        pdfb = pdfb*b0;
        M(i) = mu0*M(i-1) + (i-1)*S0*M(i-2) + pdfa - pdfb;
      }
    }
    out["index"] = index;
    out["y"] = M;
    return out;
  }else{
    arma::uvec seqq = arma::linspace<arma::uvec>(1,n,n);
    int pk = prod(kappa+1);
    M.set_size(pk);
    M.zeros();
    // #
    // #  We create two long vectors G and H to store the two different sets
    // #  of integrals with dimension n-1.
    // #
    arma::vec nn = round(pk/(kappa+1));
    arma::vec begind = arma::cumsum(nn);
    begind = join_cols(arma::zeros<arma::vec>(1),begind);
    int pk1 = begind(n);
    // # number of (n-1)-dimensional integrals
    // #  Each row of cp corresponds to a vector that allows us to map the subscripts
    // #  to the index in the long vectors G and H
    arma::mat cp(n,n,arma::fill::zeros);
    for(int i=1; i<=n; ++i){
      arma::vec kk = kappa;
      kk(i-1) = 0;
      arma::vec cump = cumprod(kk(arma::span(0,n-2))+1);
      cp.row(i-1) = join_horiz(arma::ones<arma::rowvec>(1),cump.t());
    }
    arma::vec G(pk1,arma::fill::zeros);
    arma::vec H(pk1,arma::fill::zeros);
    arma::colvec s = sqrt(arma::diagvec(S));
    Rcpp::NumericVector ha = Rcpp::wrap((a-mu)/s);
    Rcpp::NumericVector hb = Rcpp::wrap((b-mu)/s);
    arma::mat R = S%(1/(s * s.t()));
    arma::colvec pdfa = Rcpp::dnorm(ha);
    pdfa = pdfa/s;
    arma::colvec pdfb = Rcpp::dnorm(hb);
    pdfb = pdfb/s;
    arma::uvec ii(1);
    for(int i=1; i<= n; ++i){
      ii = i;
      arma::uvec ind2 = seqq.elem(find(seqq != i));
      arma::vec kappai = kappa(ind2-1);
      arma::vec ai = a(ind2-1);
      arma::vec bi = b(ind2-1);
      arma::vec mui = mu(ind2-1);
      arma::mat Si = S(ind2-1,ii-1);
      arma::mat SSi = S(ind2-1,ind2-1)-Si*Si.t()/arma::as_scalar(S(ii-1,ii-1));
      arma::uvec ind = arma::regspace<arma::uvec>(begind(i-1)+1,1,begind(i));
      if(std::isfinite(a(i-1))){
        arma::vec mai = mui + Si/S(i-1,i-1)*(a(i-1)-mu(i-1));
        G(ind-1) = pdfa(i-1)*recintab0(kappai,ai,bi,mai,SSi);
      } 
      if(std::isfinite(b(i-1))){
        arma::vec mbi = mui + Si/S(i-1,i-1)*(b(i-1)-mu(i-1));
        H(ind-1) = pdfb(i-1)*recintab0(kappai,ai,bi,mbi,SSi);
      }
    }
    arma::vec ha1 = as<arma::vec>(ha);
    arma::vec hb1 = as<arma::vec>(hb);
    M(0) = pmvt_cpp(ha1,hb1,R);
    //M(0) = 1;
    a.elem(find_nonfinite(a)).zeros();
    b.elem(find_nonfinite(b)).zeros();
    arma::rowvec cp1 = cp.row(n-1);
    arma::mat grid = mygrid(kappa) + 1;
    for(int i=2; i<=pk; ++i){
      arma::rowvec kk = grid.row(i-1);
      //int xx = sum((kk-1)%cp1)+1;
      int i1 = min(seqq.elem(find(kk>1)));
      arma::rowvec kk1 = kk;
      kk1(i1-1)--;
      int ind3 = i - cp1(i1-1);
      M(i-1) = mu(i1-1)*M(ind3-1);
      for(int j=1; j<=n; ++j){
        int kk2 = kk1(j-1)-1;
        if(kk2 > 0){
          M(i-1) += S(i1-1,j-1)*kk2*M(ind3-cp1(j-1)-1);
        }
        int ind4 = begind(j-1) + sum(cp.row(j-1)%(kk1-1)) - cp(j-1,j-1)*kk2 + 1;
        M(i-1) += S(i1-1,j-1)*(pow(a(j-1),kk2)*G(ind4-1)-pow(b(j-1),kk2)*H(ind4-1));
      }
    }
    out["index"] = grid-1;
    out["y"] = M;
    return out;
  }
}

//*************************************************************************
//I-FUNCTION STUDENT T INTERNAL (NO GRID) (old version non-ordered)
// arma::vec ifunrec0(int k, arma::vec mu, arma::mat S, double nu, int k0){
//   int n = mu.size();
//   arma::vec y;
//   if(n==1){
//     double s = sqrt(arma::as_scalar(S));
//     double mu2 = arma::as_scalar(mu)/s;
//     if(k0==1){
//       y.set_size(1);
//       y(0) = R::pt(mu2,nu,1,0);
//     }else{
//       y.set_size(k+1);
//       y.fill(arma::datum::nan);
//       y(0) = R::pt(mu2,nu,1,0);
//       double c = sqrt(nu/M_PI)/2*tgamma((nu-1)/2)/tgamma(nu/2)*pow(1+pow(mu2,2)/nu,-(nu-1)/2);
//       if(k>0 && nu>1){
//         y(1) = mu2*y(0) + c;
//       }
//       for(int i=2; i<= k; ++i){
//         if(nu<=i){break;}
//         y(i) = ((nu+1-2*i)*mu2*y(i-1)+(i-1)*(pow(mu2,2)+nu)*y(i-2))/(nu-i);
//       }
//       y = y%Rcpp::as<arma::colvec>(vecpow(rep(s,k+1),seq(0,k)));
//     }
//     return y;
//   }
//   int n1 = n - k0;
//   if(k==0 || n1==0){
//     arma::colvec s = sqrt(arma::diagvec(S));
//     mu = mu/s;
//     arma::mat R = S%(1/(s * s.t()));
//     y.set_size(1);
//     Rcpp::NumericVector mu2 = Rcpp::wrap(mu);
//     Rcpp::NumericMatrix R2 = Rcpp::wrap(R);
//     y(0) = FmvtRcpp(mu2,R2,nu);
//     return y;
//   }
//   arma::mat Cp1(n1+k+1,1,arma::fill::ones);
//   arma::mat Cp2(n1+k+1,n1,arma::fill::zeros);
//   arma::mat C = join_horiz(Cp1,Cp2);
//   for(int j=2; j<= n1+1; ++j){
//     C(arma::span(j-1,n1+k),j-1) = arma::cumsum(C(arma::span(j-2,n1+k-1),j-2));
//   }
//   y.set_size(C(n1+k,n1));
//   //y.ones();
//   y.fill(arma::datum::nan);
//   arma::vec mus(n);
//   arma::colvec s = sqrt(arma::diagvec(S));
//   mus = mu/s;
//   arma::colvec cc = sqrt(nu)*tgamma((nu-1)/2)/(2*sqrt(M_PI)*tgamma(nu/2)*s)%pow(1+pow(mus,2)/nu,-(nu-1)/2);
//   arma::mat Si = inv_sympd(S);
//   //
//   //   J is a long vector with n blocks which is used to store J^i_{kappa} for different
//   //   values of kappa.  The first k0 blocks of J has binom(n1+k-1,n1) number of elements.
//   //   For i>k0, the i-th block of J has binom(n-i+k-1,n-i) number of elements.
//   //
//   int jlen = C(n1+k-1,n1-1) + k0*C(n1+k-1,n1);
//   arma::vec J = arma::zeros(jlen);
//   arma::vec jstart = arma::zeros(n); // jstart(i)+1 tells us the location of the first element of J^i_kappa
//   arma::uvec ind = arma::linspace<arma::uvec>(1,n,n);
//   arma::uvec ii(1);
//   arma::mat v(n-1,n-1);
//   arma::mat ss(n-1,n-1);
//   arma::vec mm(n-1);
//   for(int i=1; i<= n-1; ++i){
//     arma::uvec ind2 = ind.elem(find(ind != i));
//     ii = i;
//     v = S(ind2-1,ind2-1)-S(ind2-1,ii-1)*S(ii-1,ind2-1)/arma::as_scalar(S(ii-1,ii-1));
//     ss = (nu+pow(mus(i-1),2))/(nu-1)*v;
//     mm = mu(ind2-1)-S(ind2-1,ii-1)/arma::as_scalar(S(ii-1,ii-1))*mu(i-1);
//     if(i<=k0){
//       jstart(i) = jstart(i-1)+C(n1+k-1,n1);
//       J(arma::span(jstart(i-1),jstart(i)-1)) = cc(i-1)*ifunrec0(k-1,mm,ss,nu-1,k0-1);
//     }else{
//       jstart(i) = jstart(i-1)+C(n+k-i-1,n-i);
//       J(arma::span(jstart(i-1),jstart(i)-1)) = cc(i-1)*ifunrec0(k-1,mm,ss,nu-1,i-1);
//     }
//   }
//   ss = (nu+pow(mus(n-1),2))/(nu-1)*(S(arma::span(0,n-2),arma::span(0,n-2))-S(arma::span(0,n-2),n-1)*S(n-1,arma::span(0,n-2))/arma::as_scalar(S(n-1,n-1)));
//   mm = mu(arma::span(0,n-2))-S(arma::span(0,n-2),n-1)/arma::as_scalar(S(n-1,n-1))*mu(n-1);
//   J(jlen-1) = cc(n-1)*arma::as_scalar(ifunrec0(0,mm,ss,nu-1,0));
//   //  I_{0_n} and I_{0_n+e_n} [LINE 101]
//   arma::mat R = S%(1/(s * s.t()));
//   Rcpp::NumericVector mu2 = Rcpp::wrap(mus);
//   Rcpp::NumericMatrix R2 = Rcpp::wrap(R);
//   y(0) = FmvtRcpp(mu2,R2,nu);
//   //y(0) = 1;
//   arma::mat E(n,n,arma::fill::eye);
//   E = E.cols(k0,n-1);
//   arma::uvec jstartU = arma::conv_to<arma::uvec>::from(jstart);
//   if(nu>=1){
//     y(1) = mu(n-1)*y(0) + arma::as_scalar(S(n-1,arma::span(0,n-1))*J(jstartU));
//   }
//   //   This is to compute I_{0_{n-1},i_n} [LINE 118]
//   arma::vec m1 = mu(arma::span(0,n-2))-S(arma::span(0,n-2),n-1)/arma::as_scalar(S(n-1,n-1))*mu(n-1);
//   for(int in=2; in<= k; ++in){
//     if(in>=nu){break;}
//     y(in) = ((nu-2*in+1)*mu(n-1)*y(in-1)+(in-1)*(pow(mu(n-1),2)+nu*arma::as_scalar(S(n-1,n-1)))*y(in-2))/(nu-in) +
//       arma::as_scalar(S(n-1,arma::span(0,n-2))*J(jstartU(arma::span(0,n-2))+in-1))-(in-1)/(nu-in)*arma::as_scalar(S(n-1,n-1))*arma::as_scalar(m1.t()*J(jstartU(arma::span(0,n-2))+in-2));
//   }
//   double ar;
//   arma::irowvec oneone = arma::ones<arma::irowvec>(1);
//   arma::irowvec twotwo = arma::ones<arma::irowvec>(1);
//   twotwo = 2;
//   arma::urowvec rr = arma::ones<arma::urowvec>(1);
//   arma::rowvec inda(E.n_cols);
//   int ind1 = k+1;
//   arma::uvec ind0 = arma::linspace<arma::uvec>(1,n,n);
//   arma::vec r2;
//   for(int r=n-1; r >= max(NumericVector::create(2,k0+1)); --r){
//     //     r is the position of the first nonzero element of kappa [LINE 129]
//     //     i_r=1
//     //     Compute I_{0_{r-1},1,i_{r+1},...i_n} for i_{r+1}+...+i_n=0,...,k-1
//     rr = r;
//     arma::vec C1 = inv_sympd(S(arma::span(r,n-1),arma::span(r,n-1)))*S(arma::span(r,n-1),r-1);
//     m1 = mu(r-1)-C1.t()*mu(arma::span(r,n-1));
//     arma::vec r1 = S(arma::span(0,r-1),r-1)-S(arma::span(0,r-1),arma::span(r,n-1))*C1;
//     arma::irowvec ind = join_horiz(join_horiz(arma::zeros<arma::irowvec>(r-1),oneone),join_horiz(arma::zeros<arma::irowvec>(n-r-1),-oneone));
//     int si = 0;  // si=sum(ind), it is updated instead of computed every time [CHECKED]
//     while(ind(r)<k-1||si==0){
//       for(int i=n; i >= k0+1; --i){
//         if(si<k){
//           ind(i-1) = ind(i-1)+1;
//           si = si+1;
//           ind1 = ind1 +1;
//           break;
//         }else{
//           si = si-ind(i-1);
//           ind(i-1) = 0;
//         }
//       }
//       if(si<nu){
//         inda = ind(arma::span(k0,ind.n_elem-1)) - E.row(r-1);
//         y(ind1-1) = arma::as_scalar(m1*y(arma::as_scalar(colexindi(inda,k,C))-1));
//         for(int i=r+1; i<=n; ++i){
//           y(ind1-1) = arma::as_scalar(y(ind1-1)+C1(i-r-1)*y(colexindi(inda+E.row(i-1),k,C)-1));
//         }
//         for(int i=1; i<=r; ++i){
//           ii =  i;
//           arma::mat aux1 = join_horiz(arma::conv_to<arma::mat>::from(ii),arma::conv_to<arma::mat>::from(ind(arma::span(r,n-1))));
//           arma::uvec indb = colexindj(aux1,k,C,jstart);
//           y(ind1-1) = y(ind1-1)+arma::as_scalar(r1(i-1)*J(indb-1));
//         }
//       }
//     }
//     //   i_r>1
//     //   Compute I_{0_{r-1},i_r,i_{r+1},...i_n} for i_r+...+i_n=0,...,k and i_r>1
//     if(r==n-1){
//       ar = arma::as_scalar(S(r-1,r-1))/arma::as_scalar(S(n-1,n-1));
//     }else{
//       arma::mat Sa = inv_sympd(S(arma::span(r+1,n-1),arma::span(r+1,n-1)));
//       ar = arma::as_scalar(S(r-1,r-1)-S(r-1,arma::span(r+1,n-1))*Sa*S(arma::span(r+1,n-1),r-1))/arma::as_scalar(S(r,r)-S(r,arma::span(r+1,n-1))*Sa*S(arma::span(r+1,n-1),r));
//     }
//     arma::uvec ind2 = ind0.elem(find(ind0 ==r || (ind0 >= r+2 && ind0 <= n)));
//     arma::uvec ind3 = ind0.elem(find(ind0 <= r-1)); //index 1...r-1
//     arma::mat C2 = inv_sympd(S(ind2-1,ind2-1))*S(ind2-1,rr);
//     double m2 = mu(r)-arma::as_scalar(C2.t()*mu(ind2-1));
//     r2 = S(arma::span(0,r-2),r)-S(ind3-1,ind2-1)*C2;
//     ind = join_horiz(join_horiz(arma::zeros<arma::irowvec>(r-1),twotwo),join_horiz(arma::zeros<arma::irowvec>(n-r-1),-oneone));
//     si = 1;
//     while(ind(r-1)<k||(si==1&&k>1)){
//       for(int i=n; i >= k0+1; --i){
//         if(si<k){
//           ind(i-1) = ind(i-1)+1;
//           si = si+1;
//           ind1 = ind1 +1;
//           break;
//         }else{
//           si = si-ind(i-1);
//           ind(i-1) = 0;
//         }
//       }
//       if(si<nu){
//         inda = ind(arma::span(k0,ind.n_elem-1))-E.row(r-1);
//         arma::rowvec indb = inda-E.row(r-1)+E.row(r);
//         double ccc2 = arma::as_scalar(ind(r-1));
//         double ccc3 = arma::as_scalar(ind(r));
//         double cr = (ccc2-1)/(ccc3+1)*ar;
//         y(ind1-1) = arma::as_scalar(m1*y(colexindi(inda,k,C)-1) + cr*(y(colexindi(indb+E.row(r),k,C)-1)-m2*y(colexindi(indb,k,C)-1))+ C1(0)*y(colexindi(inda+E.row(r),k,C)-1) - cr*C2(0)*y(colexindi(indb+E.row(r-1),k,C)-1));
//         for(int i=r+2; i<=n; ++i){
//           y(ind1-1) = arma::as_scalar(y(ind1-1) + C1(i-r-1)*y(colexindi(inda+E.row(i-1),k,C)-1) - cr*C2(i-r-1)*y(colexindi(indb+E.row(i-1),k,C)-1));
//         }
//         for(int i=1; i<=r-1; ++i){
//           ii =  i;
//           arma::mat aux1 = join_horiz(arma::conv_to<arma::mat>::from(ii),arma::conv_to<arma::mat>::from(inda));
//           arma::mat aux2 = join_horiz(arma::conv_to<arma::mat>::from(ii),arma::conv_to<arma::mat>::from(indb));
//           arma::uvec indd = colexindj(aux1,k,C,jstart);
//           arma::uvec inde = colexindj(aux2,k,C,jstart);
//           y(ind1-1) = y(ind1-1)+arma::as_scalar(r1(i-1)*J(indd-1)-cr*r2(i-1)*J(inde-1));
//         }
//       }
//     }
//   }
//   //   Compute I_{1,i_2,...,i_n} for i_2+...+i_n=0,...,k-1
//   if(k0 == 0){
//     arma::vec sm = Si*mu;
//     arma::irowvec ind = join_horiz(oneone,join_horiz(arma::zeros<arma::irowvec>(n-2),-oneone));
//     int si = 0;
//     while(ind(1)<k-1||si==0){
//       for(int i=n; i>=1; --i){
//         if(si<k){
//           ind(i-1) = ind(i-1)+1;
//           si = si+1;
//           ind1 = ind1 +1;
//           break;
//         }else{
//           si = si-ind(i-1);
//           ind(i-1) = 0;
//         }
//       }
//       if(si<nu){
//         inda = ind - E.row(0);
//         arma::mat aux3 = join_horiz(arma::conv_to<arma::mat>::from(oneone),arma::conv_to<arma::mat>::from(inda(arma::span(1,n-1))));
//         arma::uvec indf = colexindj(aux3,k,C,jstart);
//         y(ind1-1) = arma::as_scalar(sm(0)*y(arma::as_scalar(colexindi(inda,k,C))-1)+J(indf-1));
//         for(int i=2; i<=n; ++i){
//           y(ind1-1) = arma::as_scalar(y(ind1-1)-Si(0,i-1)*y(colexindi(inda+E.row(i-1),k,C)-1));
//         }
//         y(ind1-1) = y(ind1-1)/Si(0,0);
//       }
//     }
//     //   Finally compute I_{i1,i2,...,in} for i1=2,...,k1 and i2+...+in=0,...,k-i1
//     ind = join_horiz(twotwo,join_horiz(arma::zeros<arma::irowvec>(n-2),-oneone));
//     si = 1;
//     while(ind(0)<k||(si==1&&k>1)){
//       for(int i=n; i >= 1; --i){
//         if(si<k){
//           ind(i-1) = ind(i-1)+1;
//           si = si+1;
//           ind1 = ind1 +1;
//           break;
//         }else{
//           si = si-ind(i-1);
//           ind(i-1) = 0;
//         }
//       }
//       if(si<nu){
//         inda = ind-E.row(0);
//         arma::rowvec indb = inda-E.row(0)+E.row(1);
//         double ccc0 = arma::as_scalar(ind(0));
//         double ccc1 = arma::as_scalar(ind(1));
//         double c = (ccc0-1)/(ccc1+1);
//         y(ind1-1) = arma::as_scalar(sm(0)*y(colexindi(inda,k,C)-1) + c*(Si(1,0)*y(colexindi(indb+E.row(0),k,C)-1)-sm(1)*y(colexindi(indb,k,C)-1)));
//         for(int i=2; i<=n; ++i){
//           y(ind1-1) = arma::as_scalar(y(ind1-1) - Si(0,i-1)*y(colexindi(inda+E.row(i-1),k,C)-1) + c*Si(1,i-1)*y(colexindi(indb+E.row(i-1),k,C)-1));
//         }
//         y(ind1-1) = y(ind1-1)/Si(0,0);
//       }
//     }
//   }
//   return y;
// }

//*************************************************************************
//I-FUNCTION STUDENT T INTERNAL (NO GRID)
// [[Rcpp::export]]
arma::vec ifunrec02(int k, arma::vec mu, arma::mat S, double nu){
  int n = mu.size();
  arma::vec y;
  if(n==1){
    double s = sqrt(arma::as_scalar(S));
    double mu2 = arma::as_scalar(mu)/s;
    y.set_size(k+1);
    y.fill(arma::datum::nan);
    y(0) = R::pt(mu2,nu,1,0);
    double c = sqrt(nu/M_PI)/2*tgamma((nu-1)/2)/tgamma(nu/2)*pow(1+pow(mu2,2)/nu,-(nu-1)/2);
    if(k>0 && nu>1){
      y(1) = mu2*y(0) + c;
    }
    for(int i=2; i<= k; ++i){
      if(nu<=i){break;}
      y(i) = ((nu+1-2*i)*mu2*y(i-1)+(i-1)*(pow(mu2,2)+nu)*y(i-2))/(nu-i);
    }
    
    NumericVector base1 = Rcpp::wrap(rep(s,k+1));
    IntegerVector exp1 = seq_len(k+1);
    exp1 = exp1 - 1;
    NumericVector exp2  = as<NumericVector>(exp1);
    
    
    y = y%Rcpp::as<arma::colvec>(vecpow(base1,exp2));
    return y;
  }
  arma::mat Cp1(n+k+1,1,arma::fill::ones);
  arma::mat Cp2(n+k+1,n,arma::fill::zeros);
  arma::mat C = join_horiz(Cp1,Cp2);
  for(int j=2; j<= n+1; ++j){
    C(arma::span(j-1,n+k),j-1) = arma::cumsum(C(arma::span(j-2,n+k-1),j-2));
  }
  y.set_size(C(n+k,n));
  y.fill(arma::datum::nan);
  arma::vec mus(n);
  arma::colvec s = sqrt(arma::diagvec(S));
  mus = mu/s;
  arma::mat R = S%(1/(s * s.t()));
  y(0) = Fpmvt_cpp(mus,R,nu);
  
  arma::mat Si = inv_sympd(S);
  
  if (k>0 && nu>1){
    //
    //   J is a long vector with n blocks which is used to store J^i_{kappa} for different
    //   values of kappa.  The first k0 blocks of J has binom(n1+k-1,n1) number of elements.
    //   For i>k0, the i-th block of J has binom(n-i+k-1,n-i) number of elements.
    //
    arma::colvec cc = sqrt(nu)*tgamma((nu-1)/2)/(2*sqrt(M_PI)*tgamma(nu/2)*s)%pow(1+pow(mus,2)/nu,-(nu-1)/2);
    int jlen = C(n+k-2,n-1); // nchoosek(n+k-2,n-1)
    arma::mat J = arma::zeros(jlen,n);
    arma::uvec ind = arma::linspace<arma::uvec>(1,n,n);//ok
    arma::uvec ii(1);//ok
    arma::mat v(n-1,n-1);
    arma::mat ss(n-1,n-1);
    arma::vec mm(n-1);
    for(int i=1; i<= n; ++i){
      arma::uvec ind2 = ind.elem(find(ind != i));
      ii = i;
      v = S(ind2-1,ind2-1)-S(ind2-1,ii-1)*S(ii-1,ind2-1)/arma::as_scalar(S(ii-1,ii-1));
      ss = (nu+pow(mus(i-1),2))/(nu-1)*v;
      mm = mu(ind2-1)-S(ind2-1,ii-1)/arma::as_scalar(S(ii-1,ii-1))*mu(i-1);
      J.col(i-1) = cc(i-1)*ifunrec02(k-1,mm,ss,nu-1);
    }
    arma::uvec indrev = arma::reverse<arma::uvec>(ind);
    arma::mat auxmat1 = J.row(0);
    y(indrev) = y(0)*mu + S*auxmat1.t();
    arma::mat ind1(n,n);
    arma::mat E(n,n,arma::fill::eye);//ok
    int start0 = 0;    // Start index of |kappa|-2
    int start1 = 1;    // Start index of |kappa|-1
    arma::mat ind0 = arma::zeros(1,n);
    arma::mat ind3 = arma::fliplr(E);
    arma::vec Q = arma::zeros(C(n+k-3,n-1));
    arma::vec c = arma::solve(S,mu);
    cc = nu + arma::as_scalar(mu.t()*c);
    //   This is to compute I_{0_{n-1},i_n}
    arma::colvec onesvec(n,arma::fill::ones);
    for(int in=2; in<= k; ++in){   //|kappa|=ii
      if(in>=nu){break;}
      //  Update Q
      int intaux1 = ind0.n_rows;
      for(int j=1; j<= intaux1; ++j){
        arma::uvec i1 = start1 + colexind2(onesvec*ind0.row(j-1) + E,C);
        Q(j-1) = arma::as_scalar((cc*y(start0+j-1) - arma::as_scalar(c.t()*y(i1-1))))/(nu-in);
      }
      ind0 = ind3;
      start0 = start1;
      start1 = start1 + ind3.n_rows;
      ind3 = colex(in,n);
      //  Update I_{kappa}
      int intaux2 = ind3.n_rows;
      for(int j=1; j<= intaux2; ++j){
        int i = min(ind.elem(find(ind3.row(j-1)>0)));   // first nonzero index
        arma::mat kappa = ind3.row(j-1) - E.row(i-1);
        arma::uvec i1 = start0 + colexind2(kappa,C);
        y(start1+j-1) = arma::as_scalar(mu(i-1)*y(i1-1));
        arma::uvec ll(1);
        for(int l=1; l<=n; ++l){
          ll=l;
          if(kappa(l-1)==0){
            arma::uvec indl = ind.elem(find(ind != l));
            arma::mat kappa2 = kappa(indl - 1);
            arma::uvec i2 = C(n+in-3,n-1) + colexind2(kappa2.t(),C);
            y(start1+j-1) += arma::as_scalar(S(i-1,l-1)*J(i2-1,ll-1));
          }else{
            y(start1+j-1) += arma::as_scalar(S(i-1,l-1)*kappa(l-1)*Q(colexind2(kappa - E.row(l-1),C)-1));
          }
        }
      }
    }
  }
  return y;
}

//*************************************************************************
//I-FUNCTION STUDENT T CASE (old version non-ordered)
// Rcpp::List ifunrec1(int k, arma::vec mu, arma::mat S, double nu, int k0){
//   Rcpp::List out;
//   int n = mu.size();
//   arma::vec y;
//   arma::mat index;
//   if(n==1){
//     double s = sqrt(arma::as_scalar(S));
//     double mu2 = arma::as_scalar(mu)/s;
//     if(k0==1){
//       y.set_size(1);
//       index.set_size(1,1);
//       index.zeros();
//       y(0) = R::pt(mu2,nu,1,0);
//     }else{
//       y.set_size(k+1);
//       y.fill(arma::datum::nan);
//       index = arma::linspace<arma::mat>(0,k,k+1);
//       y(0) = R::pt(mu2,nu,1,0);
//       double c = sqrt(nu/M_PI)/2*tgamma((nu-1)/2)/tgamma(nu/2)*pow(1+pow(mu2,2)/nu,-(nu-1)/2);
//       if(k>0 && nu>1){
//         y(1) = mu2*y(0) + c;
//       }
//       for(int i=2; i<= k; ++i){
//         if(nu<=i){break;}
//         y(i) = ((nu+1-2*i)*mu2*y(i-1)+(i-1)*(pow(mu2,2)+nu)*y(i-2))/(nu-i);
//       }
//       y = y%Rcpp::as<arma::colvec>(vecpow(rep(s,k+1),seq(0,k)));
//     }
//     out["index"] = index;
//     out["y"] = y;
//     return out;
//   }
//   int n1 = n - k0;
//   if(k==0 || n1==0){
//     arma::colvec s = sqrt(arma::diagvec(S));
//     mu = mu/s;
//     arma::mat R = S%(1/(s * s.t()));
//     y.set_size(1);
//     index.set_size(1,1);
//     index.zeros();
//     Rcpp::NumericVector mu2 = Rcpp::wrap(mu);
//     Rcpp::NumericMatrix R2 = Rcpp::wrap(R);
//     y(0) = FmvtRcpp(mu2,R2,nu);
//     out["index"] = index;
//     out["y"] = y;
//     return out;
//   }
//   arma::mat Cp1(n1+k+1,1,arma::fill::ones);
//   arma::mat Cp2(n1+k+1,n1,arma::fill::zeros);
//   arma::mat C = join_horiz(Cp1,Cp2);
//   for(int j=2; j<= n1+1; ++j){
//     C(arma::span(j-1,n1+k),j-1) = arma::cumsum(C(arma::span(j-2,n1+k-1),j-2));
//   }
//   y.set_size(C(n1+k,n1));
//   index.set_size(C(n1+k,n1),n);
//   index.zeros();
//   y.fill(arma::datum::nan);
//   arma::vec mus(n);
//   arma::colvec s = sqrt(arma::diagvec(S));
//   mus = mu/s;
//   arma::colvec cc = sqrt(nu)*tgamma((nu-1)/2)/(2*sqrt(M_PI)*tgamma(nu/2)*s)%pow(1+pow(mus,2)/nu,-(nu-1)/2);
//   arma::mat Si = inv_sympd(S);
//   //
//   //   J is a long vector with n blocks which is used to store J^i_{kappa} for different
//   //   values of kappa.  The first k0 blocks of J has binom(n1+k-1,n1) number of elements.
//   //   For i>k0, the i-th block of J has binom(n-i+k-1,n-i) number of elements.
//   //
//   int jlen = C(n1+k-1,n1-1) + k0*C(n1+k-1,n1);
//   arma::vec J = arma::zeros(jlen);
//   arma::vec jstart = arma::zeros(n); // jstart(i)+1 tells us the location of the first element of J^i_kappa
//   arma::uvec ind = arma::linspace<arma::uvec>(1,n,n);
//   arma::uvec ii(1);
//   arma::mat v(n-1,n-1);
//   arma::mat ss(n-1,n-1);
//   arma::vec mm(n-1);
//   for(int i=1; i<= n-1; ++i){
//     arma::uvec ind2 = ind.elem(find(ind != i));
//     ii = i;
//     v = S(ind2-1,ind2-1)-S(ind2-1,ii-1)*S(ii-1,ind2-1)/arma::as_scalar(S(ii-1,ii-1));
//     ss = (nu+pow(mus(i-1),2))/(nu-1)*v;
//     mm = mu(ind2-1)-S(ind2-1,ii-1)/arma::as_scalar(S(ii-1,ii-1))*mu(i-1);
//     if(i<=k0){
//       jstart(i) = jstart(i-1)+C(n1+k-1,n1);
//       J(arma::span(jstart(i-1),jstart(i)-1)) = cc(i-1)*ifunrec0(k-1,mm,ss,nu-1,k0-1);
//     }else{
//       jstart(i) = jstart(i-1)+C(n+k-i-1,n-i);
//       J(arma::span(jstart(i-1),jstart(i)-1)) = cc(i-1)*ifunrec0(k-1,mm,ss,nu-1,i-1);
//     }
//   }
//   ss = (nu+pow(mus(n-1),2))/(nu-1)*(S(arma::span(0,n-2),arma::span(0,n-2))-S(arma::span(0,n-2),n-1)*S(n-1,arma::span(0,n-2))/arma::as_scalar(S(n-1,n-1)));
//   mm = mu(arma::span(0,n-2))-S(arma::span(0,n-2),n-1)/arma::as_scalar(S(n-1,n-1))*mu(n-1);
//   J(jlen-1) = cc(n-1)*arma::as_scalar(ifunrec0(0,mm,ss,nu-1,0));
//   //  I_{0_n} and I_{0_n+e_n} [LINE 101]
//   arma::mat R = S%(1/(s * s.t()));
//   Rcpp::NumericVector mu2 = Rcpp::wrap(mus);
//   Rcpp::NumericMatrix R2 = Rcpp::wrap(R);
//   y(0) = FmvtRcpp(mu2,R2,nu);
//   arma::mat E(n,n,arma::fill::eye);
//   E = E.cols(k0,n-1);
//   arma::uvec jstartU = arma::conv_to<arma::uvec>::from(jstart);
//   if(nu>=1){
//     index(1,n-1) = 1; //PRINTING INDEX
//     y(1) = mu(n-1)*y(0) + arma::as_scalar(S(n-1,arma::span(0,n-1))*J(jstartU));
//   }
//   //   This is to compute I_{0_{n-1},i_n} [LINE 118]
//   arma::vec m1 = mu(arma::span(0,n-2))-S(arma::span(0,n-2),n-1)/arma::as_scalar(S(n-1,n-1))*mu(n-1);
//   for(int in=2; in<= k; ++in){
//     if(in>=nu){break;}
//     index(in,n-1) = in; //PRINTING INDEX
//     y(in) = ((nu-2*in+1)*mu(n-1)*y(in-1)+(in-1)*(pow(mu(n-1),2)+nu*arma::as_scalar(S(n-1,n-1)))*y(in-2))/(nu-in) +
//       arma::as_scalar(S(n-1,arma::span(0,n-2))*J(jstartU(arma::span(0,n-2))+in-1))-(in-1)/(nu-in)*arma::as_scalar(S(n-1,n-1))*arma::as_scalar(m1.t()*J(jstartU(arma::span(0,n-2))+in-2));
//   }
//   double ar;
//   arma::irowvec oneone = arma::ones<arma::irowvec>(1);
//   arma::irowvec twotwo = arma::ones<arma::irowvec>(1);
//   twotwo = 2;
//   arma::urowvec rr = arma::ones<arma::urowvec>(1);
//   arma::rowvec inda(E.n_cols);
//   int ind1 = k+1;
//   arma::uvec ind0 = arma::linspace<arma::uvec>(1,n,n);
//   arma::vec r2;
//   for(int r=n-1; r >= max(NumericVector::create(2,k0+1)); --r){
//     //     r is the position of the first nonzero element of kappa [LINE 129]
//     //     i_r=1
//     //     Compute I_{0_{r-1},1,i_{r+1},...i_n} for i_{r+1}+...+i_n=0,...,k-1
//     rr = r;
//     arma::vec C1 = inv_sympd(S(arma::span(r,n-1),arma::span(r,n-1)))*S(arma::span(r,n-1),r-1);
//     m1 = mu(r-1)-C1.t()*mu(arma::span(r,n-1));
//     arma::vec r1 = S(arma::span(0,r-1),r-1)-S(arma::span(0,r-1),arma::span(r,n-1))*C1;
//     arma::irowvec ind = join_horiz(join_horiz(arma::zeros<arma::irowvec>(r-1),oneone),join_horiz(arma::zeros<arma::irowvec>(n-r-1),-oneone));
//     arma::uvec ind0 = arma::linspace<arma::uvec>(1,n,n);
//     int si = 0;  // si=sum(ind), it is updated instead of computed every time [CHECKED]
//     while(ind(r)<k-1||si==0){
//       for(int i=n; i >= k0+1; --i){
//         if(si<k){
//           ind(i-1) = ind(i-1)+1;
//           si = si+1;
//           ind1 = ind1 +1;
//           break;
//         }else{
//           si = si-ind(i-1);
//           ind(i-1) = 0;
//         }
//       }
//       if(si<nu){
//         index.row(ind1-1) = arma::conv_to<arma::mat>::from(ind(arma::span(k0,ind.n_elem-1)));  //PRINTING INDEX
//         inda = ind(arma::span(k0,ind.n_elem-1)) - E.row(r-1);
//         y(ind1-1) = arma::as_scalar(m1*y(arma::as_scalar(colexindi(inda,k,C))-1));
//         for(int i=r+1; i<=n; ++i){
//           y(ind1-1) = arma::as_scalar(y(ind1-1)+C1(i-r-1)*y(colexindi(inda+E.row(i-1),k,C)-1));
//         }
//         for(int i=1; i<=r; ++i){
//           ii =  i;
//           arma::mat aux1 = join_horiz(arma::conv_to<arma::mat>::from(ii),arma::conv_to<arma::mat>::from(ind(arma::span(r,n-1))));
//           arma::uvec indb = colexindj(aux1,k,C,jstart);
//           y(ind1-1) = y(ind1-1)+arma::as_scalar(r1(i-1)*J(indb-1));
//         }
//       }
//     }
//     //   i_r>1
//     //   Compute I_{0_{r-1},i_r,i_{r+1},...i_n} for i_r+...+i_n=0,...,k and i_r>1
//     if(r==n-1){
//       ar = arma::as_scalar(S(r-1,r-1))/arma::as_scalar(S(n-1,n-1));
//     }else{
//       arma::mat Sa = inv_sympd(S(arma::span(r+1,n-1),arma::span(r+1,n-1)));
//       ar = arma::as_scalar(S(r-1,r-1)-S(r-1,arma::span(r+1,n-1))*Sa*S(arma::span(r+1,n-1),r-1))/arma::as_scalar(S(r,r)-S(r,arma::span(r+1,n-1))*Sa*S(arma::span(r+1,n-1),r));
//     }
//     arma::uvec ind2 = ind0.elem(find(ind0 ==r || (ind0 >= r+2 && ind0 <= n)));
//     arma::uvec ind3 = ind0.elem(find(ind0 <= r-1)); //index 1...r-1
//     arma::mat C2 = inv_sympd(S(ind2-1,ind2-1))*S(ind2-1,rr);
//     double m2 = mu(r)-arma::as_scalar(C2.t()*mu(ind2-1));
//     r2 = S(arma::span(0,r-2),r)-S(ind3-1,ind2-1)*C2;
//     ind = join_horiz(join_horiz(arma::zeros<arma::irowvec>(r-1),twotwo),join_horiz(arma::zeros<arma::irowvec>(n-r-1),-oneone));
//     si = 1;
//     while(ind(r-1)<k||(si==1&&k>1)){
//       for(int i=n; i >= k0+1; --i){
//         if(si<k){
//           ind(i-1) = ind(i-1)+1;
//           si = si+1;
//           ind1 = ind1 +1;
//           break;
//         }else{
//           si = si-ind(i-1);
//           ind(i-1) = 0;
//         }
//       }
//       if(si<nu){
//         inda = ind(arma::span(k0,ind.n_elem-1))-E.row(r-1);
//         arma::rowvec indb = inda-E.row(r-1)+E.row(r);
//         double ccc2 = arma::as_scalar(ind(r-1));
//         double ccc3 = arma::as_scalar(ind(r));
//         double cr = (ccc2-1)/(ccc3+1)*ar;
//         index.row(ind1-1) = arma::conv_to<arma::mat>::from(ind(arma::span(k0,ind.n_elem-1)));  //PRINTING INDEX
//         y(ind1-1) = arma::as_scalar(m1*y(colexindi(inda,k,C)-1) + cr*(y(colexindi(indb+E.row(r),k,C)-1)-m2*y(colexindi(indb,k,C)-1))+ C1(0)*y(colexindi(inda+E.row(r),k,C)-1) - cr*C2(0)*y(colexindi(indb+E.row(r-1),k,C)-1));
//         for(int i=r+2; i<=n; ++i){
//           y(ind1-1) = arma::as_scalar(y(ind1-1) + C1(i-r-1)*y(colexindi(inda+E.row(i-1),k,C)-1) - cr*C2(i-r-1)*y(colexindi(indb+E.row(i-1),k,C)-1));
//         }
//         for(int i=1; i<=r-1; ++i){
//           ii =  i;
//           arma::mat aux1 = join_horiz(arma::conv_to<arma::mat>::from(ii),arma::conv_to<arma::mat>::from(inda));
//           arma::mat aux2 = join_horiz(arma::conv_to<arma::mat>::from(ii),arma::conv_to<arma::mat>::from(indb));
//           arma::uvec indd = colexindj(aux1,k,C,jstart);
//           arma::uvec inde = colexindj(aux2,k,C,jstart);
//           y(ind1-1) = y(ind1-1)+arma::as_scalar(r1(i-1)*J(indd-1)-cr*r2(i-1)*J(inde-1));
//         }
//       }
//     }
//   }
//   //   Compute I_{1,i_2,...,i_n} for i_2+...+i_n=0,...,k-1
//   if(k0 == 0){
//     arma::vec sm = Si*mu;
//     arma::irowvec ind = join_horiz(oneone,join_horiz(arma::zeros<arma::irowvec>(n-2),-oneone));
//     int si = 0;
//     while(ind(1)<k-1||si==0){
//       for(int i=n; i>=1; --i){
//         if(si<k){
//           ind(i-1) = ind(i-1)+1;
//           si = si+1;
//           ind1 = ind1 +1;
//           break;
//         }else{
//           si = si-ind(i-1);
//           ind(i-1) = 0;
//         }
//       }
//       
//       if(si<nu){
//         inda = ind - E.row(0);
//         arma::mat aux3 = join_horiz(arma::conv_to<arma::mat>::from(oneone),arma::conv_to<arma::mat>::from(inda(arma::span(1,n-1))));
//         arma::uvec indf = colexindj(aux3,k,C,jstart);
//         index.row(ind1-1) = arma::conv_to<arma::mat>::from(ind);  //PRINTING INDEX
//         y(ind1-1) = arma::as_scalar(sm(0)*y(arma::as_scalar(colexindi(inda,k,C))-1)+J(indf-1));
//         for(int i=2; i<=n; ++i){
//           y(ind1-1) = arma::as_scalar(y(ind1-1)-Si(0,i-1)*y(colexindi(inda+E.row(i-1),k,C)-1));
//         }
//         y(ind1-1) = y(ind1-1)/Si(0,0);
//       }
//     }
//     //   Finally compute I_{i1,i2,...,in} for i1=2,...,k1 and i2+...+in=0,...,k-i1
//     ind = join_horiz(twotwo,join_horiz(arma::zeros<arma::irowvec>(n-2),-oneone));
//     si = 1;
//     while(ind(0)<k||(si==1&&k>1)){
//       for(int i=n; i >= 1; --i){
//         if(si<k){
//           ind(i-1) = ind(i-1)+1;
//           si = si+1;
//           ind1 = ind1 +1;
//           break;
//         }else{
//           si = si-ind(i-1);
//           ind(i-1) = 0;
//         }
//       }
//       if(si<nu){
//         inda = ind-E.row(0);
//         arma::rowvec indb = inda-E.row(0)+E.row(1);
//         double ccc0 = arma::as_scalar(ind(0));
//         double ccc1 = arma::as_scalar(ind(1));
//         double c = (ccc0-1)/(ccc1+1);
//         index.row(ind1-1) = arma::conv_to<arma::mat>::from(ind);  //PRINTING INDEX
//         y(ind1-1) = arma::as_scalar(sm(0)*y(colexindi(inda,k,C)-1) + c*(Si(1,0)*y(colexindi(indb+E.row(0),k,C)-1)-sm(1)*y(colexindi(indb,k,C)-1)));
//         for(int i=2; i<=n; ++i){
//           y(ind1-1) = arma::as_scalar(y(ind1-1) - Si(0,i-1)*y(colexindi(inda+E.row(i-1),k,C)-1) + c*Si(1,i-1)*y(colexindi(indb+E.row(i-1),k,C)-1));
//         }
//         y(ind1-1) = y(ind1-1)/Si(0,0);
//       }
//     }
//   }
//   out["index"] = index;
//   out["y"] = y;
//   return out;
// }

//*************************************************************************
//I-FUNCTION STUDENT T CASE
Rcpp::List ifunrec2(int k, arma::vec mu, arma::mat S, double nu){
  Rcpp::List out;
  int n = mu.size();
  arma::vec y;
  arma::mat index;
  if(n==1){
    double s = sqrt(arma::as_scalar(S));
    double mu2 = arma::as_scalar(mu)/s;
    y.set_size(k+1);
    y.fill(arma::datum::nan);
    index = arma::linspace<arma::mat>(0,k,k+1);
    y(0) = R::pt(mu2,nu,1,0);
    double c = sqrt(nu/M_PI)/2*tgamma((nu-1)/2)/tgamma(nu/2)*pow(1+pow(mu2,2)/nu,-(nu-1)/2);
    if(k>0 && nu>1){
      y(1) = mu2*y(0) + c;
    }
    for(int i=2; i<= k; ++i){
      if(nu<=i){break;}
      y(i) = ((nu+1-2*i)*mu2*y(i-1)+(i-1)*(pow(mu2,2)+nu)*y(i-2))/(nu-i);
    }
    
    NumericVector base1 = Rcpp::wrap(rep(s,k+1));
    IntegerVector exp1 = seq_len(k+1);
    exp1 = exp1 - 1;
    NumericVector exp2  = as<NumericVector>(exp1);
    
    y = y%Rcpp::as<arma::colvec>(vecpow(base1,exp2));
    
    //y = y%Rcpp::as<arma::colvec>(vecpow(rep(s,k+1),seq(0,k)));
    
    out["index"] = index;
    out["y"] = y;
    return out;
  }
  arma::mat Cp1(n+k+1,1,arma::fill::ones);
  arma::mat Cp2(n+k+1,n,arma::fill::zeros);
  arma::mat C = join_horiz(Cp1,Cp2);
  for(int j=2; j<= n+1; ++j){
    C(arma::span(j-1,n+k),j-1) = arma::cumsum(C(arma::span(j-2,n+k-1),j-2));
  }
  y.set_size(C(n+k,n));
  y.fill(arma::datum::nan);
  index.set_size(1,n);
  index.zeros();
  arma::vec mus(n);
  arma::colvec s = sqrt(arma::diagvec(S));
  mus = mu/s;
  arma::mat R = S%(1/(s * s.t()));
  y(0) = Fpmvt_cpp(mus,R,nu);
  
  arma::mat Si = inv_sympd(S);
  
  if (k>0 && nu>1){
    //
    //   J is a long vector with n blocks which is used to store J^i_{kappa} for different
    //   values of kappa.  The first k0 blocks of J has binom(n1+k-1,n1) number of elements.
    //   For i>k0, the i-th block of J has binom(n-i+k-1,n-i) number of elements.
    //
    arma::colvec cc = sqrt(nu)*tgamma((nu-1)/2)/(2*sqrt(M_PI)*tgamma(nu/2)*s)%pow(1+pow(mus,2)/nu,-(nu-1)/2);
    int jlen = C(n+k-2,n-1); // nchoosek(n+k-2,n-1)
    arma::mat J = arma::zeros(jlen,n);
    arma::uvec ind = arma::linspace<arma::uvec>(1,n,n);//ok
    arma::uvec ii(1);//ok
    arma::mat v(n-1,n-1);
    arma::mat ss(n-1,n-1);
    arma::vec mm(n-1);
    for(int i=1; i<= n; ++i){
      arma::uvec ind2 = ind.elem(find(ind != i));
      ii = i;
      v = S(ind2-1,ind2-1)-S(ind2-1,ii-1)*S(ii-1,ind2-1)/arma::as_scalar(S(ii-1,ii-1));
      ss = (nu+pow(mus(i-1),2))/(nu-1)*v;
      mm = mu(ind2-1)-S(ind2-1,ii-1)/arma::as_scalar(S(ii-1,ii-1))*mu(i-1);
      J.col(i-1) = cc(i-1)*ifunrec02(k-1,mm,ss,nu-1);
    }
    arma::uvec indrev = arma::reverse<arma::uvec>(ind);
    arma::mat auxmat1 = J.row(0);
    y(indrev) = y(0)*mu + S*auxmat1.t();
    arma::mat ind1(n,n);
    arma::mat E(n,n,arma::fill::eye);//ok
    int start0 = 0;    // Start index of |kappa|-2
    int start1 = 1;    // Start index of |kappa|-1
    arma::mat ind0 = arma::zeros(1,n);
    arma::mat ind3 = arma::fliplr(E);
    index = join_vert(index,ind3);
    arma::vec Q = arma::zeros(C(n+k-3,n-1));
    arma::vec c = arma::solve(S,mu);
    cc = nu + arma::as_scalar(mu.t()*c);
    //   This is to compute I_{0_{n-1},i_n}
    arma::colvec onesvec(n,arma::fill::ones);
    for(int in=2; in<= k; ++in){   //|kappa|=ii
      if(in>=nu){break;}
      //  Update Q
      int intaux1 = ind0.n_rows;
      for(int j=1; j<= intaux1; ++j){
        arma::uvec i1 = start1 + colexind2(onesvec*ind0.row(j-1) + E,C);
        Q(j-1) = arma::as_scalar((cc*y(start0+j-1) - arma::as_scalar(c.t()*y(i1-1))))/(nu-in);
      }
      ind0 = ind3;
      start0 = start1;
      start1 = start1 + ind3.n_rows;
      ind3 = colex(in,n);
      index = join_vert(index,ind3);
      //  Update I_{kappa}
      int intaux2 = ind3.n_rows;
      for(int j=1; j<= intaux2; ++j){
        int i = min(ind.elem(find(ind3.row(j-1)>0)));   // first nonzero index
        arma::mat kappa = ind3.row(j-1) - E.row(i-1);
        arma::uvec i1 = start0 + colexind2(kappa,C);
        y(start1+j-1) = arma::as_scalar(mu(i-1)*y(i1-1));
        arma::uvec ll(1);
        for(int l=1; l<=n; ++l){
          ll=l;
          if(kappa(l-1)==0){
            arma::uvec indl = ind.elem(find(ind != l));
            arma::mat kappa2 = kappa(indl - 1);
            arma::uvec i2 = C(n+in-3,n-1) + colexind2(kappa2.t(),C);
            y(start1+j-1) += arma::as_scalar(S(i-1,l-1)*J(i2-1,ll-1));
          }else{
            y(start1+j-1) += arma::as_scalar(S(i-1,l-1)*kappa(l-1)*Q(colexind2(kappa - E.row(l-1),C)-1));
          }
        }
      }
    }
  }
  out["index"] = index;
  out["y"] = y;
  return out;
}

//*************************************************************************
//*************************** user functions ******************************
//*************************************************************************

//MEANVAR -> FOLDED -> NORMAL
// [[Rcpp::export]]
Rcpp::List RcppmeanvarFN(arma::vec mu,arma::mat S){
  int n = mu.size();
  arma::colvec muY(n);
  arma::mat varY(n,n);
  muY.fill(arma::datum::nan);
  varY.fill(arma::datum::nan);
  Rcpp::List out;
  arma::colvec s = sqrt(arma::diagvec(S));
  Rcpp::NumericVector h = Rcpp::wrap(mu/s);
  arma::colvec pdfh = Rcpp::dnorm(h);
  arma::colvec cdfh = Rcpp::pnorm(h);
  double dos = 2;
  arma::colvec erfv = erf(Rcpp::as<arma::mat>(h)/sqrt(dos));
  muY = s%(2*pdfh + Rcpp::as<arma::mat>(h)%erfv);
  arma::mat R = S%(1/(s * s.t()));
  
  arma::rowvec jv(n,arma::fill::ones);
  arma::mat h1 =  Rcpp::as<arma::mat>(h)*jv;
  arma::mat A = (h1-R%h1.t())%(1/sqrt(2*(1 - R%R)));
  A.diag().zeros();
  arma::mat gam = (Rcpp::as<arma::mat>(h)*pdfh.t())%erf(A);
  
  arma::vec hij(2);
  arma::mat W(2,2,arma::fill::eye);
  for(int i=0; i<n; ++i){
    varY(i,i) = pow(mu(i),2) + pow(s(i),2);
    for(int j=0; j<i; ++j){
      double r   = R(i,j);
      double eta = sqrt(1 - pow(r,2));
      W(0,1) = W(1,0) = r;
      hij(0) = h(i);
      hij(1) = h(j);
      double normmt = Fpmvt_cpp(hij,W);
      
      
      double p  = 4*normmt - 2*cdfh(i) - 2*cdfh(j) + 1;
      
      double c  = sqrt(pow(h(i),2) + pow(h(j),2) - 2*r*h(i)*h(j))/eta;
      
      varY(i,j) = s(i)*s(j)*(p*(h(i)*h(j)+r)+2*gam(i,j)+2*gam(j,i)+4*eta/sqrt(2*M_PI)*R::dnorm(c,0,1,0));
      varY(j,i) = varY(i,j);
    }
  }
  arma::mat EYY(n,n);
  EYY = varY;
  varY = EYY - muY * muY.t();
  out["mean"] = muY;
  out["EYY"] = EYY;
  out["varcov"] = varY;
  return out;
}

//*************************************************************************
//MOMENTS -> FOLDED -> NORMAL
// [[Rcpp::export]]
Rcpp::List RcppmomentsFN(arma::vec kappa, arma::vec mu, arma::mat S){
  Rcpp::List out;
  int n = mu.size();
  arma::vec a(n,arma::fill::zeros);
  arma::vec b(n);
  b.fill(arma::datum::inf);
  out = recintab1(kappa,a,b,mu,S);
  arma::vec M = out[1];
  double dos = 2;
  double ene = n;
  for(int i=0; i<=pow(dos,ene)-2; ++i){
    arma::rowvec ind1 = DecToSigns(i,n);
    M += recintab0(kappa,a,b,mu%ind1.t(),S%(ind1.t()*ind1));
  }
  M[0] = 1;
  out[1] = M;
  return out;
}

//*************************************************************************
//MEANVAR -> FOLDED -> STUDENT T
// [[Rcpp::export]]
Rcpp::List RcppmeanvarFT(arma::vec mu,arma::mat S,double nu){
  int n = mu.size();
  arma::colvec muY(n);
  arma::mat EYY(n,n);
  arma::mat varY(n,n);
  muY.fill(arma::datum::nan);
  EYY.fill(arma::datum::nan);
  varY.fill(arma::datum::nan);
  Rcpp::List out;
  if(nu <= 1) return out;
  arma::colvec s = sqrt(arma::diagvec(S));
  Rcpp::NumericVector h = Rcpp::wrap(mu/s);
  arma::colvec c0 = sqrt(nu)*tgamma((nu-1)/2)/sqrt(M_PI)/tgamma(nu/2)*pow(1+pow(h,2)/nu,-(nu-1)/2);
  arma::colvec c1 = 2*Rcpp::pt(h,nu)- 1;
  muY = s%(c0 +  Rcpp::as<arma::mat>(h)%c1);
  if(nu<=2){
    out["mean"] = muY;
    out["EYY"] = EYY;
    out["varcov"] = varY;
    return out;
  }
  arma::mat R = S%(1/(s * s.t()));
  arma::vec hij(2);
  arma::mat W(2,2,arma::fill::eye);
  W(0,0) = W(1,1) = 1;
  for(int i=0; i<n; ++i){
    varY(i,i) = pow(mu(i),2) + pow(s(i),2)*nu/(nu-2);
    for(int j=0; j<i; ++j){
      double r   = R(i,j);
      double eta = sqrt(1 - pow(r,2));
      W(0,1) = W(1,0) = r;
      hij(0) = h(i);
      hij(1) = h(j);
      double pmt = Fpmvt_cpp(hij,W,nu);
      
      double p   = 4*pmt - c1(i) - c1(j) - 1;
      double zij = sqrt(nu-1)*(h(i)-r*h(j))/eta/sqrt(nu+pow(h(j),2));
      double zji = sqrt(nu-1)*(h(j)-r*h(i))/eta/sqrt(nu+pow(h(i),2));
      double cc  = pow(1+(pow(h(i),2)+pow(h(j),2)-2*r*h(i)*h(j))/pow(eta,2)/nu,-(nu-2)/2);
      double r1  = r/(nu-2);
      double v1 =  R::pt(zij,nu-1,1,0);
      double v2 =  R::pt(zji,nu-1,1,0);
      varY(i,j) = s(i)*s(j)*(p*(h(i)*h(j)+nu*r1)+(h(i)-r1*h(j))*c0(j)*(2*v1-1)+(h(j)-r1*h(i))*c0(i)*(2*v2-1)+2*nu/(nu-2)*eta/M_PI*cc);
      varY(j,i) = varY(i,j);
    }
  }
  EYY = varY;
  varY = EYY - muY * muY.t();
  out["mean"] = muY;
  out["EYY"] = EYY;
  out["varcov"] = varY;
  return out;
}

//*************************************************************************
//MOMENTS -> FOLDED -> STUDENT T (non-ordered)
// // [[Rcpp::export]]
// Rcpp::List momentsFT(int k, arma::vec mu, arma::mat S, double nu){
//   Rcpp::List out;
//   int n = mu.size();
//   out = ifunrec2(k,mu,S,nu);
//   arma::vec M = out[1];
//   for(int i=0; i<=pow(2,n)-2; ++i){
//     arma::rowvec ind1 = DecToSigns(i,n);
//     M += ifunrec02(k,mu%ind1.t(),S%(ind1.t()*ind1),nu);
//   }
//   M[0] = 1;
//   out[1] = M;
//   return out;
// }

//*************************************************************************
//INVERTING ORDER
// [[Rcpp::export]]
arma::mat Rcpporder(int k,arma::mat W){
  int l;
  int counter = 1;
  int n = W.n_cols - 1;
  for(int j=1; j<= k; ++j){
    l = binom(n+j-1,n-1);
    arma::mat TEMP = W.rows(counter+1-1,counter+l-1);
    W.rows(counter+1-1,counter+l-1) = flipud(TEMP);
    counter += l;
  }
  return W;
}

//*************************************************************************
//MOMENTS -> FOLDED -> STUDENT T
// [[Rcpp::export]]
arma::mat RcppmomentsFT(int k, arma::vec mu, arma::mat S, double nu){
  Rcpp::List out;
  int n = mu.size();
  out = ifunrec2(k,mu,S,nu);
  arma::vec M = out[1];
  double dos = 2;
  double ene = n;
  for(int i=0; i<=pow(dos,ene)-2; ++i){
    arma::rowvec ind1 = DecToSigns(i,n);
    M += ifunrec02(k,mu%ind1.t(),S%(ind1.t()*ind1),nu);
  }
  M[0] = 1;
  arma::mat L0 = out[0];
  arma::vec L1 = M;
  arma::mat W = join_horiz(L0,L1);
  int l;
  int counter = 1;
  for(int j=1; j<= k; ++j){
    l = binom(n+j-1,n-1);
    arma::mat TEMP = W.rows(counter+1-1,counter+l-1);
    W.rows(counter+1-1,counter+l-1) = flipud(TEMP);
    counter += l;
  }
  return W;
}

//*************************************************************************
//*************************************************************************
//MEANVAR -> TRUNCATED -> NORMAL
// [[Rcpp::export]]
Rcpp::List RcppmeanvarN_b(arma::vec b,arma::vec mu,arma::mat S, double p){
  Rcpp::List out;
  int n = mu.n_elem;
  arma::colvec s = sqrt(arma::diagvec(S));
  if(n==1){
    double b0 = arma::as_scalar(b);
    double mu0 = arma::as_scalar(mu);
    double S0 = arma::as_scalar(S);
    double s1 = sqrt(S0);
    double b1 = (b0-mu0)/s1;
    double p = R::pnorm(b1,0,1,1,0);
    double muY = mu0 - R::dnorm(b1,0,1,0)*s1/p;
    if(std::isinf(b0)){b0 = 0;}
    double varY = S0 +(mu0-muY)*muY-b0*R::dnorm(b1,0,1,0)*s1/p;
    double EYY = varY + pow(muY,2);
    out["mean"] = muY;
    out["EYY"] = EYY;
    out["varcov"] = varY;
    return out;
  }else{
    arma::vec muY(n);
    arma::mat varY(n,n);
    arma::uvec seqq = arma::linspace<arma::uvec>(1,n,n);
    arma::vec b1 = b-mu;
    arma::vec hb = b1/s;
    arma::mat R = S%(1/(s * s.t()));
    //double p = Fpmvt_cpp(hb,R);
    
    //double p = 1;
    arma::vec qb = Rcppqfun_b(b1,S);
    muY = mu + S*qb/p;
    arma::mat D(n,n,arma::fill::zeros);
    arma::uvec ii(1);
    double db;
    arma::vec mbi;
    for(int i=1; i<=n; ++i){
      arma::vec wb(n-1,arma::fill::zeros);
      ii = i;
      arma::uvec ind2 = seqq.elem(find(seqq != i));
      arma::mat SSi = S(ind2-1,ind2-1)-S(ind2-1,ii-1)*S(ii-1,ind2-1)/arma::as_scalar(S(ii-1,ii-1));
      D(i-1,i-1) -= b(i-1)*qb(i-1);
      mbi = mu(ind2-1) + S(ind2-1,ii-1)/S(i-1,i-1)*b1(i-1);
      arma::vec qb2 = Rcppqfun_b(b(ind2-1)-mbi,SSi);
      db = R::dnorm(b(i-1),mu(i-1),s(i-1),0);
      wb = qb(i-1)*mbi - db*SSi*qb2;
      D(ii-1,ind2-1) = -wb.t();
    }
    varY = S + S*(D - qb*muY.t())/p;
    varY = (varY + varY.t())/2;
    
    arma::mat EYY(n,n);
    EYY = varY + muY * muY.t();
    out["mean"] = muY;
    out["EYY"] = EYY;
    out["varcov"] = varY;
    
    return out;
  }
}

//*************************************************************************
//MEANVAR -> TRUNCATED -> NORMAL
// [[Rcpp::export]]
Rcpp::List RcppmeanvarN(arma::vec a,arma::vec b,arma::vec mu,arma::mat S, double p){
  Rcpp::List out;
  int n = mu.n_elem;
  arma::colvec s = sqrt(arma::diagvec(S));
  if(n==1){
    double a0 = arma::as_scalar(a);
    double b0 = arma::as_scalar(b);
    double mu0 = arma::as_scalar(mu);
    double S0 = arma::as_scalar(S);
    double s1 = sqrt(S0);
    double a1 = (a0-mu0)/s1;
    double b1 = (b0-mu0)/s1;
    double p = R::pnorm(b1,0,1,1,0) - R::pnorm(a1,0,1,1,0);
    double muY = mu0 + (R::dnorm(a1,0,1,0) - R::dnorm(b1,0,1,0))*s1/p;
    if(std::isinf(a0)){a0 = 0;}
    if(std::isinf(b0)){b0 = 0;}
    double varY = S0 +(mu0-muY)*muY+(a0*R::dnorm(a1,0,1,0)-b0*R::dnorm(b1,0,1,0))*s1/p;
    
    double EYY = varY + pow(muY,2);
    out["mean"] = muY;
    out["EYY"] = EYY;
    out["varcov"] = varY;
    return out;
  }else{
    arma::vec muY(n);
    arma::mat varY(n,n);
    arma::uvec seqq = arma::linspace<arma::uvec>(1,n,n);
    arma::vec a1 = a-mu;
    arma::vec b1 = b-mu;
    arma::vec ha = a1/s;
    arma::vec hb = b1/s;
    arma::mat R = S%(1/(s * s.t()));
    // NumericMatrix W = Rcpp::wrap(R);
    // double p = pmvtRcpp(ha,hb,W,0);
    //double p = 1;
    //double p = pmvt_cpp(ha,hb,R);
    Rcpp::List run = Rcppqfun(a1,b1,S);
    arma::vec qa = run["qa"];
    arma::vec qb = run["qb"];
    arma::vec q = qa-qb;
    muY = mu + S*q/p;
    arma::mat D(n,n,arma::fill::zeros);
    arma::uvec ii(1);
    double da, db;
    arma::vec mai,mbi;
    for(int i=1; i<=n; ++i){
      arma::vec wa(n-1,arma::fill::zeros);
      arma::vec wb(n-1,arma::fill::zeros);
      ii = i;
      arma::uvec ind2 = seqq.elem(find(seqq != i));
      arma::mat SSi = S(ind2-1,ind2-1)-S(ind2-1,ii-1)*S(ii-1,ind2-1)/arma::as_scalar(S(ii-1,ii-1));
      if(std::isfinite(a(i-1))){
        D(i-1,i-1) = a(i-1)*qa(i-1);
        mai = mu(ind2-1) + S(ind2-1,ii-1)/S(i-1,i-1)*a1(i-1);
        Rcpp::List run1 = Rcppqfun(a(ind2-1)-mai,b(ind2-1)-mai,SSi);
        arma::vec qa1 = run1["qa"];
        arma::vec qb1 = run1["qb"];
        da = R::dnorm(a(i-1),mu(i-1),s(i-1),0);
        wa = qa(i-1)*mai + da*SSi*(qa1-qb1);
      } 
      if(std::isfinite(b(i-1))){
        D(i-1,i-1) -= b(i-1)*qb(i-1);
        mbi = mu(ind2-1) + S(ind2-1,ii-1)/S(i-1,i-1)*b1(i-1);
        Rcpp::List run2 = Rcppqfun(a(ind2-1)-mbi,b(ind2-1)-mbi,SSi);
        arma::vec qa2 = run2["qa"];
        arma::vec qb2 = run2["qb"];
        db = R::dnorm(b(i-1),mu(i-1),s(i-1),0);
        wb = qb(i-1)*mbi + db*SSi*(qa2-qb2);
      }
      D(ii-1,ind2-1) = wa.t()-wb.t();
    }
    varY = S + S*(D - q*muY.t())/p;
    varY = (varY + varY.t())/2;
    
    arma::mat EYY(n,n);
    EYY = varY + muY * muY.t();
    out["mean"] = muY;
    out["EYY"] = EYY;
    out["varcov"] = varY;
    
    return out;
  }
}

//*************************************************************************
//MEANVAR -> TRUNCATED -> NORMAL
// [[Rcpp::export]]
Rcpp::List RcppmeanvarN_ab(arma::vec a,arma::vec b,arma::vec mu,arma::mat S, double p){
  Rcpp::List out;
  int n = mu.n_elem;
  arma::colvec s = sqrt(arma::diagvec(S));
  if(n==1){
    double a0 = arma::as_scalar(a);
    double b0 = arma::as_scalar(b);
    double mu0 = arma::as_scalar(mu);
    double S0 = arma::as_scalar(S);
    double s1 = sqrt(S0);
    double a1 = (a0-mu0)/s1;
    double b1 = (b0-mu0)/s1;
    double p = R::pnorm(b1,0,1,1,0) - R::pnorm(a1,0,1,1,0);
    double muY = mu0 + (R::dnorm(a1,0,1,0) - R::dnorm(b1,0,1,0))*s1/p;
    double varY = S0 +(mu0-muY)*muY+(a0*R::dnorm(a1,0,1,0)-b0*R::dnorm(b1,0,1,0))*s1/p;
    
    double EYY = varY + pow(muY,2);
    out["mean"] = muY;
    out["EYY"] = EYY;
    out["varcov"] = varY;
    return out;
  }else{
    arma::vec muY(n);
    arma::mat varY(n,n);
    arma::uvec seqq = arma::linspace<arma::uvec>(1,n,n);
    arma::vec a1 = a-mu;
    arma::vec b1 = b-mu;
    arma::vec ha = a1/s;
    arma::vec hb = b1/s;
    arma::mat R = S%(1/(s * s.t()));
    // NumericMatrix W = Rcpp::wrap(R);
    // double p = pmvtRcpp(ha,hb,W,0);
    //double p = 1;
    //double p = pmvt_cpp(ha,hb,R);
    Rcpp::List run = Rcppqfun(a1,b1,S);
    arma::vec qa = run["qa"];
    arma::vec qb = run["qb"];
    arma::vec q = qa-qb;
    muY = mu + S*q/p;
    arma::mat D(n,n,arma::fill::zeros);
    arma::uvec ii(1);
    double da, db;
    arma::vec mai,mbi;
    for(int i=1; i<=n; ++i){
      arma::vec wa(n-1,arma::fill::zeros);
      arma::vec wb(n-1,arma::fill::zeros);
      ii = i;
      arma::uvec ind2 = seqq.elem(find(seqq != i));
      arma::mat SSi = S(ind2-1,ind2-1)-S(ind2-1,ii-1)*S(ii-1,ind2-1)/arma::as_scalar(S(ii-1,ii-1));
      
        D(i-1,i-1) = a(i-1)*qa(i-1);
        mai = mu(ind2-1) + S(ind2-1,ii-1)/S(i-1,i-1)*a1(i-1);
        Rcpp::List run1 = Rcppqfun(a(ind2-1)-mai,b(ind2-1)-mai,SSi);
        arma::vec qa1 = run1["qa"];
        arma::vec qb1 = run1["qb"];
        da = R::dnorm(a(i-1),mu(i-1),s(i-1),0);
        wa = qa(i-1)*mai + da*SSi*(qa1-qb1);
      
        D(i-1,i-1) -= b(i-1)*qb(i-1);
        mbi = mu(ind2-1) + S(ind2-1,ii-1)/S(i-1,i-1)*b1(i-1);
        Rcpp::List run2 = Rcppqfun(a(ind2-1)-mbi,b(ind2-1)-mbi,SSi);
        arma::vec qa2 = run2["qa"];
        arma::vec qb2 = run2["qb"];
        db = R::dnorm(b(i-1),mu(i-1),s(i-1),0);
        wb = qb(i-1)*mbi + db*SSi*(qa2-qb2);
      
      D(ii-1,ind2-1) = wa.t()-wb.t();
    }
    varY = S + S*(D - q*muY.t())/p;
    varY = (varY + varY.t())/2;
    
    arma::mat EYY(n,n);
    EYY = varY + muY * muY.t();
    out["mean"] = muY;
    out["EYY"] = EYY;
    out["varcov"] = varY;
    
    return out;
  }
}

//*************************************************************************
//*************************************************************************
//ONLY MEAN -> TRUNCATED -> NORMAL
// [[Rcpp::export]]
Rcpp::List RcpponlymeanN(arma::vec a,arma::vec b,arma::vec mu,arma::mat S, double p){
  Rcpp::List out;
  int n = mu.n_elem;
  arma::colvec s = sqrt(arma::diagvec(S));
  if(n==1){
    double a0 = arma::as_scalar(a);
    double b0 = arma::as_scalar(b);
    double mu0 = arma::as_scalar(mu);
    double S0 = arma::as_scalar(S);
    double s1 = sqrt(S0);
    double a1 = (a0-mu0)/s1;
    double b1 = (b0-mu0)/s1;
    double p = R::pnorm(b1,0,1,1,0) - R::pnorm(a1,0,1,1,0);
    double muY = mu0 + (R::dnorm(a1,0,1,0) - R::dnorm(b1,0,1,0))*s1/p;
    out["mean"] = muY;
    return out;
  }else{
    arma::vec muY(n);
    arma::uvec seqq = arma::linspace<arma::uvec>(1,n,n);
    arma::vec a1 = a-mu;
    arma::vec b1 = b-mu;
    
    arma::vec ha = a1/s;
    arma::vec hb = b1/s;
    arma::mat R = S%(1/(s * s.t()));
    //double p = pmvt_cpp(ha,hb,R);
    //double p = 1;
    Rcpp::List run = Rcppqfun(a1,b1,S);
    arma::vec qa = run["qa"];
    arma::vec qb = run["qb"];
    arma::vec q = qa-qb;
    muY = mu + S*q/p;
    out["mean"] = muY;
    return out;
  }
}

//*************************************************************************
//MOMENTS -> TRUNCATED -> NORMAL
// [[Rcpp::export]]
Rcpp::List RcppmomentsN(arma::vec kappa,arma::vec a,arma::vec b, arma::vec mu, arma::mat S){
  Rcpp::List out = recintab1(kappa,a,b,mu,S);
  arma::vec M = out[1];
  //M = M/M(0);
  out[1] = M/M(0);
  return out;
}

