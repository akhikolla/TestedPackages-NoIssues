// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

//typedef Eigen::SparseMatrix<double> SpMat;
//typedef Eigen::SparseMatrix<double>::InnerIterator InIterMat;
//typedef Eigen::SparseVector<double> SpVec;
//typedef Eigen::SparseVector<double> SpVecf;
//typedef Eigen::SparseVector<int> SpVeci;
//typedef SpVec::InnerIterator InIterVec;
//typedef Eigen::SparseMatrix<double> SpMat;


//*****  Center and standardize  *****/
// [[Rcpp::export]]
List stdQ(Eigen::MatrixXd X){
  Eigen::VectorXd mX=X.colwise().mean(), sdX(X.cols()), sdXi;
  X.rowwise()-=mX.transpose();
  sdX=X.colwise().norm()/sqrt((double)X.rows());
  sdXi=1.0/sdX.array();
  X=X*sdXi.asDiagonal();
  return List::create(Named("x")=X, Named("sd")=sdX);
}

/*****  Log-pl of eta,  ties  *****/
// [[Rcpp::export]]
double plQ(Eigen::VectorXd& xb, Eigen::VectorXi& nevent,Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n){
  int i, j, q;
  double ll=0.0, SSi;
  if(loc1.size()==0){
    return 0.0;
  }else{
    Eigen::VectorXd exb = (xb.array()).exp();
    SSi=exb.sum();
    for(i=0;i<n;++i){
      for(j=loc1(i)-1, q=0;q<nevent1(i);j++, q++){ll+=xb(j);}
      ll-=nevent1(i)*log(SSi);
      for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){SSi-=exb(j);}
    }
    return(ll);
  }
}

/*****  1st order derivatives of log-pl of eta,  ties *****/
/*****  1st order derivatives of negavitve-log-pl of eta = -pl1*****/
// [[Rcpp::export]]
Eigen::VectorXd dfQ(Eigen::VectorXd& xb, Eigen::VectorXd& tevent, int& N, Eigen::VectorXi& nevent, 
                    Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n){
  int i, j, q;
  double denSi, c1=0.0;
  Eigen::VectorXd pl1(N);
  Eigen::VectorXd exb = (xb.array()).exp();
  denSi=exb.sum();
  for(i=0;i<n;++i){
    c1+=(nevent1(i)/denSi);
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){
      denSi-=exb(j);
      pl1(j)=tevent(j)-exb(j)*c1;
    }
  }
  return(pl1);
}

/*****  2nd order derivatives of log-pl of eta,  ties *****/
/*****  2nd order derivatives of negavitve-log-pl of eta = -pl2*****/
// [[Rcpp::export]]
Eigen::VectorXd d2Q(Eigen::VectorXd& xb, Eigen::VectorXd& tevent, int& N, Eigen::VectorXi& nevent, 
                    Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n){
  int i, j, q;
  double denSi, c1=0.0;
  Eigen::VectorXd pl2(N);
  Eigen::VectorXd exb = (xb.array()).exp();
  denSi=exb.sum();
  for(i=0;i<n;++i){
    c1+=(nevent1(i)/denSi);
    //c2+=(nevent1(i)/pow(denSi, 2)); this part is relatively small comparing to the first part, so can be ignored
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){
      denSi-=exb(j);
      pl2(j)=exb(j)*c1;
    }
  }
  return(pl2);
}

/***** Lambda path (max)  *****/
// [[Rcpp::export]]
double max_lambda(List& dt){
  //information for Cox
  Eigen::MatrixXd X;
  Eigen::VectorXd tevent, m, l_new;
  Eigen::VectorXi nevent, nevent1, loc1, g, K1;
  int p, N, n;
  // Extract information
  X = dt[0]; tevent = dt[2]; N = dt[3]; nevent = dt[4]; 
  nevent1 = dt[5]; loc1 = dt[6]; n = dt[7]; g = dt[8]; 
  m = dt[9]; K1 = dt[11]; p =dt[13];
  
  Eigen::VectorXd r = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd rb = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd xb = Eigen::VectorXd::Zero(N);
  //number of groups
  int ig = K1.size()-1;
  Eigen::VectorXd lmd = Eigen::VectorXd::Zero(ig);
  
  r = Eigen::VectorXd::Zero(N);
  rb = Eigen::VectorXd::Zero(p);
  //calculate first-order derivative r
  r = dfQ(xb, tevent, N, nevent, nevent1, loc1, n);
  for(int i=0; i<p; i++){rb(i) = (r.transpose()).dot(X.col(i));}
  
  for(int i=0; i<ig; i++){
    int kg = K1(i+1)-K1(i);
    l_new = Eigen::VectorXd::Zero(kg);
    for(int j=K1(i); j<K1(i+1); j++){
      l_new(j-K1(i)) = rb(j);
    }
    lmd(i) = l_new.norm()/(N*sqrt((double)(kg)));
  }
  return (lmd.maxCoeff());
}

//*****  convert Eigen::VectorXd to NumericVector  *****/
// [[Rcpp::export]]
NumericVector convd2n(Eigen::VectorXd X) {
  Rcpp::NumericVector a(wrap(X));
  return a; 
}

//*****  check non-numeric values  *****/
// [[Rcpp::export]]
int isNA(NumericVector x) {
  int n = x.size(), res = 0;
  LogicalVector out(n);
  
  for (int i = 0; i < n; ++i) {
    out[i] = NumericVector::is_na(x[i]);
    if(out[i]){
      res = 1;
      break;
    }
  }
  return res;
}

//*****  check Infinite values  *****/
// [[Rcpp::export]]
int isInf(NumericVector x) {
  int n = x.size(), res = 0;
  LogicalVector out(n);
  
  for (int i = 0; i < n; ++i) {
    out[i] = Rcpp::traits::is_infinite<REALSXP>(x[i]);
    if(out[i]){
      res = 1;
      break;
    }
  }
  return res;
}

/*****  Compare two strings *****/
// [[Rcpp::export]]
int compare(const char* cx, const char* cy){
  // 0 = the same; #0 = different
  return std::strcmp(cx,cy);
}

//*****  group-descent algorithm  *****/
// [[Rcpp::export]]
List grpCoxQ(List& dt, const char* penalty, Eigen::VectorXd lambda, int nlambda, 
             double gamma, double thresh, int maxit){
  //information for Cox
  Eigen::MatrixXd X;
  Eigen::VectorXd tevent, m;
  Eigen::VectorXi nevent, nevent1, loc1, g, K1;
  int p, N, n;
  // Extract information
  X = dt[0]; tevent = dt[2]; N = dt[3]; nevent = dt[4]; 
  nevent1 = dt[5]; loc1 = dt[6]; n = dt[7]; g = dt[8]; 
  m = dt[9]; K1 = dt[11]; p =dt[13];
  
  //intermediate quantities
  //log-likelihood values
  double llQ;
  //beta from previous iteration and current beta
  Eigen::VectorXd beta0 = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd betac = Eigen::VectorXd::Zero(p);
  //xb=eta and first-order derivative value
  Eigen::VectorXd xb = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd r = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd rb = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd s2 = Eigen::VectorXd::Zero(N);
  //number of groups
  int ig = K1.size()-1;
  //active vector
  Eigen::VectorXi e = Eigen::VectorXi::Zero(ig);
  Eigen::VectorXd l_new;
  //others
  int il, it=0, kg, iadd;
  double tau, maxChange, lambdaj, norm_lnew, shift;
  
  //outcome
  Eigen::MatrixXd Betas = Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd locbeta = Eigen::VectorXd::Zero(nlambda);
  //number of iterations of each beta value
  Eigen::VectorXi iter = Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXi flag = Eigen::VectorXi::Zero(nlambda);
  
  //Initialization
  // locbeta(0)
  if(loc1.size()==0){
    llQ = 0.0;
  }else{
    xb = Eigen::VectorXd::Zero(N);
    llQ = plQ(xb,nevent,nevent1,loc1,n);
  }
  locbeta(0) = llQ;
  
  //lambda path
  for(il=1; il<nlambda; il++){
    beta0 = Betas.col(il-1);
    if(loc1.size()==0){
      betac = Eigen::VectorXd::Zero(p);
      llQ = 0.0;
    }else{
      s2= Eigen::VectorXd::Zero(N);
      s2 = d2Q(xb, tevent, N, nevent, nevent1, loc1, n)/N;
      //checking
      if(N>p){
        tau = (N/p)*N*s2.maxCoeff();
      }else{
        if(compare(penalty, "glasso")==0){
          tau = p*s2.maxCoeff(); //(p/N)*N
        }else{
          tau = N*s2.maxCoeff();
        }
      }
      
      //Start
      betac = Eigen::VectorXd::Zero(p);
      
      r=Eigen::VectorXd::Zero(N);
      e=Eigen::VectorXi::Zero(ig);
      //scan for active
      for(int i=0; i<ig; i++){
        if(e(i)==0){
          for(int j=K1(i); j<K1(i+1); j++){
            if(beta0(j) != 0){e(i) = 1; break;}
          }
        }
      }
      
      it = 0;
      //while loop
      while(it<maxit){
        while(it<maxit){
          iter(il)++;
          it++;
          
          xb = Eigen::VectorXd::Zero(N);
          for(int i=0; i<p; i++){xb += (X.col(i))*beta0(i);}
          r = Eigen::VectorXd::Zero(N);
          rb = Eigen::VectorXd::Zero(p);
          //calculate first-order derivative r
          r = dfQ(xb, tevent, N, nevent, nevent1, loc1, n)/N;
          //check for Inf: // warning("Overflow problem...");
          if((isNA(convd2n(r))==1) || (isInf(convd2n(r))==1)){flag(il)=1;break;} 
          
          for(int i=0; i<p; i++){rb(i) = (r.transpose()).dot(X.col(i));}
          //check for Inf: // warning("Overflow problem...");
          if((isNA(convd2n(rb))==1) || (isInf(convd2n(rb))==1)){flag(il)=1;break;} 
          
          //update beta
          maxChange=0;
          for(int i=0; i<ig; i++){
            if(e(i)){
              kg = K1(i+1)-K1(i);
              lambdaj = lambda(il)*sqrt((double)(kg));
              l_new = Eigen::VectorXd::Zero(kg);
              for(int j=K1(i); j<K1(i+1); j++){
                l_new(j-K1(i)) = rb(j)/tau + beta0(j);
              }
              
              norm_lnew = l_new.norm();
              //update beta
              // group LASSO
              if(compare(penalty, "glasso")==0){
                if(norm_lnew > lambdaj/tau){
                  for(int j=K1(i); j<K1(i+1); j++){
                    betac(j) = (1-lambdaj/(tau*norm_lnew)) * l_new(j-K1(i));
                    shift = betac(j) - beta0(j);
                    if(std::abs(shift) > maxChange) maxChange=std::abs(shift);
                    r -= shift*(X.col(j))/N;
                    xb += shift*(X.col(j));
                  }
                }
              }
              // group SCAD
              if(compare(penalty, "gSCAD")==0){
                if(norm_lnew <= (lambdaj*(1+1/tau))){
                  if(norm_lnew > (lambdaj/tau)){
                    for(int j=K1(i); j<K1(i+1); j++){
                      betac(j) = (1-lambdaj/(tau*norm_lnew)) * l_new(j-K1(i));
                      shift = betac(j) - beta0(j);
                      if(std::abs(shift) > maxChange) maxChange=std::abs(shift);
                      r -= shift*(X.col(j))/N;
                      xb += shift*(X.col(j));
                    }
                  }
                }
                if((norm_lnew > (lambdaj*(1+1/tau))) & (norm_lnew <= (lambdaj*gamma))){
                  for(int j=K1(i); j<K1(i+1); j++){
                    betac(j) = tau*(gamma-1)/(tau*(gamma-1)-1)*(1-lambdaj*gamma/(tau*(gamma-1)*norm_lnew)) * l_new(j-K1(i));
                    shift = betac(j) - beta0(j);
                    if(std::abs(shift) > maxChange) maxChange=std::abs(shift);
                    r -= shift*(X.col(j))/N;
                    xb += shift*(X.col(j));
                  }
                }
                if(norm_lnew > (lambdaj*gamma)){
                  for(int j=K1(i); j<K1(i+1); j++){
                    betac(j) = l_new(j-K1(i));
                    shift = betac(j) - beta0(j);
                    if(std::abs(shift) > maxChange) maxChange=std::abs(shift);
                    r -= shift*(X.col(j))/N;
                    xb += shift*(X.col(j));
                  }
                }
              }
              // group MCP
              if(compare(penalty, "gMCP")==0){
                if(norm_lnew <= (lambdaj*gamma)){
                  if(norm_lnew > (lambdaj/tau)){
                    for(int j=K1(i); j<K1(i+1); j++){
                      betac(j) = tau*gamma/(tau*gamma-1)*(1-lambdaj/(tau*norm_lnew)) * l_new(j-K1(i));
                      shift = betac(j) - beta0(j);
                      if(std::abs(shift) > maxChange) maxChange=std::abs(shift);
                      r -= shift*(X.col(j))/N;
                      xb += shift*(X.col(j));
                    }
                  }
                }
                if(norm_lnew > (lambdaj*gamma)){
                  for(int j=K1(i); j<K1(i+1); j++){
                    betac(j) = l_new(j-K1(i));
                    shift = betac(j) - beta0(j);
                    if(std::abs(shift) > maxChange) maxChange=std::abs(shift);
                    r -= shift*(X.col(j))/N;
                    xb += shift*(X.col(j));
                  }
                }
              }
            }//end-e
          }//end update beta
          
          //check convergence
          for(int j=0; j<p; j++) beta0(j)=betac(j);
          if(maxChange<thresh){break;}
        }//end-while-loop
        
        //scan for active 
        iadd = 0;
        for(int i=0; i<ig; i++){
          if(e(i)==0){
            kg = K1(i+1)-K1(i);
            lambdaj = lambda(il)*sqrt((double)(kg));
            l_new = Eigen::VectorXd::Zero(kg);
            for(int j=K1(i); j<K1(i+1); j++){
              l_new(j-K1(i)) = rb(j)/tau + beta0(j);
            }
            
            norm_lnew = l_new.norm();
            //update beta
            // group LASSO
            if(compare(penalty, "glasso")==0){
              if(norm_lnew > lambdaj/tau){
                for(int j=K1(i); j<K1(i+1); j++){
                  betac(j) = (1-lambdaj/(tau*norm_lnew)) * l_new(j-K1(i));
                }
              }
            }
            // group SCAD
            if(compare(penalty, "gSCAD")==0){
              if(norm_lnew <= (lambdaj*(1+1/tau))){
                if(norm_lnew > (lambdaj/tau)){
                  for(int j=K1(i); j<K1(i+1); j++){
                    betac(j) = (1-lambdaj/(tau*norm_lnew)) * l_new(j-K1(i));
                  }
                }
              }
              if((norm_lnew > (lambdaj*(1+1/tau))) & (norm_lnew <= (lambdaj*gamma))){
                for(int j=K1(i); j<K1(i+1); j++){
                  betac(j) = tau*(gamma-1)/(tau*(gamma-1)-1)*(1-lambdaj*gamma/(tau*(gamma-1)*norm_lnew)) * l_new(j-K1(i));
                }
              }
              if(norm_lnew > (lambdaj*gamma)){
                for(int j=K1(i); j<K1(i+1); j++){
                  betac(j) = l_new(j-K1(i));
                }
              }
            }
            // group MCP
            if(compare(penalty, "gMCP")==0){
              if(norm_lnew <= (lambdaj*gamma)){
                if(norm_lnew > (lambdaj/tau)){
                  for(int j=K1(i); j<K1(i+1); j++){
                    betac(j) = tau*gamma/(tau*gamma-1)*(1-lambdaj/(tau*norm_lnew)) * l_new(j-K1(i));
                  }
                }
              }
              if(norm_lnew > (lambdaj*gamma)){
                for(int j=K1(i); j<K1(i+1); j++){
                  betac(j) = l_new(j-K1(i));
                }
              }
            }
            //scan active again
            for(int j=K1(i); j<K1(i+1); j++){
              if(betac(j) != 0){
                e(i) = 1;
                iadd = 1;
                break;
              }
            }
          }//end-e
        }
        if(iadd==0){break;}
        
        for(int j=0; j<p; j++){beta0(j) = betac(j);}  
      }//end-while-loop
      
      xb = Eigen::VectorXd::Zero(N);
      for(int i=0; i<p; i++){xb += (X.col(i))*betac(i);}
      llQ = plQ(xb,nevent,nevent1,loc1,n);
    }
    
    if((iter(il) == maxit) || (llQ == locbeta(0))){flag(il)=2;}
    Betas.col(il) = betac;
    locbeta(il) = llQ;
  }//end-lambda-loop
  
  return(List::create(Named("beta")=Betas, Named("flag")=flag, Named("ll")=locbeta, Named("nlambda")=nlambda, Named("iter")=iter));
}

//*****  cross-validation  *****/
// [[Rcpp::export]]
List cvgrpCoxQ(List& dt, const char* penalty, Eigen::VectorXd lambda, int nlambda, double gamma,
               double thresh, int maxit, List& dF){
  //information for Cox
  Eigen::MatrixXd X;
  Eigen::VectorXd tevent, m;
  Eigen::VectorXi nevent, nevent1, loc1, g, K1;
  int p, N, n;
  // Extract information
  X = dt[0]; tevent = dt[2]; N = dt[3]; nevent = dt[4]; 
  nevent1 = dt[5]; loc1 = dt[6]; n = dt[7]; g = dt[8]; 
  m = dt[9]; K1 = dt[11]; p =dt[13];
  
  //information dF
  Eigen::MatrixXd XF;
  Eigen::VectorXd teventF, mF;
  Eigen::VectorXi neventF, nevent1F, loc1F;
  int NF, nF;
  // Extract information
  XF = dF[0]; teventF = dF[2]; NF = dF[3]; neventF = dF[4]; 
  nevent1F = dF[5]; loc1F = dF[6]; nF = dF[7];
  
  //intermediate quantities
  //log-likelihood values
  double llQ, llF;
  //beta from previous iteration and current beta
  Eigen::VectorXd beta0 = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd betac = Eigen::VectorXd::Zero(p);
  //xb=eta and first-order derivative value
  Eigen::VectorXd xb = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd xbF = Eigen::VectorXd::Zero(NF);
  Eigen::VectorXd r = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd rb = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd s2 = Eigen::VectorXd::Zero(N);
  //number of groups
  int ig = K1.size()-1;
  //active vector
  Eigen::VectorXi e = Eigen::VectorXi::Zero(ig);
  Eigen::VectorXd l_new;
  //others
  int il, it=0, kg, iadd;
  double tau, maxChange, lambdaj, norm_lnew, shift;
  
  //outcome
  Eigen::MatrixXd Betas = Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd locbeta = Eigen::VectorXd::Zero(nlambda);
  Eigen::VectorXd locbetaF = Eigen::VectorXd::Zero(nlambda);
  //number of iterations of each beta value
  Eigen::VectorXi iter = Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXi flag = Eigen::VectorXi::Zero(nlambda);
  
  //Initialization
  // locbeta(0)
  if(loc1.size()==0){
    llQ = 0.0;
  }else{
    xb = Eigen::VectorXd::Zero(N);
    llQ = plQ(xb,nevent,nevent1,loc1,n);
  }
  locbeta(0) = llQ;
  if(loc1F.size()==0){
    llF = 0.0;
  }else{
    xbF = Eigen::VectorXd::Zero(NF);
    llF = plQ(xbF,neventF,nevent1F,loc1F,nF);
  }
  locbetaF(0) = llF;
  
  //lambda path
  for(il=1; il<nlambda; il++){
    //assign beta0Q 
    beta0 = Betas.col(il-1);
    if(loc1.size()==0){
      //results of q-transition
      for(int i=0; i<p; i++){betac(i) = 0.0;}
      llQ = 0.0;
    }else{
      s2= Eigen::VectorXd::Zero(N);
      s2 = d2Q(xb, tevent, N, nevent, nevent1, loc1, n)/N;
      //checking
      if(N>p){
        tau = (N/p)*N*s2.maxCoeff();
      }else{
        if(compare(penalty, "glasso")==0){
          tau = p*s2.maxCoeff(); //(p/N)*N
        }else{
          tau = N*s2.maxCoeff();
        }
      }
      
      //Start
      betac = Eigen::VectorXd::Zero(p);
      
      r=Eigen::VectorXd::Zero(N);
      e=Eigen::VectorXi::Zero(ig);
      //scan for active
      for(int i=0; i<ig; i++){
        if(e(i)==0){
          for(int j=K1(i); j<K1(i+1); j++){
            if(beta0(j) != 0){e(i) = 1; break;}
          }
        }
      }
      
      it = 0;
      //while loop
      while(it<maxit){
        while(it<maxit){
          iter(il)++;
          it++;
          
          xb = Eigen::VectorXd::Zero(N);
          for(int i=0; i<p; i++){xb += (X.col(i))*beta0(i);}
          r = Eigen::VectorXd::Zero(N);
          rb = Eigen::VectorXd::Zero(p);
          //calculate first-order derivative r
          r = dfQ(xb, tevent, N, nevent, nevent1, loc1, n)/N;
          //check for Inf: // warning("Overflow problem...");
          if((isNA(convd2n(r))==1) || (isInf(convd2n(r))==1)){flag(il)=1;break;} 
          
          for(int i=0; i<p; i++){rb(i) = (r.transpose()).dot(X.col(i));}
          //check for Inf: // warning("Overflow problem...");
          if((isNA(convd2n(rb))==1) || (isInf(convd2n(rb))==1)){flag(il)=1;break;} 
          
          //update beta
          maxChange=0;
          for(int i=0; i<ig; i++){
            if(e(i)){
              kg = K1(i+1)-K1(i);
              lambdaj = lambda(il)*sqrt((double)(kg));
              l_new = Eigen::VectorXd::Zero(kg);
              for(int j=K1(i); j<K1(i+1); j++){
                l_new(j-K1(i)) = rb(j)/tau + beta0(j);
              }
              
              norm_lnew = l_new.norm();
              //update beta
              // group LASSO
              if(compare(penalty, "glasso")==0){
                if(norm_lnew > lambdaj/tau){
                  for(int j=K1(i); j<K1(i+1); j++){
                    betac(j) = (1-lambdaj/(tau*norm_lnew)) * l_new(j-K1(i));
                    shift = betac(j) - beta0(j);
                    if(std::abs(shift) > maxChange) maxChange=std::abs(shift);
                    r -= shift*(X.col(j))/N;
                    xb += shift*(X.col(j));
                  }
                }
              }
              // group SCAD
              if(compare(penalty, "gSCAD")==0){
                if((norm_lnew <= (lambdaj*(1+1/tau))) && (tau*(gamma-1)>1)){
                  if(norm_lnew > (lambdaj/tau)){
                    for(int j=K1(i); j<K1(i+1); j++){
                      betac(j) = (1-lambdaj/(tau*norm_lnew)) * l_new(j-K1(i));
                      shift = betac(j) - beta0(j);
                      if(std::abs(shift) > maxChange) maxChange=std::abs(shift);
                      r -= shift*(X.col(j))/N;
                      xb += shift*(X.col(j));
                    }
                  }
                }
                if((norm_lnew > (lambdaj*(1+1/tau))) && (norm_lnew <= (lambdaj*gamma))){
                  for(int j=K1(i); j<K1(i+1); j++){
                    betac(j) = tau*(gamma-1)/(tau*(gamma-1)-1)*(1-lambdaj*gamma/(tau*(gamma-1)*norm_lnew)) * l_new(j-K1(i));
                    shift = betac(j) - beta0(j);
                    if(std::abs(shift) > maxChange) maxChange=std::abs(shift);
                    r -= shift*(X.col(j))/N;
                    xb += shift*(X.col(j));
                  }
                }
                if(norm_lnew > (lambdaj*gamma)){
                  for(int j=K1(i); j<K1(i+1); j++){
                    betac(j) = l_new(j-K1(i));
                    shift = betac(j) - beta0(j);
                    if(std::abs(shift) > maxChange) maxChange=std::abs(shift);
                    r -= shift*(X.col(j))/N;
                    xb += shift*(X.col(j));
                  }
                }
              }
              // group MCP
              if(compare(penalty, "gMCP")==0){
                if((norm_lnew <= (lambdaj*gamma)) & (tau*gamma>1)){
                  if(norm_lnew > (lambdaj/tau)){
                    for(int j=K1(i); j<K1(i+1); j++){
                      betac(j) = tau*gamma/(tau*gamma-1)*(1-lambdaj/(tau*norm_lnew)) * l_new(j-K1(i));
                      shift = betac(j) - beta0(j);
                      if(std::abs(shift) > maxChange) maxChange=std::abs(shift);
                      r -= shift*(X.col(j))/N;
                      xb += shift*(X.col(j));
                    }
                  }
                }
                if(norm_lnew > (lambdaj*gamma)){
                  for(int j=K1(i); j<K1(i+1); j++){
                    betac(j) = l_new(j-K1(i));
                    shift = betac(j) - beta0(j);
                    if(std::abs(shift) > maxChange) maxChange=std::abs(shift);
                    r -= shift*(X.col(j))/N;
                    xb += shift*(X.col(j));
                  }
                }
              }
            }//end-e
          }//end update beta
          
          //check convergence
          for(int j=0; j<p; j++) beta0(j)=betac(j);
          if(maxChange<thresh){break;}
        }//end-while-loop
        
        //scan for active 
        iadd = 0;
        for(int i=0; i<ig; i++){
          if(e(i)==0){
            kg = K1(i+1)-K1(i);
            lambdaj = lambda(il)*sqrt((double)(kg));
            l_new = Eigen::VectorXd::Zero(kg);
            for(int j=K1(i); j<K1(i+1); j++){
              l_new(j-K1(i)) = rb(j)/tau + beta0(j);
            }
            
            norm_lnew = l_new.norm();
            //update beta
            // group LASSO
            if(compare(penalty, "glasso")==0){
              if(norm_lnew > lambdaj/tau){
                for(int j=K1(i); j<K1(i+1); j++){
                  betac(j) = (1-lambdaj/(tau*norm_lnew)) * l_new(j-K1(i));
                }
              }
            }
            // group SCAD
            if(compare(penalty, "gSCAD")==0){
              if((norm_lnew <= (lambdaj*(1+1/tau))) && (tau*(gamma-1)>1)){
                if(norm_lnew > (lambdaj/tau)){
                  for(int j=K1(i); j<K1(i+1); j++){
                    betac(j) = (1-lambdaj/(tau*norm_lnew)) * l_new(j-K1(i));
                  }
                }
              }
              if((norm_lnew > (lambdaj*(1+1/tau))) && (norm_lnew <= (lambdaj*gamma))){
                for(int j=K1(i); j<K1(i+1); j++){
                  betac(j) = tau*(gamma-1)/(tau*(gamma-1)-1)*(1-lambdaj*gamma/(tau*(gamma-1)*norm_lnew)) * l_new(j-K1(i));
                }
              }
              if(norm_lnew > (lambdaj*gamma)){
                for(int j=K1(i); j<K1(i+1); j++){
                  betac(j) = l_new(j-K1(i));
                }
              }
            }
            // group MCP
            if(compare(penalty, "gMCP")==0){
              if((norm_lnew <= (lambdaj*gamma)) & (tau*gamma>1)){
                if(norm_lnew > (lambdaj/tau)){
                  for(int j=K1(i); j<K1(i+1); j++){
                    betac(j) = tau*gamma/(tau*gamma-1)*(1-lambdaj/(tau*norm_lnew)) * l_new(j-K1(i));
                  }
                }
              }
              if(norm_lnew > (lambdaj*gamma)){
                for(int j=K1(i); j<K1(i+1); j++){
                  betac(j) = l_new(j-K1(i));
                }
              }
            }
            //scan active again
            for(int j=K1(i); j<K1(i+1); j++){
              if(betac(j) != 0){
                e(i) = 1;
                iadd = 1;
                break;
              }
            }
          }//end-e
        }
        if(iadd==0){break;}
        
        for(int j=0; j<p; j++){beta0(j) = betac(j);}  
      }//end-while-loop
      
      xb = Eigen::VectorXd::Zero(N);
      for(int i=0; i<p; i++){xb += (X.col(i))*betac(i);}
      llQ = plQ(xb,nevent,nevent1,loc1,n);
    }
    
    if((iter(il) == maxit) || (llQ == locbeta(0))){flag(il)=2;}
    locbeta(il) = llQ;
    Betas.col(il) = betac;
    // compute the cross-validation log-likelihood of the whole dataset
    xbF = Eigen::VectorXd::Zero(NF);
    for(int i=0; i<p; i++){
      xbF += XF.col(i)*betac(i);
    }
    // log-likelihood value
    llF = plQ(xbF,neventF,nevent1F,loc1F,nF);
    locbetaF(il) = llF;
    
  }//end-lambda-loop
  
  return(List::create(Named("beta")=Betas, Named("flag")=flag, Named("ll")=locbeta, Named("lf")=locbetaF, Named("nlambda")=nlambda, Named("iter")=iter));
}
