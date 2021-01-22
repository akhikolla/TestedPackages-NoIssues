
#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>
using namespace Rcpp;


double min(double x, double y){
  if(x<y){
    return(x);
  }else{
    return(y);
  }
}

double uber(double x, double y){
  return(tgamma(x+1)/(tgamma(y+1)*tgamma(x-y+1)));
}

double beta(double p, double q){
  return(tgamma(p)*tgamma(q)/tgamma(p+q));
}

double dnbinom(int k, double mu, double eta){
  return(::Rf_dnbinom_mu(k, eta, mu, 0));
}

double dbinom(int k, int n, double p){
  return(::Rf_dbinom(k, n, p, 0));
}

double dnbinomCond(int j, int k, double mu, double size, double rho){
  int y;
  double sumh;

  sumh=0;
  for(y=0;y<(min(j,k)+1);y++){
    sumh += uber(j, y)*beta(rho*size+y, (1-rho)*size+j-y)/beta(rho*size, (1-rho)*size)*pow(size/(size+mu), (1-rho)*size)*tgamma(k+(1-rho)*size-y)/tgamma(k-y+1)/tgamma((1-rho)*size)*pow(mu/(size+mu), k-y);
  }
  return sumh;
}

// [[Rcpp::export]]
double minFunc(NumericVector x, NumericVector daten, int dataNA){

  int t;
  double logL, mu, size, rho;

  logL=0;
  mu = x[0];
  size=x[1];
  rho=x[2];

  logL += log(dnbinom(daten[0], mu, size));
  for(t=0;t<(dataNA-1);t++){
    logL += log(dnbinomCond(daten[t], daten[t+1], mu, size, rho));
  }

  return -logL;
}

// [[Rcpp::export]]
double minFuncMult(NumericVector x, NumericMatrix daten, NumericVector dataNA, int n){

  int j,t;
  double logL, mu, size, rho;

  logL=0;
  mu = x[0];
  size=x[1];
  rho=x[2];

  for(j=0;j<n;j++){
    logL += log(dnbinom(daten(j,0), mu, size));
    for(t=0;t<(dataNA[j]-1);t++){
      logL += log(dnbinomCond(daten(j,t), daten(j,t+1), mu, size, rho));
    }

  }

  return -logL;
}

// [[Rcpp::export]]
double minFuncBlinded(NumericVector x, NumericMatrix daten, NumericVector dataNA, NumericVector n, double delta){
  int j, t, nAll;
  double logL, muC, muE, size, rho, k;

  logL=0;
  nAll = n[0]+n[1];
  k=n[1]/n[0];
  muE=x[0]*(1+k)/(k+1/delta);
  muC=x[0]*(1+k)/(1+k*delta);
  size=x[1];
  rho=x[2];

  for(j=0;j<nAll;j++){
    logL += log(1/(1+k)*(k*dnbinom(daten(j,0), muE, size)+dnbinom(daten(j,0), muC, size)));
    for(t=0;t<(dataNA[j]-1);t++){
      logL += log(1/(1+k)*(k*dnbinomCond(daten(j,t), daten(j,t+1), muE, size, rho) + dnbinomCond(daten(j,t), daten(j,t+1), muC, size, rho)));
    }
  }
  return -logL;
}



double ldnbinom(int x, double mu, double eta){
  if(x < 100){
    return(x*(log(mu)-log(mu+eta))+log(tgamma(x+eta))-log(tgamma(x+1))-log(tgamma(eta))+eta*(log(eta)-log(mu+eta)));
  }else{
    return(x*(log(mu)-log(mu+eta))+(eta-1)*log((double) x)-log(tgamma(eta))+eta*(log(eta)-log(mu+eta)));
  }
}
double ldnbinomDmu(int x, double mu, double eta){
  return(x*(1/mu-1/(mu+eta))-eta/(mu+eta));
}
double ldnbinomDmuDmu(int x, double mu, double eta){
  return(x*(1/pow(mu+eta,2)-1/pow(mu,2))+eta/pow(mu+eta,2));
}
double digamma(double k){
  return(Rf_digamma(k));
}
double ldnbinomDeta(int x, double mu, double eta){
  return(-(x+eta)/(mu+eta)+log(eta/(mu+eta))+digamma(x+eta)-digamma(eta)+1);
}
double trigamma(double k){
  return(Rf_trigamma(k));
}
double ldnbinomDetaDeta(int x, double mu, double eta){
  return((x-mu)/pow(mu+eta,2)+1/eta-1/(mu+eta)+trigamma(x+eta)-trigamma(eta));
}
double ldnbinomDmuDeta(int x, double mu, double eta){
  return((x-mu)/pow(mu+eta,2));
}

/*Terms required for expected Hessian and outer gradient product matrix*/


double dnbinomPair(int kt, int ks, double mut, double mus, double eta, double rho){
  //ACHTUNG!! rho=cor(kt,ks)
  int s, t;
  double term1, term2, term3, term4, erg = 0;

  if(rho==0){
    erg = dnbinom(kt, mut, eta)*dnbinom(ks, mus, eta);
  }else if(rho==1){
    erg = dnbinom(kt+ks, mus+mut, eta)*dbinom(ks, ks+kt, mus/(mus+mut));
  }else{
    for(s=0;s<=ks;s++){
      for(t=0;t<=kt;t++){
        term1 = dnbinom(s, mus*(1-rho), eta*(1-rho));
        term2 = dnbinom(t, mut*(1-rho), eta*(1-rho));
        term3 = dnbinom(ks+kt-s-t, (mus+mut)*rho, eta*rho);
        term4 = dbinom(ks-s, ks+kt-s-t, mus/(mus+mut));
        erg+=term1*term2*term3*term4;
      }
    }
  }
  return(erg);
}
double ldnbinomDmuDmuExp(double mu, double eta){
  return(1/(mu+eta)-1/mu);
}
double ExpTerm1(double mu, double eta, int approx = 50){
  int i;
  double erg=0;
  for(i=0;i<approx;i++){
    erg+=digamma(i+eta)*dnbinom(i, mu, eta);
  }
  return(erg);
}
double ExpTerm2(double mut, double mus, double eta, double rho, int approx = 50){
  //ACHTUNG!! rho=cor(kt,ks)
  int i,j;
  double erg=0;
  for(i=0;i<approx;i++){
    for(j=0;j<approx;j++){
      erg+=j*digamma(i+eta)*dnbinomPair(j,i, mut, mus, eta, rho);
    }
  }
  return(erg);
}
double ExpTerm3(double mut, double mus, double eta, double rho, int approx = 50){
  //ACHTUNG!! rho=cor(kt,ks)
  int i,j;
  double erg=0;
  for(i=0;i<approx;i++){
    for(j=0;j<approx;j++){
      erg+=digamma(j+eta)*digamma(i+eta)*dnbinomPair(j,i, mut, mus, eta, rho);
    }
  }
  return(erg);
}



/* Three trends available. New trends require extensions in this part*/
double trend(NumericVector lambda, int i, int t, int type){
  //i=1 corresponds to treatment group, i=2 to control group
  //Type 1 constant trend lambda=(lambda_1, lambda_2)
  if(type == 1){
    if(i==1){
      return(exp(lambda(0)+lambda(1)));
    }else{
      return(exp(lambda(0)));
    }
    //Type 2 exponential trend lambda=c(lambda_0,lambda_1,lambda_2)
  }else if(type==2){
    if(i==1){
      return(exp(lambda(0)+t*lambda(1)+t*lambda(2)));
    }else{
      return(exp(lambda(0)+t*lambda(1)));
    }
  }else{
    return(0);
  }
}
NumericVector trendGrad(NumericVector lambda, int i, int t, int type){
  //i=1 corresponds to treatment group, i=2 to control group
  //Type 1 constant trend lambda=(lambda_1, lambda_2)
  if(type == 1){
    NumericVector grad(2);
    if(i==1){
      grad(0) = exp(lambda(0)+lambda(1));
      grad(1) = exp(lambda(0)+lambda(1));
      return(grad);
    }else{
      grad(0) = exp(lambda(0));
      grad(1) = 0;
      return(grad);
    }
    //Type 2 exponential trend lambda=c(lambda_1, lambda_2, lambda_3)
  }else if(type==2){
    NumericVector grad(3);
    if(i==1){
      grad(0) = exp(lambda(0)+t*(lambda(1)+lambda(2)));
      grad(1) = exp(lambda(0)+t*(lambda(1)+lambda(2)))*t;
      grad(2) = exp(lambda(0)+t*(lambda(1)+lambda(2)))*t;
      return(grad);
    }else{
      grad(0) = exp(lambda(0)+t*lambda(1));
      grad(1) = exp(lambda(0)+t*lambda(1))*t;
      grad(2) = 0;
      return(grad);
    }
  }else{
    return(0);
  }
}
NumericMatrix trendHess(NumericVector lambda, int i, int t, int type){
  //i=1 corresponds to treatment group, i=2 to control group

  //Type 1 constant trend lambda=(lambda_1, lambda_2)
  if(type == 1){
    NumericMatrix hess(2,2);
    if(i==1){
      hess(0,0)=exp(lambda(0)+lambda(1));
      hess(0,1)=exp(lambda(0)+lambda(1));
      hess(1,0)=exp(lambda(0)+lambda(1));
      hess(1,1)=exp(lambda(0)+lambda(1));
      return(hess);
    }else{
      hess(0,0)=exp(lambda(0));
      hess(0,1)=0;
      hess(1,0)=0;
      hess(1,1)=0;
      return(hess);
    }
    //Type 2 exponential trend lambda=c(lambda_1,lambda_2,lambda_3)
  }else if(type==2){
    NumericMatrix hess(3,3);
    if(i==1){
      hess(0,0)=exp(lambda(0)+t*(lambda(1)+lambda(2)));
      hess(0,1)=t*exp(lambda(0)+t*lambda(1)+t*lambda(2));hess(1,0)=hess(0,1);
      hess(0,2)=t*exp(lambda(0)+t*lambda(1)+t*lambda(2));hess(2,0)=hess(0,2);
      hess(1,1)=pow((double) t,(double) 2)*exp(lambda(0)+t*lambda(1)+t*lambda(2));
      hess(1,2)=pow((double) t,(double) 2)*exp(lambda(0)+t*lambda(1)+t*lambda(2));hess(2,1)=hess(1,2);
      hess(2,2)=pow((double) t,(double) 2)*exp(lambda(0)+t*lambda(1)+t*lambda(2));
      return(hess);
    }else{
      hess(0,0)=exp(lambda(0)+t*lambda(1));
      hess(0,1)=t*exp(lambda(0)+t*lambda(1));hess(1,0)=hess(0,1);
      hess(0,2)=0;hess(2,0)=hess(0,2);
      hess(1,1)=pow((double) t,(double) 2)*exp(lambda(0)+t*lambda(1));
      hess(1,2)=0;hess(2,1)=hess(1,2);
      hess(2,2)=0;
      return(hess);
    }
  }else{
    return(0);
  }
}

/*Everything required for ML-estimation in this part*/
// [[Rcpp::export]]
double mlFirst(NumericVector y, NumericMatrix groupE, NumericMatrix groupC, int nE, int nC, NumericVector tpE, NumericVector tpC, int type){
  /* y corresponds to the parameter which are to be estimated. The last entry
  corresponds to the size, while the others correspond to the individual time points.
  Further structure of y: Last entry = size. Before each entry depends on the
  chosen trend */

  int j,t,ny; /*j for patients; t for time points*/
  double logL; /*log likelihood*/

  ny = y.length();

  //First group
  logL=0;
  for(j=0;j<nE;j++){
    for(t=0;t<tpE(j);t++){
      logL += ldnbinom(groupE(j,t), trend(y, 1, t, type), y(ny-1));
    }
  }
  //Second group
  for(j=0;j<nC;j++){
    for(t=0;t<tpC(j);t++){
      logL += ldnbinom(groupC(j,t), trend(y, 2, t, type), y(ny-1));
    }
  }
  return(-1/((float) nE+ (float) nC)*logL);
}

// [[Rcpp::export]]
double mlFirstOneGroup(NumericVector y, NumericMatrix groupC, int nC, NumericVector tpC, int type){
  /* y corresponds to the parameter which are to be estimated. The last entry
  corresponds to the size, while the others correspond to the individual time points.
  Further structure of y: Last entry = size. Before each entry depends on the
  chosen trend */

  int j,t,ny; /*j for patients; t for time points*/
  double logL; /*log likelihood*/

  ny = y.length();

  logL=0;
  for(j=0;j<nC;j++){
    for(t=0;t<tpC(j);t++){
      logL += ldnbinom(groupC(j,t), trend(y, 2, t, type), y(ny-1));
    }
  }
  return(-1/((float) nC)*logL);
}

// [[Rcpp::export]]
double mlSecond(double rho, NumericVector y, NumericMatrix groupE, NumericMatrix groupC, int nE, int nC, NumericVector tpE, NumericVector tpC, int type){

  int j,t,s,nY,r,l; /*j for patients; s,t for time points; k,l for specific calculations*/
  double logL,term1,term2,term3,term4,pairwiseProb; /*log likelihood; term1-4 specific terms in sum*/

  nY = y.length();

  //First group
  logL=0;
  for(j=0;j<nE;j++){
    for(s=0;s<(tpE(j)-1);s++){
      for(t=s+1;t<tpE(j);t++){
        pairwiseProb = 0;
        for(r=0;r<=groupE(j,s);r++){
          for(l=0;l<=groupE(j,t);l++){
            term1 = dnbinom(r, trend(y, 1, s, type)*(1-pow(rho, abs(s-t))), y(nY-1)*(1-pow(rho, abs(s-t))));
            term2 = dnbinom(l, trend(y, 1, t, type)*(1-pow(rho, abs(s-t))), y(nY-1)*(1-pow(rho, abs(s-t))));
            term3 = dnbinom(groupE(j,t)+groupE(j,s) - r - l, (trend(y, 1, s, type)+trend(y, 1, t, type))*pow(rho, abs(s-t)), y(nY-1)*pow(rho, abs(s-t)));
            term4 = dbinom(groupE(j,s)-r, groupE(j,s)+groupE(j,t)-r-l, trend(y, 1, s, type)/(trend(y, 1, s, type)+trend(y, 1, t, type)));

            pairwiseProb += term1*term2*term3*term4;
          }
        }
        logL += log(pairwiseProb);
      }
    }
  }

  //Second Group
  for(j=0;j<nC;j++){
    for(s=0;s<(tpC(j)-1);s++){
      for(t=s+1;t<tpC(j);t++){
        pairwiseProb = 0;
        for(r=0;r<=groupC(j,s);r++){
          for(l=0;l<=groupC(j,t);l++){
            term1 = dnbinom(r, trend(y, 2, s, type)*(1-pow(rho, abs(s-t))), y(nY-1)*(1-pow(rho, abs(s-t))));
            term2 = dnbinom(l, trend(y, 2, t, type)*(1-pow(rho, abs(s-t))), y(nY-1)*(1-pow(rho, abs(s-t))));
            term3 = dnbinom(groupC(j,t)+groupC(j,s) - r - l, (trend(y, 2, s, type)+trend(y, 2, t, type))*pow(rho, abs(s-t)), y(nY-1)*pow(rho, abs(s-t)));
            term4 = dbinom(groupC(j,s)-r, groupC(j,s)+groupC(j,t)-r-l, trend(y, 2, s, type)/(trend(y, 2, s, type)+trend(y, 2, t, type)));

            pairwiseProb += term1*term2*term3*term4;

          }
        }
        logL += log(pairwiseProb);
      }
    }
  }
  return(-1/((float) nE+ (float) nC)*logL);
}

// [[Rcpp::export]]
double mlSecondOneGroup(double rho, NumericVector y, NumericMatrix groupC, int nC, NumericVector tpC, int type){

  int j,t,s,nY,r,l; /*j for patients; s,t for time points; k,l for specific calculations*/
  double logL,term1,term2,term3,term4,pairwiseProb; /*log likelihood; term1-4 specific terms in sum*/

  nY = y.length();

  logL=0;
  for(j=0;j<nC;j++){
    for(s=0;s<(tpC(j)-1);s++){
      for(t=s+1;t<tpC(j);t++){
        pairwiseProb = 0;
        for(r=0;r<=groupC(j,s);r++){
          for(l=0;l<=groupC(j,t);l++){
            term1 = dnbinom(r, trend(y, 2, s, type)*(1-pow(rho, abs(s-t))), y(nY-1)*(1-pow(rho, abs(s-t))));
            term2 = dnbinom(l, trend(y, 2, t, type)*(1-pow(rho, abs(s-t))), y(nY-1)*(1-pow(rho, abs(s-t))));
            term3 = dnbinom(groupC(j,t)+groupC(j,s) - r - l, (trend(y, 2, s, type)+trend(y, 2, t, type))*pow(rho, abs(s-t)), y(nY-1)*pow(rho, abs(s-t)));
            term4 = dbinom(groupC(j,s)-r, groupC(j,s)+groupC(j,t)-r-l, trend(y, 2, s, type)/(trend(y, 2, s, type)+trend(y, 2, t, type)));

            pairwiseProb += term1*term2*term3*term4;

          }
        }
        logL += log(pairwiseProb);
      }
    }
  }
  return(-1/((float) nC)*logL);
}

// [[Rcpp::export]]
NumericVector mlFirstGrad(NumericVector y, NumericMatrix groupE, NumericMatrix groupC, int nE, int nC, NumericVector tpE, NumericVector tpC, int type){
  /* y corresponds to the parameter which are to be estimated. The last entry
  corresponds to the size, while the others correspond to the individual time points.
  Further structure of y: Last entry = size. Before each entry depends on the
  chosen trend */
  int j,t,k,nG,nL;
  //nG is length of gradient corresponding to number of parameters; nL is number of parameters in lambda
  nG=y.length();
  nL=nG-1;
  NumericVector Jhelp(nG);
  NumericVector trendG(nL);
  NumericVector lambda(nL);
  double eta;

  for(k=0;k<nL;k++){
    lambda(k)=y(k);
  }
  eta=y(nG-1);

  for(j=0;j<nE;j++){
    for(t=0;t<tpE(j);t++){
      trendG = trendGrad(lambda, 1, t, type);
      for(k=0;k<nL;k++){
        Jhelp(k)+=ldnbinomDmu(groupE(j,t), trend(lambda, 1, t, type), eta)*trendG(k);
      }
      Jhelp(nG-1)+=ldnbinomDeta(groupE(j,t), trend(lambda, 1, t, type), eta);
    }
  }
  for(j=0;j<nC;j++){
    for(t=0;t<tpC(j);t++){
      trendG = trendGrad(lambda, 2, t, type);
      for(k=0;k<nL;k++){
        Jhelp(k)+=ldnbinomDmu(groupC(j,t), trend(lambda, 2, t, type), eta)*trendG(k);
      }
      Jhelp(nG-1)+=ldnbinomDeta(groupC(j,t), trend(lambda, 2, t, type), eta);
    }
  }
  return(-1/((float) nE + (float) nC)*Jhelp);
}

// [[Rcpp::export]]
NumericMatrix mlFirstHObs(NumericVector y, NumericMatrix groupE, NumericMatrix groupC, int nE, int nC, NumericVector tpE, NumericVector tpC, int type){
  int j,t,k,l,nG,nL;
  double divN;
  //nG is length of gradient corresponding to number of parameters; nL is number of parameters in lambda
  nG=y.length();
  nL=nG-1;

  NumericMatrix H(nG,nG);
  NumericMatrix trendH(nL,nL);
  NumericVector trendG(nL);
  NumericVector lambda(nL);
  double eta;

  for(k=0;k<nL;k++){
    lambda(k)=y(k);
  }
  eta=y(nG-1);

  for(j=0;j<nE;j++){
    for(t=0;t<tpE(j);t++){
      //Fill Hessian with DlambdaDlambda
      trendG = trendGrad(lambda, 1, t, type);
      trendH = trendHess(lambda, 1, t, type);
      for(k=0;k<nL;k++){
        for(l=0;l<nL;l++){
          H(k,l)+=ldnbinomDmuDmu(groupE(j,t), trend(lambda, 1, t, type), eta)*trendG(k)*trendG(l)+ldnbinomDmu(groupE(j,t), trend(lambda, 1, t, type), eta)*trendH(k,l);
        }
      }
      //Fill Hessian with DlambdaDeta
      for(k=0;k<nL;k++){
        H(k, nL)+=ldnbinomDmuDeta(groupE(j,t), trend(lambda, 1, t, type), eta)*trendG(k);
        H(nL, k)=H(k,nL);
      }
      //..and DetaDeta
      H(nL, nL)+=ldnbinomDetaDeta(groupE(j,t), trend(lambda, 1, t, type), eta);
    }
  }
  for(j=0;j<nC;j++){
    for(t=0;t<tpC(j);t++){
      //Fill Hessian with DlambdaDlambda
      trendG = trendGrad(lambda, 2, t, type);
      trendH = trendHess(lambda, 2, t, type);
      for(k=0;k<nL;k++){
        for(l=0;l<nL;l++){
          H(k,l)+=ldnbinomDmuDmu(groupC(j,t), trend(lambda, 2, t, type), eta)*trendG(k)*trendG(l)+ldnbinomDmu(groupC(j,t), trend(lambda, 2, t, type), eta)*trendH(k,l);
        }
      }
      //Fill Hessian with DlambdaDeta
      for(k=0;k<nL;k++){
        H(k, nL)+=ldnbinomDmuDeta(groupC(j,t), trend(lambda, 2, t, type), eta)*trendG(k);
        H(nL, k)=H(k,nL);
      }
      //..and DetaDeta
      H(nL, nL)+=ldnbinomDetaDeta(groupC(j,t), trend(lambda, 2, t, type), eta);
    }
  }
  divN = (float)nE + (float)nC;
  for(j=0;j<nG;j++){
    for(k=0;k<nG;k++){
      H(j,k)=H(j,k)/divN;
    }
  }
  return(H);
}

// [[Rcpp::export]]
NumericMatrix mlFirstHExp(NumericVector y, double kf, int tp, int type){
  int t,k,l,nG,nL;
  //nG is length of gradient corresponding to number of parameters; nL is number of parameters in lambda
  nG=y.length();
  nL=nG-1;

  NumericMatrix HE(nG,nG);
  NumericMatrix HC(nG,nG);
  NumericMatrix H(nG, nG);
  NumericVector trendG(nL);
  NumericVector lambda(nL);
  double eta;

  for(k=0;k<nL;k++){
    lambda(k)=y(k);
  }
  eta=y(nG-1);

  //Fill Hessian with DlambdaDlambda
  for(t=0;t<tp;t++){
    trendG = trendGrad(lambda, 1, t, type);
    for(k=0;k<nL;k++){
      for(l=0;l<nL;l++){
        HE(k,l)+=ldnbinomDmuDmuExp(trend(lambda, 1, t, type), eta)*trendG(k)*trendG(l);
      }
    }
  }
  //Fill Hessian with DlambdaDeta
  for(k=0;k<nL;k++){
    HE(k, nL)=0;
    HE(nL, k)=0;
  }
  //..and DetaDeta
  for(t=0;t<tp;t++){
    HE(nL, nL)+=1/eta-1/(trend(lambda, 1, t, type)+eta)+ExpTerm1(trend(lambda, 1, t, type), eta)-trigamma(eta);
  }

  //Fill Hessian with DlambdaDlambda
  for(t=0;t<tp;t++){
    trendG = trendGrad(lambda, 2, t, type);
    for(k=0;k<nL;k++){
      for(l=0;l<nL;l++){
        HC(k,l)+=ldnbinomDmuDmuExp(trend(lambda, 2, t, type), eta)*trendG(k)*trendG(l);
      }
    }
  }
  //Fill Hessian with DlambdaDeta
  for(k=0;k<nL;k++){
    HC(k, nL)=0;
    HC(nL, k)=0;
  }
  //..and DetaDeta
  for(t=0;t<tp;t++){
    HC(nL, nL)+=1/eta-1/(trend(lambda, 2, t, type)+eta)+ExpTerm1(trend(lambda, 2, t, type), eta)-trigamma(eta);
  }

  for(k=0;k<nG;k++){
    for(l=0;l<nG;l++){
      H(k,l)=-kf/(kf+1)*HE(k,l)-1/(kf+1)*HC(k,l);
    }
  }

  return(H);
}

// [[Rcpp::export]]
NumericMatrix mlFirstJObs(NumericVector y, NumericMatrix groupE, NumericMatrix groupC, int nE, int nC, NumericVector tpE, NumericVector tpC, int type){
  int j,t,k,l,nG,nL;
  double divN;

  nG=y.length();
  nL=nG-1;

  NumericMatrix J(nG,nG);
  NumericVector Jhelp(nG);
  NumericVector trendG(nL);
  NumericVector lambda(nL);
  double eta;

  for(k=0;k<nL;k++){
    lambda(k)=y(k);
  }
  eta=y(nG-1);

  for(j=0;j<nE;j++){
    for(k=0;k<nG;k++){
      Jhelp(k)=0;
    }
    for(t=0;t<tpE(j);t++){
      trendG = trendGrad(lambda, 1, t, type);
      for(k=0;k<nL;k++){
        Jhelp(k)+=ldnbinomDmu(groupE(j,t), trend(lambda, 1, t, type), eta)*trendG(k);
      }
      Jhelp(nL)+=ldnbinomDeta(groupE(j,t), trend(lambda, 1, t, type), eta);
    }
    for(k=0;k<nG;k++){
      for(l=k;l<nG;l++){
        J(k,l)+=Jhelp(k)*Jhelp(l);
        J(l,k)=J(k,l);
      }
    }
  }
  for(j=0;j<nC;j++){
    for(k=0;k<nG;k++){
      Jhelp(k)=0;
    }
    for(t=0;t<tpC(j);t++){
      trendG = trendGrad(lambda, 2, t, type);
      for(k=0;k<nL;k++){
        Jhelp(k)+=ldnbinomDmu(groupC(j,t), trend(lambda, 2, t, type), eta)*trendG(k);
      }
      Jhelp(nL)+=ldnbinomDeta(groupC(j,t), trend(lambda, 2, t, type), eta);
    }
    for(k=0;k<nG;k++){
      for(l=k;l<nG;l++){
        J(k,l)+=Jhelp(k)*Jhelp(l);
        J(l,k)=J(k,l);
      }
    }
  }
  divN=((float)nE+(float)nC);
  for(j=0;j<nG;j++){
    for(k=0;k<nG;k++){
      J(j,k)=J(j,k)/divN;
    }
  }
  return(J);
}

// [[Rcpp::export]]
NumericMatrix mlFirstJExp(NumericVector y, double rho, double kf, int tp, int type, int approx = 50){
  //Dieses rho entspricht lediglich dem autokorrelationskoeffizienten
  int t,s,k,l,nG,nL;

  nG=y.length();
  nL=nG-1;

  NumericMatrix J(nG,nG);
  NumericMatrix JE(nG, nG);
  NumericMatrix JC(nG, nG);

  NumericVector trendGEt(nL);
  NumericVector trendGEs(nL);
  NumericVector trendGCt(nL);
  NumericVector trendGCs(nL);

  NumericVector ExpTerm1AllE(tp);
  NumericVector ExpTerm1AllC(tp);
  NumericMatrix ExpTerm2AllE(tp, tp);
  NumericMatrix ExpTerm2AllC(tp, tp);
  NumericMatrix ExpTerm3AllE(tp, tp);
  NumericMatrix ExpTerm3AllC(tp, tp);

  NumericVector lambda(nL);
  double eta, fEt, fEs, fCt, fCs;

  //Calculate helpmatrices
  for(k=0;k<nL;k++){
    lambda(k)=y(k);
  }
  eta=y(nG-1);

  for(t=0;t<tp;t++){
    ExpTerm1AllE(t) = ExpTerm1(trend(lambda, 1, t, type), eta, approx);
    ExpTerm1AllC(t) = ExpTerm1(trend(lambda, 2, t, type), eta, approx);
  }

  for(t=0;t<tp;t++){
    for(s=0;s<tp;s++){
      ExpTerm2AllE(t, s)=ExpTerm2(trend(lambda, 1, t, type), trend(lambda, 1, s, type), eta, pow(rho, abs(s-t)), approx);
      ExpTerm2AllC(t, s)=ExpTerm2(trend(lambda, 2, t, type), trend(lambda, 2, s, type), eta, pow(rho, abs(s-t)), approx);
      if(s >= t){
        ExpTerm3AllE(t, s)=ExpTerm3(trend(lambda, 1, t, type), trend(lambda, 1, s, type), eta, pow(rho, abs(s-t)), approx);
        ExpTerm3AllC(t, s)=ExpTerm3(trend(lambda, 2, t, type), trend(lambda, 2, s, type), eta, pow(rho, abs(s-t)), approx);
      }else{
        ExpTerm3AllE(t, s) = ExpTerm3AllE(s, t);
        ExpTerm3AllC(t, s) = ExpTerm3AllC(s, t);
      }
    }
  }

  for(t=0;t<tp;t++){
    trendGEt = trendGrad(lambda, 1, t, type);
    trendGCt = trendGrad(lambda, 2, t, type);
    for(s=0;s<tp;s++){
      trendGEs = trendGrad(lambda, 1, s, type);
      trendGCs = trendGrad(lambda, 2, s, type);
      fEt = trend(lambda, 1, t, type); fEs = trend(lambda, 1, s, type);
      fCt = trend(lambda, 2, t, type); fCs = trend(lambda, 2, s, type);
      //Fill J with Deta X Deta
      if(s==t){
        JE(nL, nL) += ExpTerm1AllE(s)*(1+log(eta/(fEt+eta))-digamma(eta)-eta/(fEt+eta)) + ExpTerm1AllE(t)*(1+log(eta/(fEs+eta))-digamma(eta)-eta/(fEs+eta)) - ExpTerm2AllE(t, s)/(fEt+eta) - ExpTerm2AllE(s, t)/(fEs+eta)+ExpTerm3AllE(t, s) + (log(eta/(fEt+eta))-digamma(eta))*(log(eta/(fEs+eta))-digamma(eta))+fEt/(eta*(fEt+eta));
        JC(nL, nL) += ExpTerm1AllC(s)*(1+log(eta/(fCt+eta))-digamma(eta)-eta/(fCt+eta)) + ExpTerm1AllC(t)*(1+log(eta/(fCs+eta))-digamma(eta)-eta/(fCs+eta)) - ExpTerm2AllC(t, s)/(fCt+eta) - ExpTerm2AllC(s, t)/(fCs+eta)+ExpTerm3AllC(t, s) + (log(eta/(fCt+eta))-digamma(eta))*(log(eta/(fCs+eta))-digamma(eta))+fCt/(eta*(fCt+eta));
      }else{
        JE(nL, nL) += ExpTerm1AllE(s)*(1+log(eta/(fEt+eta))-digamma(eta)-eta/(fEt+eta)) + ExpTerm1AllE(t)*(1+log(eta/(fEs+eta))-digamma(eta)-eta/(fEs+eta)) - ExpTerm2AllE(t, s)/(fEt+eta) - ExpTerm2AllE(s, t)/(fEs+eta)+ExpTerm3AllE(t, s) + (log(eta/(fEt+eta))-digamma(eta))*(log(eta/(fEs+eta))-digamma(eta))+pow(rho, abs(t-s))*fEt*fEs/(eta*(fEt+eta)*(fEs+eta));
        JC(nL, nL) += ExpTerm1AllC(s)*(1+log(eta/(fCt+eta))-digamma(eta)-eta/(fCt+eta)) + ExpTerm1AllC(t)*(1+log(eta/(fCs+eta))-digamma(eta)-eta/(fCs+eta)) - ExpTerm2AllC(t, s)/(fCt+eta) - ExpTerm2AllC(s, t)/(fCs+eta)+ExpTerm3AllC(t, s) + (log(eta/(fCt+eta))-digamma(eta))*(log(eta/(fCs+eta))-digamma(eta))+pow(rho, abs(t-s))*fCt*fCs/(eta*(fCt+eta)*(fCs+eta));
      }
      for(k=0;k<nL;k++){
        //Fill J with Dlambda X Deta
        if(s==t){
          JE(k, nL) += -1/(fEt+eta)*trendGEs(k) + (ExpTerm2AllE(s, t)*eta/(fEs*(fEs+eta))-ExpTerm1AllE(t)*eta/(fEs+eta))*trendGEs(k);
          JC(k, nL) += -1/(fCt+eta)*trendGCs(k) + (ExpTerm2AllC(s, t)*eta/(fCs*(fCs+eta))-ExpTerm1AllC(t)*eta/(fCs+eta))*trendGCs(k);
          JE(nL, k) = JE(k, nL);
          JC(nL, k) = JC(k, nL);
        }else{
          JE(k, nL) += -pow(rho, abs(t-s))*fEt/((fEt+eta)*(fEs+eta))*trendGEs(k) + (ExpTerm2AllE(s, t)*eta/(fEs*(fEs+eta))-ExpTerm1AllE(t)*eta/(fEs+eta))*trendGEs(k);
          JC(k, nL) += -pow(rho, abs(t-s))*fCt/((fCt+eta)*(fCs+eta))*trendGCs(k) + (ExpTerm2AllC(s, t)*eta/(fCs*(fCs+eta))-ExpTerm1AllC(t)*eta/(fCs+eta))*trendGCs(k);
          JE(nL, k) = JE(k, nL);
          JC(nL, k) = JC(k, nL);
        }
        for(l=0;l<nL;l++){
          //Fill J with Dlambda X Dlambda
          if(s==t){
            JE(k,l) += eta/(fEt*(fEt+eta))*trendGEt(k)*trendGEt(l);
            JC(k,l) += eta/(fCt*(fCt+eta))*trendGCt(k)*trendGCt(l);
          }else{
            JE(k,l) += pow(rho, abs(t-s))*eta/((fEt + eta)*(fEs + eta))*trendGEs(l)*trendGEt(k);
            JC(k,l) += pow(rho, abs(t-s))*eta/((fCt + eta)*(fCs + eta))*trendGCs(l)*trendGCt(k);
          }
        }
      }
    }
  }
  for(k=0;k<nG;k++){
    for(l=0;l<nG;l++){
      J(k,l)=kf/(kf+1)*JE(k,l)+1/(kf+1)*JC(k,l);
    }
  }
  return(J);
}

/* Blinded Estimation of parameters*/
// [[Rcpp::export]]
double mlFirstBlinded(NumericVector y, NumericMatrix group, int n, NumericVector tp, int type, double theta, double k){
  //Type 1 constant trend y=(lambda_1, eta)
  //Type 2 exponential trend y=c(lambda_1,lambda_2, eta)

  int j,t; /*j for patients; t for time points*/
int nY;
double logL, eta; /*log likelihood*/

nY = y.length();
NumericVector lambda(nY);
eta=y(nY-1);

if(type==1){
  lambda(0)=y(0);
  lambda(1)=theta;
}else if(type ==2){
  lambda(0)=y(0);
  lambda(1)=y(1);
  lambda(2)=theta;
}

logL=0;
for(j=0;j<n;j++){
  for(t=0;t<tp(j);t++){
    logL += log(k/(1+k)*dnbinom(group(j,t), trend(lambda, 1, t, type), eta)+1/(1+k)*dnbinom(group(j,t), trend(lambda, 2, t, type), eta));
  }
}
return(-1/((float) n)*logL);
}

// [[Rcpp::export]]
double mlSecondBlinded(double rho, NumericVector y, NumericMatrix group, int n, NumericVector tp, int type, double k){

  int j,t,s,nY,r,l; /*j for patients; s,t for time points; k,l for specific calculations*/
double logL,term1,term2,term3,term4,pairwiseProb; /*log likelihood; term1-4 specific terms in sum*/

nY = y.length();

logL=0;
for(j=0;j<n;j++){
  for(s=0;s<(tp(j)-1);s++){
    for(t=s+1;t<tp(j);t++){
      pairwiseProb = 0;
      for(r=0;r<=group(j,s);r++){
        for(l=0;l<=group(j,t);l++){
          term1 = dnbinom(r, trend(y, 1, s, type)*(1-pow(rho, abs(s-t))), y(nY-1)*(1-pow(rho, abs(s-t))));
          term2 = dnbinom(l, trend(y, 1, t, type)*(1-pow(rho, abs(s-t))), y(nY-1)*(1-pow(rho, abs(s-t))));
          term3 = dnbinom(group(j,t)+group(j,s) - r - l, (trend(y, 1, s, type)+trend(y, 1, t, type))*pow(rho, abs(s-t)), y(nY-1)*pow(rho, abs(s-t)));
          term4 = dbinom(group(j,s)-r, group(j,s)+group(j,t)-r-l, trend(y, 1, s, type)/(trend(y, 1, s, type)+trend(y, 1, t, type)));

          pairwiseProb += k/(1+k)*term1*term2*term3*term4;

          term1 = dnbinom(r, trend(y, 2, s, type)*(1-pow(rho, abs(s-t))), y(nY-1)*(1-pow(rho, abs(s-t))));
          term2 = dnbinom(l, trend(y, 2, t, type)*(1-pow(rho, abs(s-t))), y(nY-1)*(1-pow(rho, abs(s-t))));
          term3 = dnbinom(group(j,t)+group(j,s) - r - l, (trend(y, 2, s, type)+trend(y, 2, t, type))*pow(rho, abs(s-t)), y(nY-1)*pow(rho, abs(s-t)));
          term4 = dbinom(group(j,s)-r, group(j,s)+group(j,t)-r-l,trend(y, 2, s, type)/(trend(y, 2, s, type)+trend(y, 2, t, type)));

          pairwiseProb += 1/(1+k)*term1*term2*term3*term4;
        }
      }
      logL += log(pairwiseProb);
    }
  }
}
return(-logL);
}
