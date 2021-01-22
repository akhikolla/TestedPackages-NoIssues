#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector cusum(NumericVector x){
  int n = x.size();
  int iter;
  NumericVector out(n);
  NumericVector I_plus(n);
  NumericVector I_minus(n);
  double sum_of_x, i_a, I_plusinv,fct,npw2,n_inv,n_a=n;
  i_a=double(0);
  I_plusinv=double(0);
  fct=0;
  n_inv=1/n_a;
  npw2=n_a*n_a;
  sum_of_x=0;
  for (iter=1;iter < n; iter++) sum_of_x +=x[iter];
  I_minus[0] = 1/sqrt(npw2-n_a) * sum_of_x;
  I_plus[0]=sqrt(1-n_inv)*x[0];
  out[0]=I_plus[0] - I_minus[0];
  
  for(iter=1;iter<n-1;iter++){
    i_a=(double)iter;
    I_plusinv=1/(i_a+1);
    fct=sqrt((n_a-i_a-1)*i_a*I_plusinv/(n_a-i_a));
    I_plus[iter]=x[iter]*sqrt(I_plusinv - n_inv)+I_plus[iter-1]*fct;
    I_minus[iter]=I_minus[iter-1]/fct - x[iter]/sqrt(npw2*I_plusinv - n_a);
    out[iter]=(I_plus[iter] - I_minus[iter]);
  }
  for (iter=0;iter < n-1; iter++) out[iter]=out[iter]/pow(sum_of_x/n_a,1.0);//divide by the mean
  return(out);
}

//[[Rcpp::export]]
NumericVector finner_prod_maxp(NumericVector x,double p){
  int n = x.size();
  int max_b=0;
  int iter,s0,e0;
  double aux,maxipi,max_inner;
  max_inner=0;
  s0=floor((1-p)*n);
  e0=ceil(p*n);
  maxipi=0;
  for(iter=s0-1;iter<e0-1;iter++){
    aux=fabs(x[iter]);
    if(aux>maxipi){
      max_b=iter+1;
      maxipi=aux;
    }
  }
  max_inner=x[max_b-1];
  NumericVector ret(2);
  ret(0)=max_inner;
  ret(1)=max_b+1;
  return ret;
  }

//[[Rcpp::export]]
NumericVector across_fip(NumericMatrix X, NumericVector tau, double p, NumericVector epp, double p1, double Ts){
  int n_row=X.nrow();
  int n_col=X.ncol();
  bool TF = false;
  NumericVector aux_vec(n_row-1);
  NumericVector out(n_row -1);
  for (int j=0;j<n_col-1;j++){
    aux_vec=cusum(X(_,j));
    for (int i=0;i<n_row-1;i++){
      TF=p1<epp(j);
      if (TF==true){
        out(i)=out(i)+0;
        continue;
      }
      if (fabs(aux_vec(i))>tau[j]*Ts) out(i)=out(i)+fabs(aux_vec(i));
    }
  }
  return(finner_prod_maxp(out,p));
}

//[[Rcpp::export]]
List multi_across_fip(NumericMatrix X, int M, double min_draw, NumericVector tau,NumericVector p,NumericVector epp,double Ts){
  int len = X.nrow();
  int d = X.ncol();
  double p1,p2;
  NumericVector out1(M),out2(M),out3(M),out4(M),aux_fip;
  List ret;
  for (int m=0; m<M; m++){
    p1=ceil(R::runif(1,(len-2)+1));
    p2=ceil(R::runif(p1,len));
    if(m==0){
      p1=1+epp(d-1);
      p2=len-epp(d-1);
    }
    if ((p2-p1)<min_draw) continue;
    out3[m]=p1;
    out4[m]=p2;
    if(m==0){
      aux_fip=across_fip(X(Range(p1,p2),_),tau,p(0),epp,p1,Ts);;
    } else aux_fip=across_fip(X(Range(p1,p2),_),tau,p(1),epp,p1,Ts);
    out1[m]=aux_fip(1)+p1-1;
    out2[m]=aux_fip(0);
  }
  ret["sup_b"]=out1;
  ret["max_inner"]=out2;
  ret["p1_mat"]=out3;
  ret["p2_mat"]=out4;
  return ret;
}
