#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]

NumericMatrix g1(NumericVector T, int order)
{
    NumericMatrix M(order,1);
    NumericVector tempT=NumericVector::create(0,T[2]/pow(T[1],2),-1/T[1],0);

    if(order==2) {
        M(0,0)=tempT[1];
        M(1,0)=tempT[2];
    }
    if(order==3) {
        M(0,0)=0;
        M(1,0)=tempT[1];
        M(2,0)=tempT[2];
    }
    if(order==4){
        M(0,0)=tempT[0];
        M(1,0)=tempT[1];
        M(2,0)=tempT[2];
        M(3,0)=tempT[3];
    }
    return M;
}

// [[Rcpp::export]]

NumericMatrix g2(NumericVector T, double dt, int order)
{
    NumericMatrix M(order,1);
    NumericVector tempT(T.size());

    if(T[1]<0)
        tempT=NumericVector::create(1/((T[0]-dt)*T[1]),(-log((T[0]-dt)/dt)+T[2])/pow(T[1],2),-1/T[1],0);
    else
        tempT=NumericVector::create(-1/((T[0]+dt)*T[1]),(-log(-dt/(T[0]+dt))+T[2])/pow(T[1],2),-1/T[1],0);

    if(order==2) {
        M(0,0)=tempT[1];
        M(1,0)=tempT[2];
    }
    if(order==3) {
        M(0,0)=tempT[0];
        M(1,0)=tempT[1];
        M(2,0)=tempT[2];
    }
    if(order==4){
        M(0,0)=tempT[0];
        M(1,0)=tempT[1];
        M(2,0)=tempT[2];
        M(3,0)=tempT[3];
    }
    return M;
}

// [[Rcpp::export]]

NumericMatrix f234(NumericVector T, double x, int order)
{
    NumericMatrix M(order,1);

    if(order==2) {
        M(0,0)=exp((T[1]*x+T[2])/2)/(1+exp(T[1]*x+T[2]));
        M(1,0)=x*exp((T[1]*x+T[2])/2)/(1+exp(T[1]*x+T[2]));
    }
    if(order==3) {
        M(0,0)=1/(1+exp(T[1]*x+T[2]));
        M(1,0)=(-T[0]*x*exp(T[1]*x+T[2]))/pow(1+exp(T[1]*x+T[2]),2);
        M(2,0)=(-T[0]*exp(T[1]*x+T[2]))/pow(1+exp(T[1]*x+T[2]),2);
    }
    if(order==4) {
        M(0,0)=1/(1+exp(T[1]*x+T[2]));
        M(1,0)=(-T[0]*x*exp(T[1]*x+T[2]))/pow(1+exp(T[1]*x+T[2]),2);
        M(2,0)=(-T[0]*exp(T[1]*x+T[2]))/pow(1+exp(T[1]*x+T[2]),2);
        M(3,0)=1;
    }
    return M;
}

// [[Rcpp::export]]

double SDM(NumericMatrix M)
{
    double sum=0;
    for(int i=0;i<M.nrow();i++)
        sum=sum+M(i,i);

    return sum;
}

// [[Rcpp::export]]

NumericMatrix Trans(NumericMatrix M)
{
    int nrow=M.ncol(),ncol=M.nrow();
    NumericMatrix N(nrow,ncol);

    for(int i=0;i<nrow;i++)
        for(int j=0;j<ncol;j++)
        N(i,j)=M(j,i);

    return N;
}

// [[Rcpp::export]]

NumericMatrix Minus(NumericMatrix M1, NumericMatrix M2)
{
    int nrow=M1.nrow(),ncol=M1.ncol();
    NumericMatrix M3(nrow,ncol);

    for(int i=0;i<nrow;i++)
        for(int j=0;j<ncol;j++)
            M3(i,j)=M1(i,j)-M2(i,j);

    return M3;
}

// [[Rcpp::export]]

NumericMatrix Plus(NumericMatrix M1, NumericMatrix M2)
{
    int nrow=M1.nrow(),ncol=M1.ncol();
    NumericMatrix M3(nrow,ncol);

    for(int i=0;i<nrow;i++)
        for(int j=0;j<ncol;j++)
            M3(i,j)=M1(i,j)+M2(i,j);

    return M3;
}

// [[Rcpp::export]]

NumericMatrix Multiple(NumericMatrix M1, NumericMatrix M2)
{
    int nrow=M1.nrow(),ncol=M2.ncol(),P=M1.ncol();
    NumericMatrix M3(nrow,ncol);
    double sum;

    for(int i=0;i<nrow;i++)
        for(int j=0;j<ncol;j++) {
            sum=0;
            for(int k=0;k<P;k++)
                sum=sum+M1(i,k)*M2(k,j);
            M3(i,j)=sum;
        }

    return M3;
}

// [[Rcpp::export]]

NumericMatrix infor234(NumericVector T, double x, int order)
{
    NumericMatrix M=f234(T,x,order);

    return Multiple(M,Trans(M));
}

// [[Rcpp::export]]

double d1(NumericVector T, double x, double xl, NumericMatrix inv, int order)
{
    NumericMatrix M=Multiple(inv,Minus(infor234(T,x,order),infor234(T,xl,order)));
    double s=SDM(M);

    return s;
}


// [[Rcpp::export]]

double d2(NumericVector T, double x, double xl, NumericMatrix inv, int order)
{
    double s1=Multiple(Multiple(Multiple(Multiple(Trans(g1(T,order)),inv),Minus(infor234(T,x,order),infor234(T,xl,order))),inv),g1(T,order))(0,0);
    double s2=pow(Multiple(Multiple(Trans(g1(T,order)),inv),g1(T,order))(0,0),-1);

    return -s1*s2;
}


// [[Rcpp::export]]

double d3(NumericVector T, double x, double xl, NumericMatrix inv, double dt, int order)
{
    double s1=Multiple(Multiple(Multiple(Multiple(Trans(g2(T,dt,order)),inv),Minus(infor234(T,x,order),infor234(T,xl,order))),inv),g2(T,dt,order))(0,0);
    double s2=pow(Multiple(Multiple(Trans(g2(T,dt,order)),inv),g2(T,dt,order))(0,0),-1);

    return -s1*s2;
}

// [[Rcpp::export]]

double dd1(NumericVector T, double x1, double x2, double xl, NumericMatrix inv, int order)
{
    NumericMatrix M=Multiple(Multiple(Multiple(inv,Minus(infor234(T,x2,order),infor234(T,xl,order))),inv),Minus(infor234(T,x1,order),infor234(T,xl,order)));
    double s=SDM(M);

    return -s;
}

// [[Rcpp::export]]

double dd2(NumericVector T, double x1, double x2, double xl, NumericMatrix inv, int order)
{
    double s1=Multiple(Multiple(Multiple(Multiple(Multiple(Multiple(Trans(g1(T,order)),inv),Minus(infor234(T,x2,order),infor234(T,xl,order))),inv),Minus(infor234(T,x1,order),infor234(T,xl,order))),inv),g1(T,order))(0,0);
    double s2=Multiple(Multiple(Multiple(Multiple(Multiple(Multiple(Trans(g1(T,order)),inv),Minus(infor234(T,x1,order),infor234(T,xl,order))),inv),Minus(infor234(T,x2,order),infor234(T,xl,order))),inv),g1(T,order))(0,0);
    double s3=Multiple(Multiple(Trans(g1(T,order)),inv),g1(T,order))(0,0);
    double s4=Multiple(Multiple(Multiple(Multiple(Trans(g1(T,order)),inv),Minus(infor234(T,x1,order),infor234(T,xl,order))),inv),g1(T,order))(0,0);
    double s5=Multiple(Multiple(Multiple(Multiple(Trans(g1(T,order)),inv),Minus(infor234(T,x2,order),infor234(T,xl,order))),inv),g1(T,order))(0,0);
    double s6=pow(Multiple(Multiple(Trans(g1(T,order)),inv),g1(T,order))(0,0),-2);

    return (s1+s2*s3-s4*s5)*s6;
}

// [[Rcpp::export]]

double dd3(NumericVector T, double x1, double x2, double xl, NumericMatrix inv, double dt, int order)
{
    double s1=Multiple(Multiple(Multiple(Multiple(Multiple(Multiple(Trans(g2(T,dt,order)),inv),Minus(infor234(T,x2,order),infor234(T,xl,order))),inv),Minus(infor234(T,x1,order),infor234(T,xl,order))),inv),g2(T,dt,order))(0,0);
    double s2=Multiple(Multiple(Multiple(Multiple(Multiple(Multiple(Trans(g2(T,dt,order)),inv),Minus(infor234(T,x1,order),infor234(T,xl,order))),inv),Minus(infor234(T,x2,order),infor234(T,xl,order))),inv),g2(T,dt,order))(0,0);
    double s3=Multiple(Multiple(Trans(g2(T,dt,order)),inv),g2(T,dt,order))(0,0);
    double s4=Multiple(Multiple(Multiple(Multiple(Trans(g2(T,dt,order)),inv),Minus(infor234(T,x1,order),infor234(T,xl,order))),inv),g2(T,dt,order))(0,0);
    double s5=Multiple(Multiple(Multiple(Multiple(Trans(g2(T,dt,order)),inv),Minus(infor234(T,x2,order),infor234(T,xl,order))),inv),g2(T,dt,order))(0,0);
    double s6=pow(Multiple(Multiple(Trans(g2(T,dt,order)),inv),g2(T,dt,order))(0,0),-2);

    return (s1+s2*s3-s4*s5)*s6;
}

// [[Rcpp::export]]

double ds1(NumericVector T, double x, NumericMatrix inv, int order)
{
    return Multiple(Multiple(Trans(f234(T,x,order)),inv),f234(T,x,order))(0,0)/order;
}

// [[Rcpp::export]]

double ds2(NumericVector T, double x, NumericMatrix inv, int order)
{
    double s1=Multiple(Multiple(Trans(f234(T,x,order)),inv),g1(T,order))(0,0);
    double s2=Multiple(Multiple(Trans(g1(T,order)),inv),f234(T,x,order))(0,0);
    double s3=pow(Multiple(Multiple(Trans(g1(T,order)),inv),g1(T,order))(0,0),-1);

    return s1*s2*s3;
}

// [[Rcpp::export]]

double ds3(NumericVector T, double x, NumericMatrix inv, double dt, int order)
{
    double s1=Multiple(Multiple(Trans(f234(T,x,order)),inv),g2(T,dt,order))(0,0);
    double s2=Multiple(Multiple(Trans(g2(T,dt,order)),inv),f234(T,x,order))(0,0);
    double s3=pow(Multiple(Multiple(Trans(g2(T,dt,order)),inv),g2(T,dt,order))(0,0),-1);

    return s1*s2*s3;
}

// [[Rcpp::export]]

NumericMatrix sMultiple(double s, NumericMatrix M)
{
    int nrow=M.nrow(),ncol=M.ncol();
    NumericMatrix N(nrow,ncol);

    for(int i=0;i<nrow;i++)
        for(int j=0;j<ncol;j++)
            N(i,j)=s*M(i,j);

    return N;
}

// [[Rcpp::export]]

NumericMatrix upinfor(NumericVector W, NumericVector T, NumericVector X, int order)
{
    NumericMatrix last_infor=infor234(T,X[X.size()-1],order);
    NumericMatrix infor=sMultiple(1-sum(W),last_infor);

    for(int i=0;i<X.size()-1;i++)
        infor=Plus(infor,sMultiple(W[i],infor234(T,X[i],order)));

    return infor;
}


// [[Rcpp::export]]

NumericVector c1_weight_1(NumericVector W, NumericVector T, NumericVector X, NumericMatrix inv, int order)
{
    int p=W.size(),k=X.size();
    NumericMatrix V=Multiple(g1(T,order),Trans(g1(T,order)));
    NumericVector f1(p);

    for(int i=0;i<p;i++)
        f1[i]=-SDM(Multiple(Multiple(Multiple(inv,Minus(infor234(T,X[i],order),infor234(T,X[k-1],order))),inv),V));

    return f1;
}

// [[Rcpp::export]]

NumericMatrix c1_weight_2(NumericVector W, NumericVector T, NumericVector X, NumericMatrix inv, int order)
{
    int p=W.size(),k=X.size();
    NumericMatrix V=Multiple(g1(T,order),Trans(g1(T,order)));
    NumericMatrix f2(p,p);
    NumericMatrix M1,M2;

    for(int i=0;i<p;i++)
        for(int j=0;j<p;j++) {
            M1=Multiple(Multiple(Multiple(Multiple(inv,Minus(infor234(T,X[j],order),infor234(T,X[k-1],order))),inv),Minus(infor234(T,X[i],order),infor234(T,X[k-1],order))),inv);
            M2=Multiple(Multiple(Multiple(Multiple(inv,Minus(infor234(T,X[i],order),infor234(T,X[k-1],order))),inv),Minus(infor234(T,X[j],order),infor234(T,X[k-1],order))),inv);
            f2(i,j)=SDM(Multiple(Plus(M1,M2),V));
        }

    return f2;
}

// [[Rcpp::export]]

NumericVector c_weight_1(NumericVector W, NumericVector T, NumericVector X, NumericMatrix inv, double dt, int order)
{
    int p=W.size(),k=X.size();
    NumericMatrix V=Multiple(g2(T,dt,order),Trans(g2(T,dt,order)));
    NumericVector f1(p);

    for(int i=0;i<p;i++)
        f1[i]=d3(T,X[i],X[k-1],inv,dt,order);

    return f1;
}

// [[Rcpp::export]]

NumericMatrix c_weight_2(NumericVector W, NumericVector T, NumericVector X, NumericMatrix inv, double dt, int order)
{
    int p=W.size(),k=X.size();
    NumericMatrix V=Multiple(g2(T,dt,order),Trans(g2(T,dt,order)));
    NumericMatrix f2(p,p);

    for(int i=0;i<p;i++)
        for(int j=0;j<p;j++)
            f2(i,j)=dd3(T,X[i],X[j],X[k-1],inv,dt,order);

    return f2;
}

// [[Rcpp::export]]

NumericVector D_weight_1(NumericVector W, NumericVector T, NumericVector X, NumericMatrix inv, int order)
{
    int p=W.size(),k=X.size();
    NumericVector f1(p);

    for(int i=0;i<p;i++)
        f1[i]=SDM(Multiple(inv,Minus(infor234(T,X[i],order),infor234(T,X[k-1],order))));

    return f1;
}

// [[Rcpp::export]]

NumericMatrix D_weight_2(NumericVector W, NumericVector T, NumericVector X, NumericMatrix inv, int order)
{
    int p=W.size(),k=X.size();
    NumericMatrix f2(p,p);

    for(int i=0;i<p;i++)
        for(int j=0;j<p;j++)
            f2(i,j)=-SDM(Multiple(Multiple(Multiple(inv,Minus(infor234(T,X[j],order),infor234(T,X[k-1],order))),inv),Minus(infor234(T,X[i],order),infor234(T,X[k-1],order))));

    return f2;
}

// [[Rcpp::export]]

NumericVector M_weight_1(NumericVector W, NumericVector T, NumericVector X, NumericMatrix inv, double dt, int order, NumericVector lambda)
{
    int p=W.size(),k=X.size();
    NumericVector f1(p);

    for(int i=0;i<p;i++)
        f1[i]=lambda[0]*d1(T,X[i],X[k-1],inv,order)/order-lambda[1]*d2(T,X[i],X[k-1],inv,order)-(1-lambda[0]-lambda[1])*d3(T,X[i],X[k-1],inv,dt,order);

    return f1;
}

// [[Rcpp::export]]

NumericMatrix M_weight_2(NumericVector W, NumericVector T, NumericVector X, NumericMatrix inv, double dt, int order, NumericVector lambda)
{
    int p=W.size(),k=X.size();
    NumericMatrix f2(p,p);

    for(int i=0;i<p;i++)
        for(int j=0;j<p;j++)
            f2(i,j)=lambda[0]*dd1(T,X[i],X[j],X[k-1],inv,order)/order-lambda[1]*dd2(T,X[i],X[j],X[k-1],inv,order)-(1-lambda[0]-lambda[1])*dd3(T,X[i],X[j],X[k-1],inv,dt,order);

    return f2;
}



