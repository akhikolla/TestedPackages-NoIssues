// [[Rcpp::depends(iBATCGH]]
#include <MHXIProbit.h>

// [[Rcpp::interfaces(cpp)]]
double MHxi2(int* mXI,arma::mat A,arma::mat Picum,double* mDist,double disfix,double* mSum,int* mR,double* mY,double c,double d,double delta,arma::colvec mu,arma::colvec sigma,double* mX,arma::mat tran,arma::colvec Pi,arma::imat *xiout,arma::colvec *sumout,double den,double alpha0,double alpha1,arma::mat Hs,int choice,int sel,int m,int mg,int g,int nOss)
{
arma::icolvec R(mR,mg,false);
arma::colvec sumold(mSum,(m+1),false);
arma::mat X(mX,nOss,m,false);
arma::imat xi(mXI,nOss,m,false);
arma::colvec distance(mDist,m,false);
arma::imat xinew(nOss,m);
xinew=xi;
int* mXInew=xinew.memptr();

//Proposed xi
double ran;

if (choice==0){
					GetRNGstate();
                    ran= unif_rand();
					PutRNGstate();
                    if(ran<Picum[0]){xinew(sel,0)=1;}
                    if(Picum[0]<=ran && ran<Picum[1]){xinew(sel,0)=2;}
                    if(Picum[1]<=ran && ran<Picum[2]){xinew(sel,0)=3;}
                    if(ran>Picum[2]){xinew(sel,0)=4;}
}
if(choice!=0){
                 GetRNGstate();
				 ran= unif_rand();
				 PutRNGstate();
                 if(ran<A(xi(sel,choice-1)-1,0)){xinew(sel,choice)=1;}
                 if(A(xi(sel,choice-1)-1,0)<=ran && ran<A(xi(sel,choice-1)-1,1)){xinew(sel,choice)=2;}
                 if(A(xi(sel,choice-1)-1,1)<=ran && ran<A(xi(sel,choice-1)-1,2)){xinew(sel,choice)=3;}
                 if(ran>A(xi(sel,choice-1)-1,2)){xinew(sel,choice)=4;}
}
if(xinew(sel,choice)!=xi(sel,choice)){
//Compute new gamma, omega1, omega2
arma::colvec sum(m+1);
sum=sumold;
double* mSumNew=sum.memptr();

if(choice>0){
              arma::uvec stemp=find(xinew.col(choice)==xinew.col(choice-1));
              sum[choice]=stemp.n_rows*((exp(1-(distance[choice-1]/disfix))-1)/DEN)/nOss;
}

if(choice<m-1){
               arma::uvec stemp=find(xinew.col(choice)==xinew.col(choice+1));
               sum[choice+1]=stemp.n_rows*((exp(1-(distance[choice]/disfix))-1)/DEN)/nOss;
}

//Compute log\pi(R|xinew)-log\pi(R|xiold) and log\pi(Y|xinew)-log\pi(Y|xiold)

double piRxi=0;
double piYxi=0;
int jm;

for (int j = 0; j < g; j++){
    jm=j*m;
    //Compute log(p(R|xinew))-log(p(R|xi))
    piRxi=piRxi+pRXI2(choice,jm,m,alpha0,alpha1,mSum,mSumNew,mR,mg);
    //Compute log\pi(Y|xinew)-log\pi(Y|xiold)    
    if(R[jm+choice]==1){piYxi=piYxi+mml(mR,mXInew,mY,c,d,delta,jm+choice,Hs,nOss,m,g,mg)-mml(mR,mXI,mY,c,d,delta,jm+choice,Hs,nOss,m,g,mg);}
}

//Compute log\pi(X|xinew)-log\pi(X|xiold) and log\pi(xinew|xiold)-log\pi(xiold|xinew)

double piXxi=0;
double pixi=0;

//Compute log\pi(X|xinew)-log\pi(X|xiold)
piXxi=piXxi+ldnwpi(X(sel,choice),mu[xinew(sel,choice)-1],sigma[xinew(sel,choice)-1])-ldnwpi(X(sel,choice),mu[xi(sel,choice)-1],sigma[xi(sel,choice)-1]);
//Compute log\pi(xinew|xiold)-log\pi(xiold|xinew)
if(choice<(m-1)){pixi=pixi+log(tran(xinew(sel,choice)-1,xinew(sel,(choice+1))-1))-log(tran(xi(sel,choice)-1,xi(sel,(choice+1))-1));}

//Accept/Reject the proposed value
double pi2=piRxi+piYxi+piXxi+pixi;
GetRNGstate();
double random=unif_rand();
PutRNGstate();
double move = (0<pi2)?0:pi2;
int acc=0;
if(random<exp(move)||move==0){
	*xiout=xinew;
	*sumout=sum;
	acc=1;
}

return acc;
}
else{
return 0;     
}
}

