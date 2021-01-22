// [[Rcpp::depends(iBATCGH]]
#include <MHXImixture.h>

// [[Rcpp::interfaces(cpp)]]
double MHxi(int* mXI,arma::mat A,arma::mat Picum,double* mDist,double disfix,double alpha,double* mGam,double* mOm1,double* mOm2,double* mSum,double e,double f,int* mR,double* mY,double c,double d,double delta,arma::colvec mu,arma::colvec sigma,double* mX,arma::mat tran,arma::colvec Pi,arma::imat *xiout,arma::colvec *gammaout,arma::colvec *omega1out,arma::colvec *omega2out, arma::colvec *sumout,double den,double add0,double add1,arma::mat Hs,int choice,int sel,int m,int mg,int g,int nOss, int indep)
{

arma::icolvec R(mR,mg,false);
arma::mat X(mX,nOss,m,false);
arma::colvec gammaold(mGam,m,false);
arma::colvec omega1old(mOm1,m,false);
arma::colvec omega2old(mOm2,m,false);
arma::colvec sumold(mSum,(m+1),false);
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
arma::colvec gamma(m);
arma::colvec omega1(m);
arma::colvec omega2(m);
arma::colvec sum(m+1);
sum=sumold;
gamma=gammaold;
omega1=omega1old;
omega2=omega2old;
double* mGamNew=gamma.memptr();
double* mOm1New=omega1.memptr();
double* mOm2New=omega2.memptr();

if(choice>0){
              arma::uvec stemp=find(xinew.col(choice)==xinew.col(choice-1));
              sum[choice]=stemp.n_rows*((exp(1-(distance[choice-1]/disfix))-1)/DEN)/nOss;
}

if(choice<m-1){
               arma::uvec stemp=find(xinew.col(choice)==xinew.col(choice+1));
               sum[choice+1]=stemp.n_rows*((exp(1-(distance[choice]/disfix))-1)/DEN)/nOss;
}

if(choice==0){
gamma[choice]=alpha/(alpha+sum[choice]+sum[choice+1]);
gamma[choice+1]=alpha/(alpha+sum[choice+1]+sum[choice+2]);
omega1[choice+1]=sum[choice+1]/(alpha+sum[choice+1]+sum[choice+2]);
omega2[choice]=sum[choice+1]/(alpha+sum[choice]+sum[choice+1]);
omega2[choice+1]=sum[choice+2]/(alpha+sum[choice+1]+sum[choice+2]);              
}

if(choice==m-1){
gamma[choice]=alpha/(alpha+sum[choice]+sum[choice+1]);
gamma[choice-1]=alpha/(alpha+sum[choice-1]+sum[choice]);
omega1[choice]=sum[choice]/(alpha+sum[choice]+sum[choice+1]);
omega1[choice-1]=sum[choice-1]/(alpha+sum[choice-1]+sum[choice]);
omega2[choice-1]=sum[choice]/(alpha+sum[choice-1]+sum[choice]);              
}

if (choice>0 && choice<m-1){
gamma[choice]=alpha/(alpha+sum[choice]+sum[choice+1]);
gamma[choice-1]=alpha/(alpha+sum[choice-1]+sum[choice]);
gamma[choice+1]=alpha/(alpha+sum[choice+1]+sum[choice+2]);
omega1[choice]=sum[choice]/(alpha+sum[choice]+sum[choice+1]);
omega1[choice-1]=sum[choice-1]/(alpha+sum[choice-1]+sum[choice]);
omega1[choice+1]=sum[choice+1]/(alpha+sum[choice+1]+sum[choice+2]);
omega2[choice]=sum[choice+1]/(alpha+sum[choice]+sum[choice+1]);
omega2[choice-1]=sum[choice]/(alpha+sum[choice-1]+sum[choice]);
omega2[choice+1]=sum[choice+2]/(alpha+sum[choice+1]+sum[choice+2]);
}

//Compute log\pi(R|xinew)-log\pi(R|xiold) and log\pi(Y|xinew)-log\pi(Y|xiold)

double piRxi=0;
double piYxi=0;
int jm;

for (int j = 0; j < g; j++){
    jm=j*m;
    //Compute log(p(R|xinew))-log(p(R|xi))
	if(indep != 0){
		piRxi=piRxi+pRXI(choice,jm,m,add0,add1,mGam,mOm1,mOm2,mGamNew,mOm1New,mOm2New,mR,mg);
    }
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
	*gammaout=gamma;
	*omega1out=omega1;
	*omega2out=omega2;
	*sumout=sum;
	acc=1;
}

return acc;
}
else{
return 0;
}
}

