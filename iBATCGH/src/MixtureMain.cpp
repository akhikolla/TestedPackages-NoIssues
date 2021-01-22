// [[Rcpp::depends(iBATCGH]]
#include <MixtureMain.h>

// [[Rcpp::export]]
Rcpp::List iBAT(arma::mat Y,arma::mat X,arma::colvec distance,double disfix,int intercept,arma::imat xi,arma::icolvec R,arma::mat tran,arma::colvec mu,arma::colvec sigma,double cmu,double c,double delta,double d,double e,double f,double alpha,arma::colvec deltak,arma::colvec tauk,arma::colvec upp_bounds,arma::colvec low_bounds,arma::colvec alpha_IG,arma::colvec beta_IG,arma::colvec low_IG,arma::rowvec a,int niter,int burnin,int Cout,double phi,float pR,int selectioncgh,float pXI, int indep){
//Declare and initialize variables
double den=(exp(1.0)-1.0);//Fixed quantity
double betaden = Rf_beta(e,f);//Fixed quantity

int ns=tran.n_rows;
arma::mat A(ns,ns);
arma::vec Pi(ns);
arma::vec Picum(ns);
arma::mat tranew(ns,ns);
arma::vec Pinew(ns);
arma::uvec xieq;

//Armadillo pointers
double* mX;
double* mY;
double* mGam;
double* mOm1;
double* mOm2;
double* mSum;
double* mDist;
int* mXI;
int* mR;
mR=R.memptr();
mX=X.memptr();
mY=Y.memptr();
mXI=xi.memptr();
mDist=distance.memptr();

int accR=0;
int accXI=0;
int accA=0;

Pi=initial(tran);//Probability for the fisrt CGH
A=arma::cumsum(tran,1);
Picum=arma::cumsum(Pi);

Rcpp::List output(4*niter+6);
int m=xi.n_cols;
int nOss=xi.n_rows;
int g=Y.n_cols;
int mg=m*g;
if(R.n_rows==1){R.resize(mg);R.zeros();mR=R.memptr();}
if(selectioncgh==-1){selectioncgh=0.1*nOss;}

arma::imat xi1(nOss,m);
arma::imat xi3(nOss,m);
arma::imat xi4(nOss,m);
xi1.zeros();
xi3.zeros();
xi4.zeros();
arma::colvec gamma(m);
arma::colvec omega1(m);
arma::colvec omega2(m);
arma::colvec sum(m+1);
mGam=gamma.memptr();
mOm1=omega1.memptr();
mOm2=omega2.memptr();
mSum=sum.memptr();
double piA;
arma::uvec xipos1;
arma::uvec xipos2;
arma::uvec xipos3;
arma::uvec xipos4;
arma::vec n_k(4);
arma::vec theta2(4);
arma::vec weight_means(4);
arma::vec V_k(4);


//Compute H_n
arma::mat os(nOss,nOss);
arma::mat Is(nOss,nOss);
arma::mat Hs(nOss,nOss);
Is.eye();
os.ones();

if(intercept==1){Hs=Is - os/(nOss+cmu);}
if(intercept!=1){Hs=Is;}

//Counts the number of transition given estimated xi
arma::mat counts(4,4);
counts.zeros();
for (int i=0; i<nOss; i++){
    for(int j=0; j<(m-1);j++){
            counts(xi(i,j)-1,xi(i,j+1)-1)++;
    }
}

//Find the number of rows not in neutral state for each column
arma::vec countnotwo(m);

for(int i=0; i<m; i++){
arma::uvec notwo = find(xi.col(i)!=2);
countnotwo[i]=notwo.n_rows;
}

double add0 = Rf_beta(e,f+1)/betaden;//Fixed value
double add1 = Rf_beta(e+1,f)/betaden;//Fixed value

//Compute initial values of gamma, omega1, omega2 e sum
omegagamma(distance,xi,disfix,alpha,&gamma,&omega1,&omega2,&sum);
arma::vec vetac(niter);

for (int i=0; i<niter; i++){
if (i%Cout==0 && i!=0){Rcpp::Rcout << "Iteration " << i << " done."<< std::endl;}//Print every 1,000 iter
//MH for the association matrix R
//Extract genes to be changed
arma::icolvec choiceg(g);
choiceg.zeros();
int nR=0;
int aux=0;
int sommo=0;
do{
GetRNGstate();
aux=(int) Rf_rgeom(pR);
PutRNGstate();
sommo=sommo+aux;
choiceg[nR]=sommo;
sommo++;
nR=nR+1;
}while(sommo<g+1);
nR=nR-1;

if(nR<1){
aux=floor(unif_rand()*g);
choiceg[0]=aux;
nR=1;          
}
choiceg.resize(nR);

arma::vec accRvec(nR);
accRvec.zeros();
//For each gene selected perform MH step
for (int ggg=0; ggg<nR; ggg++){
accRvec[ggg]=MHR(phi,mR,mX,mY,c,d,delta,mGam,mOm1,mOm2,e,f,&R,add0,add1,mXI,Hs,countnotwo,nR,choiceg[ggg],selectioncgh,m,mg,g,nOss,indep);
}
arma::uvec nonesR=find(accRvec==1);
if(nonesR.n_rows>0){accR=accR+1;}

// MH for the matrix of latent states xi
//Select samples to be changed
arma::icolvec sel(nOss);
sel.zeros();
int nxi=0;
aux=0;
sommo=0;
do{
GetRNGstate();
aux=(int) Rf_rgeom(pXI);
PutRNGstate();
sommo=sommo+aux;
sel[nxi]=sommo;
sommo++;
nxi=nxi+1;
}while(sommo<nOss+1);
nxi=nxi-1;
if(nxi<1){
aux=floor(unif_rand()*nOss);
sel[0]=aux;
nxi=1;         
}
sel.resize(nxi);

arma::vec accXIvec(nxi);
accXIvec.zeros();
int choice = floor(unif_rand()*m);
//For each sample selected perform MH step
for(int sss=0; sss<nxi; sss++){
if(choice!=(m-1)){counts(xi(sel[sss],choice)-1,xi(sel[sss],choice+1)-1) -= 1;}
if(xi(sel[sss],choice)==2){countnotwo[choice] -=1;}
accXIvec[sss]=MHxi(mXI,A,Picum,mDist,disfix,alpha,mGam,mOm1,mOm2,mSum,e,f,mR,mY,c,d,delta,mu,sigma,mX,tran,Pi,&xi,&gamma,&omega1,&omega2,&sum,den,add0,add1,Hs,choice,sel[sss],m,mg,g,nOss,indep);
if(choice!=(m-1)){counts(xi(sel[sss],choice)-1,xi(sel[sss],choice+1)-1) += 1;}
if(xi(sel[sss],choice)==2){countnotwo[choice] +=1;}
}
arma::uvec nonesXI=find(accXIvec==1);
if(nonesXI.n_rows>0){accXI=accXI+1;}

//Metropolis for the transition matrix
tranew.row(0)=dir((a+counts.row(0)));
tranew.row(1)=dir((a+counts.row(1)));
tranew.row(2)=dir((a+counts.row(2)));
tranew.row(3)=dir((a+counts.row(3)));
Pinew=initial(tranew);
piA=0;
xieq=find(xi.col(0)==1);
piA=piA+xieq.n_rows*(log(Pinew[0])-log(Pi[0]));
xieq=find(xi.col(0)==2);
piA=piA+xieq.n_rows*(log(Pinew[1])-log(Pi[1]));
xieq=find(xi.col(0)==3);
piA=piA+xieq.n_rows*(log(Pinew[2])-log(Pi[2]));
xieq=find(xi.col(0)==4);
piA=piA+xieq.n_rows*(log(Pinew[3])-log(Pi[3]));
GetRNGstate();
double random=unif_rand();
PutRNGstate();
double move = (0<piA)?0:piA;

if(random<exp(move)||move==0){
	tran=tranew;
	Pi=Pinew;
    A=arma::cumsum(tran,1);
    Picum=arma::cumsum(Pi);
    accA++;
}

//Store the results at each iteration
arma::uvec Rpos=find(R==1);
output[i]=Rpos;
xipos1=find(xi==1);
if(i>=burnin){xi1.elem(xipos1)=xi1.elem(xipos1)+1;}
xipos3=find(xi==3);
if(i>=burnin){xi3.elem(xipos3)=xi3.elem(xipos3)+1;}
xipos4=find(xi==4);
if(i>=burnin){xi4.elem(xipos4)=xi4.elem(xipos4)+1;}
output[i+niter]=tran;

// Gibbs for the state specific mean and variance
xipos2=find(xi==2);
n_k[0]=xipos1.n_rows;
n_k[1]=xipos2.n_rows;
n_k[2]=xipos3.n_rows;
n_k[3]=xipos4.n_rows;


for (int j=0; j<4; j++){
    if(j==3){low_bounds[3]=mu[2]+3*sigma[2];}
    if(n_k[j]==0){
				  GetRNGstate();
                  mu[j]=truncnorm(deltak[j],tauk[j],low_bounds[j],upp_bounds[j]);
                  sigma[j]=pow((truncgamma(alpha_IG[j],1/(beta_IG[j]),pow(low_IG[j],-2))),-0.5);
                  PutRNGstate();
    }
    if(n_k[j]!=0){
                  theta2[j]=pow(tauk[j],-2)+(n_k[j]*pow(sigma[j],-2));
                  if (j==0){weight_means[j]=( deltak[j]*pow(tauk[j],-2) + mean(X.elem(xipos1))* (n_k[j]*pow(sigma[j],-2))  )/theta2[j];}
                  if (j==1){weight_means[j]=( deltak[j]*pow(tauk[j],-2) + mean(X.elem(xipos2))* (n_k[j]*pow(sigma[j],-2))  )/theta2[j];}
                  if (j==2){weight_means[j]=( deltak[j]*pow(tauk[j],-2) + mean(X.elem(xipos3))* (n_k[j]*pow(sigma[j],-2))  )/theta2[j];}
                  if (j==3){weight_means[j]=( deltak[j]*pow(tauk[j],-2) + mean(X.elem(xipos4))* (n_k[j]*pow(sigma[j],-2))  )/theta2[j];}  
                  GetRNGstate();
				  mu[j]=truncnorm(weight_means[j],1/sqrt(theta2[j]),low_bounds[j],upp_bounds[j]);
                  if (j==0){V_k[j]=arma::as_scalar(arma::trans(X.elem(xipos1)-mu[j])*(X.elem(xipos1)-mu[j]));}
                  if (j==1){V_k[j]=arma::as_scalar(arma::trans(X.elem(xipos2)-mu[j])*(X.elem(xipos2)-mu[j]));}
                  if (j==2){V_k[j]=arma::as_scalar(arma::trans(X.elem(xipos3)-mu[j])*(X.elem(xipos3)-mu[j]));}
                  if (j==3){V_k[j]=arma::as_scalar(arma::trans(X.elem(xipos4)-mu[j])*(X.elem(xipos4)-mu[j]));}
                  sigma[j]=pow((truncgamma(alpha_IG[j]+n_k[j]/2,1/(beta_IG[j]+V_k[j]/2),pow(low_IG[j],-2))),-0.5);
                  PutRNGstate();
    }
}
//Store mu and sigma
output[i+2*niter]=mu;
output[i+3*niter]=sigma;
}
//Store xi and acceptance rates
output[4*niter]=xi1;
output[4*niter+1]=xi3;
output[4*niter+2]=xi4;
output[4*niter+3]=accR;
output[4*niter+4]=accXI;
output[4*niter+5]=accA;
return output;
}
