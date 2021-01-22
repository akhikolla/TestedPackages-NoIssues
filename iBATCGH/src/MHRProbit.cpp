// [[Rcpp::depends(iBATCGH]]
#include <MHRProbit.h>

// [[Rcpp::interfaces(cpp)]]
double MHR2(double phi,int* mR,double* mX,double* mY,double c,double d,double delta,double* mSum,arma::icolvec *Rout,double alpha0,double alpha1,int* mXI,arma::mat Hs, arma::vec countnotwo, int nR, int choiceg, int selectioncgh,int m,int mg, int g, int nOss)
{
arma::icolvec R(mR,mg,false);
arma::icolvec Rnew(mg);
Rnew=R;
int* mRnew=Rnew.memptr();
int acc=0;
double piRxi=0;
int chg;
int choice;
int choice1;
int choice2;
arma::ivec flag(mg);
int* mFlag=flag.memptr(); 
    
//Choose between Adding/Deleting or Swap step 
	GetRNGstate();
	double x= unif_rand();
	PutRNGstate();
chg=choiceg*m;

//zeros are the zero position in R, ones are the one position in R, countz number of zeros and counto number of ones
arma::uvec ones = find(R.rows(chg,(chg+m-1)));
arma::uvec zeros = find(R.rows(chg,(chg+m-1))==0);
int countz=zeros.n_rows;
int counto=ones.n_rows;

//zerosok are the zero position in R that can be candidate for being changed into a one position
arma::uvec zerosok(countz);
zerosok.zeros();
int countzok=0;
for(int j=0; j<countz; j++){
        if(countnotwo[zeros[j]]>selectioncgh){
                     zerosok[countzok]=zeros[j];
                     countzok++;
        }
}
//Put the positions that could be changed into a single vector
arma::uvec posok = arma::join_cols(ones,zerosok);

if(ones.n_rows!=0 || countzok!=0){

//AD step
if(countz==m || countz==0 || x<phi || countzok==0){
          //Select position to be changed
		  GetRNGstate();
          choice=posok[floor(unif_rand()*(countzok+counto))];
	      PutRNGstate();
          //Change value of the position selected
          Rnew[chg+choice]=abs(Rnew[chg+choice]-1);
          //Compute log\pi(Rnew|...)-log\pi(Rold|...)
          piRxi=piRxi+pRAD2(choice,chg,m,alpha0,alpha1,mSum,mRnew,mR,mg);
          //Add the log ratio of the likelihoods: log(likelihood|Rnew)- log(likelihood|Rold)      
          piRxi=piRxi+mml(mRnew,mXI,mY,c,d,delta,chg+choice,Hs,nOss,m,g,mg)-mml(mR,mXI,mY,c,d,delta,chg+choice,Hs,nOss,m,g,mg);
}

//Swap step
else{
     //AD[i]=0;
     //Select positions to be changed
	 GetRNGstate();
     choice1=zerosok[floor(unif_rand()*(countzok))];
     choice2=ones[floor(unif_rand()*(counto))];
	 PutRNGstate();
     //Change values of the positions selected
     Rnew[chg+choice1]=abs(Rnew[chg+choice1]-1);
     Rnew[chg+choice2]=abs(Rnew[chg+choice2]-1);
     flag.zeros();
     //Compute log(\pi(Rnew|...))-log(\pi(Rold|...))
     piRxi=piRxi+pRS2(choice1,chg,m,alpha0,alpha1,mSum,mRnew,mR,mFlag,&flag,mg);
     piRxi=piRxi+pRS2(choice2,chg,m,alpha0,alpha1,mSum,mRnew,mR,mFlag,&flag,mg);
     //Add the log ratio of the likelihoods: log(likelihood|Rnew)- log(likelihood|Rold)
     piRxi=piRxi+mml(mRnew,mXI,mY,c,d,delta,chg+choice1,Hs,nOss,m,g,mg)-mml(mR,mXI,mY,c,d,delta,chg+choice1,Hs,nOss,m,g,mg);
}  
}

//Accept/Reject the proposed value/values
GetRNGstate();
double random=unif_rand(); 
PutRNGstate();
double move = (0<piRxi)?0:piRxi;
if(random<exp(move)||move==0){
                              *Rout=Rnew;
                              acc=1;}

return acc;//Return one if accepted, zero if not (acceptance rate)

}

