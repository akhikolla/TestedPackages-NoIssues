// [[Rcpp::depends(iBATCGH]]
#include <MHRmixture.h>

// [[Rcpp::interfaces(cpp)]]
double MHR(double phi,int* mR,double* mX,double* mY,double c,double d,double delta,double* mGam,double* mOm1,double* mOm2,double e,double f,arma::icolvec *Rout,double add0,double add1,int* mXI,arma::mat Hs,arma::vec countnotwo,int nR,int choiceg,int selectioncgh,int m,int mg, int g, int nOss, int indep)
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
//Rcpp::Rcout << "Line 2 done."<< std::endl;    
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
          //Change the value of the selected position
          Rnew[chg+choice]=abs(Rnew[chg+choice]-1);
          //Compute log\pi(Rnew|...)-log\pi(Rold|...)
	if(indep != 0){
          piRxi=piRxi+pRAD(choice,chg,m,add0,add1,mGam,mOm1,mOm2,mRnew,mR,mg);
	}
	if(indep==0){
		piRxi=piRxi+log(((Rnew[chg+choice]==1)?add1:add0))-log(((Rnew[chg+choice]==1)?add0:add1));
	}
          //Add the log ratio of the likelihoods: log(likelihood|Rnew)- log(likelihood|Rold)      
          piRxi=piRxi+mml(mRnew,mXI,mY,c,d,delta,chg+choice,Hs,nOss,m,g,mg)-mml(mR,mXI,mY,c,d,delta,chg+choice,Hs,nOss,m,g,mg);
}

//Swap step
else{
     //Select positions to be changed
	 GetRNGstate();
     choice1=zerosok[floor(unif_rand()*(countzok))];
     choice2=ones[floor(unif_rand()*(counto))];
	 PutRNGstate();
     //Change the value of the selected positions
     Rnew[chg+choice1]=abs(Rnew[chg+choice1]-1);
     Rnew[chg+choice2]=abs(Rnew[chg+choice2]-1);
     flag.zeros();
     //Compute log\pi(Rnew|...)-log\pi(Rold|...)
	if(indep != 0){
		piRxi=piRxi+pRS(choice1,chg,m,add0,add1,mGam,mOm1,mOm2,mRnew,mR,mFlag,&flag,mg);
		piRxi=piRxi+pRS(choice2,chg,m,add0,add1,mGam,mOm1,mOm2,mRnew,mR,mFlag,&flag,mg);
	}
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

return acc;//return one if accepted zero if not (acceptance rate)

}

