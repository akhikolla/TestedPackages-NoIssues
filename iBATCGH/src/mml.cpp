#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(cpp)]]
double mml(int* mR,int* mX,double* mY, double c, double d, double delta, int elem,arma::mat Hs,int nOss,int p,int g,int mg) 
{
arma::icolvec R(mR,mg,false);
arma::mat Y(mY,nOss,g,false);
arma::imat X(mX,nOss,p,false);

int ind=elem/p;
double output;

//log-marginal model likelihood of Yi

//step one: Extract the elements of R involved (g-th regression) and compute the number of ones (included regressors)
arma::icolvec RYi=R.subvec(ind*p,(ind+1)*p-1);
int ki=sum(RYi);

//step two: Extract g-th column of Y
arma::colvec Yi=Y.col(ind);

//step three: extract the included columns of X
arma::mat Zi(nOss,ki);
int temp = 0;
arma::mat Ui(ki,ki);
double qi;
//If it is not the null model
if (ki>0){

   for(int i = 0 ; i < p ; i++){
           if(RYi[i]==1){
                  for(int j = 0 ; j < nOss ; j++){
                          Zi(j,temp)=X(j,i);       
                  }
           temp++;
           }    
   }    

//step four: Compute U_g
Ui=(arma::trans(Zi)*Hs*Zi);
for (int i=0; i<ki; i++){
     Ui(i,i) += c;
}

//step five: Compute q_g
qi=arma::as_scalar((arma::trans(Yi)*Hs*Yi)-(arma::trans(Yi)*Hs*Zi*arma::inv(Ui)*arma::trans(Zi)*Hs*Yi));

//step six: Compute likelihood
output=ki*0.5*log(c)-0.5*(log(arma::det(Ui)))-((nOss+delta)*0.5)*log(d+qi);
}
//If there are no regressor included
else{

//step four: Compute U_g=0
//step five: Compute q_g
double qi=arma::as_scalar((arma::trans(Yi)*Hs*Yi));
//step six: Compute likelihood
output=-((nOss+delta)*0.5)*log(d+qi);
}

return output;
}

