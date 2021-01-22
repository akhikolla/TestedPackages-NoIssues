#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]

double Qfunc(double x, int qq=999) {
  using namespace arma;
  using namespace std;
  double Qf=1;
  double pi = M_PI;
  
  switch (qq)
  {
  case 1:
    Qf = 1;
    break;
  case 2:
    Qf =1.0/3.0 + 2*pow(cos(pi*x),2.0)/3.0;
    break;
  case 3:
    Qf = 2.0/15.0 + 11*pow(cos(pi*x),2.0)/15.0 + 2*pow(cos(pi*x),4.0)/15.0;
    break;
  case 4:
    Qf =(17+180*pow(cos(pi*x),2)+114*pow(cos(pi*x),4)+4*pow(cos(pi*x),6))/315.0;
    break;
  case 5:
    Qf =(62+1072*pow(cos(pi*x),2)+1452*pow(cos(pi*x),4)+247*pow(cos(pi*x),6)+2*pow(cos(pi*x),8))/2835.0;
    break;
  case 6:
    Qf =(1382+35396*pow(cos(pi*x),2)+83021*pow(cos(pi*x),4)+34096*pow(cos(pi*x),6)+2026*pow(cos(pi*x),8)+4*pow(cos(pi*x),10))/155295.0;
    break;
  }
  
  return Qf;
}


using namespace Rcpp;
// [[Rcpp::export]]
arma::vec eigenvalues(int nn, int qq=999){
  using namespace arma;
  double pi = M_PI;
  arma::vec svec(nn-qq);
  arma::vec atten(nn-qq);
  arma::vec atten1(2*nn);
  arma::vec ev(nn);
#pragma omp parallel for 
  for (int i=0;i<nn-qq;i++){
    svec(i)=pow(pi*((i+1)+1.0/2*((qq+1)%2)+floor((qq-1)/2)),(2*qq))/nn;}
#pragma omp parallel for 
  for (int j=0;j<2*nn;j++){
    atten1(j)= pow(sin((j+1)*pi/(2*nn))/((j+1)*pi/(2*nn)),2*qq) / Qfunc((j+1)/(2.0*nn),qq);}
#pragma omp parallel for 
  for (int j=0;j<nn-qq;j++){
    atten(j)= atten1(j+qq);}
  svec = svec%atten;
  ev = zeros<vec>(nn);
#pragma omp parallel for 
  for (int i=qq;i<nn;i++)
  {
    ev(i)=svec(i-qq);
  }

  return ev;
  
}


// [[Rcpp::export]]
arma::mat Turn(arma::mat M){
  int nn =arma::rank(M);
  arma::mat Mo,R,CMo;
  arma::qr(Mo,R,M);
  arma::mat A=M+Mo;
  arma::mat B=M-Mo;
  arma::vec a(nn),b(nn);//contain the norm of the coloumns of A and B matrices
  arma::vec x(nn),y(nn); //auxillary vectors
  arma::mat D(nn,2);// same as the R code
  arma::vec d(nn);// same as the R code
  arma::mat dd; // diag(d)
#pragma omp parallel for
  for(int i=0;i<nn;i++){
#pragma omp parallel for
    for(int j=0;j<nn;j++){
      x(j)=A(j,i);
      y(j)=B(i,j);
    }
    a(i) = norm(x,2);
    b(i) = norm(y,2);
  }
  
#pragma omp parallel for
  for(int i=0;i<nn;i++){
    
    D(i,0)=a(i);
    D(i,1)=b(i);
  }
#pragma omp parallel for
  for(int i=0;i<nn;i++){
    if(D(i,0) <= D(i,1)){d(i)=-1;}
    else d(i)=1;
  }
  dd=diagmat(d);
  CMo = Mo*dd;
  return CMo;
}

using namespace Rcpp;
// [[Rcpp::export]]
extern "C" SEXP drbasis(int nn, int qq=999){
  using namespace arma;
  using namespace std;
  Rcpp::List out;
  const double E = exp(1);
  arma::vec t(nn);
  arma::vec k(nn);
  typedef complex<double> dcomp;
  double Pi = M_PI;
  complex<double> comp0 (0,1);
  complex<double> comp1 (-1,0);
  complex<double> comp2 (0,-0.5);
  complex<double> comp3 (1,-1);
  complex<double> comp4 (-1,1);
  complex<double> comp5 (0,0.5);
  complex<double> comp6 (0,0.25);
  int const0 = 1;
  double const1 = -sqrt(3);
  double const2 = 2*sqrt(3);
  double const3 = -sqrt(5);
  double const4 = 6*sqrt(5);
  double const5 = -6*sqrt(5);
  double const6 = -sqrt(7);
  double const7 = 12*sqrt(7);
  double const8 = -30*sqrt(7);
  double const9 = 20*sqrt(7);
  double const10 = 3;
  double const11 = -60;
  double const12 = 270;
  double const13 = -420;
  double const14 = 210;
  double const15 = -sqrt(11);
  double const16 = 30*sqrt(11);
  double const17 = -210*sqrt(11);
  double const18 = 560*sqrt(11);
  double const19 = 630*sqrt(11);
  double const20 = 252*sqrt(11);
  
  
#pragma omp parallel for
  for (int j=0;j<nn;j++){t(j)= j/double((nn-1));}
#pragma omp parallel for
  for (int j=0;j<nn;j++){k(j)= j+1;}

  switch (qq)
  {
    {case 1:
      arma::mat M1(nn,nn);
#pragma omp parallel for
      for (int i=0;i<nn; i++){
#pragma omp parallel for
        for (int j=0;j<nn;j++){
          if (j==0){ M1(i,j)=const0;}
          else {M1(i,j) = sqrt(2)*cos((k(j)-1)*Pi*t(i));}
          M1(i,j) =M1(i,j)/sqrt(nn);
        }
      }
      arma::mat vectors=M1; arma::mat vectorsQR=Turn(M1); arma::vec values=eigenvalues(nn,1);  
      out=Rcpp::List::create(
        Rcpp::Named("eigenvectors") = vectors,
        Rcpp::Named("eigenvectorsQR")   = vectorsQR,
        Rcpp::Named("eigenvalues") = values,
        Rcpp::Named("x") = t);
      break;}
    
      {case 2:
        arma::mat M2(nn,nn);
#pragma omp parallel for
        for (int i=0;i<nn; i++){
#pragma omp parallel for
          for (int j=0;j<nn;j++){
            if (j==0){ M2(i,j)=const0;}
            else if (j==1){M2(i,j)= const1 + const2*t(i);}
            else {dcomp ans = pow(comp1,1+k(j))/pow(E,(-1.5+k(j))*
                                  Pi*(1-t(i))) + pow(E,-((-1.5+k(j))*Pi*t(i))) + sqrt(2)*
                                  cos(Pi/4 +(-1.5 + k(j))*Pi*t(i));
              M2(i,j) =ans.real();}
            M2(i,j) =M2(i,j)/sqrt(nn);
          }
        }
        arma::mat vectors=M2; arma::mat vectorsQR=Turn(M2); arma::vec values=eigenvalues(nn,2);
        out=Rcpp::List::create(
          Rcpp::Named("eigenvectors") = vectors,
          Rcpp::Named("eigenvectorsQR")   = vectorsQR,
          Rcpp::Named("eigenvalues") = values,
          Rcpp::Named("x") = t) ;
        break;}
    
        {case 3:
          arma::mat M3(nn,nn);

#pragma omp parallel for
          for (int i=0;i<nn; i++){
#pragma omp parallel for
            for (int j=0;j<nn;j++){
              if (j==0){ M3(i,j)=1;}
              else if (j==1){M3(i,j)= const1 + const2*t(i);}
              else if (j==2){M3(i,j)= const3 + const4*t(i) + const5* pow(t(i),2);}
              else {dcomp ans = (sqrt(1.5)-comp0/sqrt(2))*
                ((pow(comp1,1+k(j))/pow(E,(pow(comp1,1.0/6))*(-2.0+k(j))*
                Pi*(1-t(i))))+pow(E,-(pow(comp1,1.0/6))*((-2.0+k(j))*Pi*t(i))))+
                (sqrt(1.5)+comp0/sqrt(2))*(pow(comp1,1+k(j))/pow(E,(comp2 +sqrt(3)/2)*
                (-2.0+k(j))*Pi*(1-t(i)))+pow(E,-((comp2 +sqrt(3)/2)*(-2.0+k(j))*
                Pi*t(i))))-sqrt(2)*sin((-2.0 + k(j))*Pi*t(i));
                M3(i,j) =ans.real();}
              M3(i,j) =M3(i,j)/sqrt(nn);
            }
          }
          arma::mat vectors=M3; arma::mat vectorsQR=Turn(M3); arma::vec values=eigenvalues(nn,3);
          out=Rcpp::List::create(
            Rcpp::Named("eigenvectors") = vectors,
            Rcpp::Named("eigenvectorsQR")   = vectorsQR,
            Rcpp::Named("eigenvalues") = values,
            Rcpp::Named("x") = t);
          break;}
    
          {case 4:
            arma::mat M4(nn,nn);
#pragma omp parallel for
            for (int i=0;i<nn; i++){
#pragma omp parallel for
              for (int j=0;j<nn;j++){
                if (j==0){ M4(i,j)=const0;}
                else if (j==1){M4(i,j)= const1 + const2*t(i);}
                else if (j==2){M4(i,j)= const3 + const4*t(i) + const5* pow(t(i),2);}
                else if (j==3){M4(i,j)= const6 + const7*t(i) + const8* pow(t(i),2) + const9* pow(t(i),3);}
                else {dcomp ans =(1 + sqrt(2))*(pow(comp1,1+k(j))/pow(E,(-2.5+k(j))*
                                  Pi*(1-t(i)))+ pow(E,-((-2.5+k(j))*     Pi*t(i))))+(1/sqrt(2) - 
                                  comp0*(1 + 1/sqrt(2)))*(pow(comp1,1+k(j))/pow(E,pow(comp1,0.25)*(-2.5+k(j))*Pi*(1-t(i))) +
                                  pow(E,-(pow(comp1,0.25)*(-2.5+k(j))*Pi*t(i))))+ (1/sqrt(2)+ comp0*(1 + 1/sqrt(2)))*
                                  (pow(comp1,1+k(j))/pow(E,(comp3 *(-2.5 +k(j))*Pi*(1-t(i)))/sqrt(2)) +
                                  pow(E,(comp4 *(-2.5 +k(j))*Pi*t(i))/sqrt(2))) - sqrt(2)*cos(Pi/4 - (-2.5 + k(j))*Pi*t(i));
                  M4(i,j) = ans.real();}
                M4(i,j) = M4(i,j)/sqrt(nn);
              }
            }
            arma::mat vectors=M4; arma::mat vectorsQR=Turn(M4); arma::vec values=eigenvalues(nn,4);
            out=Rcpp::List::create(
              Rcpp::Named("eigenvectors") = vectors,
              Rcpp::Named("eigenvectorsQR")   = vectorsQR,
              Rcpp::Named("eigenvalues") = values,
              Rcpp::Named("x") = t) ;
            break;}
    
            {case 5:
              arma::mat M5(nn,nn);
#pragma omp parallel for
              for (int i=0;i<nn; i++){
#pragma omp parallel for
                for (int j=0;j<nn;j++){
                  if (j==0){ M5(i,j)=const0;}
                  else if (j==1){M5(i,j)= const1 + const2*t(i);}
                  else if (j==2){M5(i,j)= const3 + const4*t(i) + const5* pow(t(i),2);}
                  else if (j==3){M5(i,j)= const6 + const7*t(i) + const8* pow(t(i),2) + 
                    const9* pow(t(i),3);}
                  else if (j==4){M5(i,j)= const10 + const11*t(i) + const12* pow(t(i),2) + 
                    const13* pow(t(i),3) + const14* pow(t(i),4);}
                  else {dcomp ans =(sqrt(2)*(1 + sqrt(5)/2)- (comp5*(sqrt(10 - 2*sqrt(5)) + 
                                    sqrt(2*(5 + sqrt(5)))))/sqrt(2))*(pow(comp1,1 + k(j))/pow(E,pow(comp1,0.1)*
                                    (dcomp(-3 + k(j)))*Pi*(1 - t(i))) + pow(E,-(pow(comp1,0.1)*(dcomp(-3 + k(j)))*
                                    Pi*t(i)))) + (-(1/sqrt(2)) - (comp5*(sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + sqrt(5)))))/
                                      sqrt(2))*(pow(comp1,1 + k(j))/pow(E,pow(comp1,0.3)*(dcomp (-3 + k(j)))*Pi*(1 - t(i)))
                                                  + pow(E,-(pow(comp1,0.3)*(dcomp (-3 + k(j)))*Pi*t(i))))+ (sqrt(2)*(1 + sqrt(5)/2) + 
                                                    (comp5*(sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + sqrt(5)))))/sqrt(2)) *(pow(comp1,1 + k(j))/
                                                      pow(E,(sqrt(0.625 + sqrt(5)/8) - comp6*(-1 + sqrt(5)))*(dcomp (-3 + k(j)))*Pi*
                                                        (1 - t(i))) +pow(E,-((sqrt(0.625 + sqrt(5)/8.) - comp6*(-1 + sqrt(5)))*(dcomp (-3 + k(j)))*
                                                        Pi*t(i)))) +(-(1/sqrt(2)) + (comp5*(sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + sqrt(5)))))/sqrt(2)) 
                    *(pow(comp1,1 + k(j))/pow(E,(dcomp(sqrt(0.625 - sqrt(5)/8)) - comp6*(1 + sqrt(5)))*
                      (dcomp (-3 + k(j)))*Pi*(1 - t(i))) +pow(E,-((sqrt(0.625 - sqrt(5)/8) - comp6*(1 + sqrt(5)))*
                      (dcomp (-3 + k(j)))*Pi*t(i))))- sqrt(2)*cos((-3 + k(j))*Pi*t(i));
                    M5(i,j) =ans.real();}
                  M5(i,j) =M5(i,j)/sqrt(nn);
                }
              }
              arma::mat vectors=M5; arma::mat vectorsQR=Turn(M5); arma::vec values=eigenvalues(nn,5);
              out=Rcpp::List::create(
                Rcpp::Named("eigenvectors") = vectors,
                Rcpp::Named("eigenvectorsQR")   = vectorsQR,
                Rcpp::Named("eigenvalues") = values,
                Rcpp::Named("x") = t) ;
              break;}
    
              {case 6:
                arma::mat M6(nn,nn);
#pragma omp parallel for
                for (int i=0;i<nn; i++){
#pragma omp parallel for
                  for (int j=0;j<nn;j++){
                    if (j==0){ M6(i,j)=const0;}
                    else if (j==1){M6(i,j)= const1 + const2*t(i);}
                    else if (j==2){M6(i,j)= const3 + const4*t(i) + const5* pow(t(i),2);}
                    else if (j==3){M6(i,j)= const6 + const7*t(i) + const8* pow(t(i),2) 
                      + const9* pow(t(i),3);}
                    else if (j==4){M6(i,j)= const10 + const11*t(i) + const12* pow(t(i),2) 
                      + const13* pow(t(i),3) + const14* pow(t(i),4);}
                    else if (j==5){M6(i,j)= const15 + const16*t(i) + const17* pow(t(i),2) 
                      + const18* pow(t(i),3) + const19* pow(t(i),4) + const20*pow(t(i),5);}
                    else {dcomp ans =(3 + 2*sqrt(3))*(pow(comp1,1 + k(j))/pow(E,(-3.5 + k(j))*
                    Pi*(1 - t(i))) + pow(E,-((-3.5 + k(j))*Pi*t(i)))) +((1 + sqrt(3))/2 - 
                    comp5*(5.0 + 3*sqrt(3)))*(pow(comp1,1 + k(j))/pow(E,pow(comp1,1.0/6.0)*
                    (-3.5 + k(j))*Pi*(1 - t(i))) +pow(E,-(pow(comp1,1.0/6.0)*(-3.5 + k(j))*
                    Pi*t(i)))) + ((-3 - sqrt(3))/2 - comp5*(1 + sqrt(3)))*(pow(comp1,1 + k(j))/
                      pow(E,pow(comp1,1.0/3.0)*(-3.5 + k(j))*Pi*(1 - t(i))) + pow(E,-(pow(comp1,1.0/3.0)*
                        (-3.5 + k(j))*Pi*t(i)))) +((-3 - sqrt(3))/2 + comp5*(1 + sqrt(3)))* 
                        (pow(comp1,1 + k(j))/pow(E,(0.5 - comp5*sqrt(3))*(-3.5 + k(j))*Pi*(1 - t(i)))+ 
                        pow(E,-((0.5 - comp5*sqrt(3))*(-3.5 + k(j))*Pi*t(i)))) + ((1 + sqrt(3))/2 + 
                        comp5*(5.0 + 3*sqrt(3)))*(pow(comp1,1 + k(j))/pow(E,(comp2 + sqrt(3)/2)*
                        (-3.5 + k(j))*Pi*(1 - t(i))) +pow(E,-((comp2 + sqrt(3)/2)*(-3.5 + k(j))*Pi*t(i)))) 
                      - sqrt(2)*cos(Pi/4 + (-3.5 + k(j))*Pi*t(i));
                      M6(i,j) =ans.real();}
                    M6(i,j) =M6(i,j)/sqrt(nn);
                  }
                }
                arma:: mat vectors=M6; arma::mat vectorsQR=Turn(M6); arma::vec values=eigenvalues(nn,6);
                out=Rcpp::List::create(
                  Rcpp::Named("eigenvectors") = vectors,
                  Rcpp::Named("eigenvectorsQR")   = vectorsQR,
                  Rcpp::Named("eigenvalues") = values,
                  Rcpp::Named("x") = t) ;
                break;}
    
                {case 999:
                  arma::mat M1(nn,nn), M2(nn,nn), M3(nn,nn), M4(nn,nn), M5(nn,nn), M6(nn,nn);

#pragma omp parallel for
                  for (int i=0;i<nn; i++){
#pragma omp parallel for
                    for (int j=0;j<nn;j++){
                      if (j==0){ M1(i,j)=const0;}
                      else { M1(i,j) = sqrt(2)*cos((k(j)-1)*Pi*t(i));}
                      M1(i,j) =M1(i,j)/sqrt(nn);
                    }
                  }
#pragma omp parallel for
                  for (int i=0;i<nn; i++){
#pragma omp parallel for
                    for (int j=0;j<nn;j++){
                      if (j==0){ M2(i,j)=const0;}
                      else if (j==1){M2(i,j)= const1 + const2*t(i);}
                      else {dcomp ans = pow(comp1,1+k(j))/pow(E,(-1.5+k(j))*Pi*(1-t(i))) +  
                        pow(E,-((-1.5+k(j))*Pi*t(i))) + sqrt(2)*cos(Pi/4 +(-1.5 + k(j))*Pi*t(i));
                        M2(i,j) =ans.real();}
                      M2(i,j) =M2(i,j)/sqrt(nn);
                    }
                  }
#pragma omp parallel for
                  for (int i=0;i<nn; i++){
#pragma omp parallel for
                    for (int j=0;j<nn;j++){
                      if (j==0){ M3(i,j)=1;}
                      else if (j==1){M3(i,j)= const1 + const2*t(i);}
                      else if (j==2){M3(i,j)= const3 + const4*t(i) + const5* pow(t(i),2);}
                      else {dcomp ans = (sqrt(1.5)-comp0/sqrt(2))*((pow(comp1,1+k(j))/
                                         pow(E,(pow(comp1,1.0/6))*(-2.0+k(j))*Pi*(1-t(i))))+pow(E,-(pow(comp1,1.0/6))*
                                           ((-2.0+k(j))*Pi*t(i))))+(sqrt(1.5)+comp0/sqrt(2))*(pow(comp1,1+k(j))/
                                             pow(E,(comp2 +sqrt(3)/2)*(-2.0+k(j))*Pi*(1-t(i)))+pow(E,-((comp2 +sqrt(3)/2)*
                                               (-2.0+k(j))*Pi*t(i))))-sqrt(2)*sin((-2.0 + k(j))*Pi*t(i));
                        M3(i,j) =ans.real();}
                      M3(i,j) =M3(i,j)/sqrt(nn);
                    }
                  }
#pragma omp parallel for
                  for (int i=0;i<nn; i++){
#pragma omp parallel for
                    for (int j=0;j<nn;j++){
                      if (j==0){ M4(i,j)=const0;}
                      else if (j==1){M4(i,j)= const1 + const2*t(i);}
                      else if (j==2){M4(i,j)= const3 + const4*t(i) + const5* pow(t(i),2);}
                      else if (j==3){M4(i,j)= const6 + const7*t(i) + const8* pow(t(i),2) + 
                        const9* pow(t(i),3);}
                      else {dcomp ans =(1 + sqrt(2))*(pow(comp1,1+k(j))/pow(E,(-2.5+k(j))* 
                                        Pi*(1-t(i)))+ pow(E,-((-2.5+k(j))* Pi*t(i)))) +(1/sqrt(2) - 
                                        comp0*(1 + 1/sqrt(2)))*(pow(comp1,1+k(j))/pow(E,pow(comp1,0.25)*(-2.5+k(j))*
                                        Pi*(1-t(i))) +pow(E,-(pow(comp1,0.25)*(-2.5+k(j))*Pi*t(i))))+ (1/sqrt(2)+ 
                                        comp0*(1 + 1/sqrt(2)))*(pow(comp1,1+k(j))/pow(E,(comp3 *(-2.5 +k(j))*Pi*
                                        (1-t(i)))/sqrt(2)) +pow(E,(comp4 *(-2.5 +k(j))*Pi*t(i))/sqrt(2))) - 
                                        sqrt(2)*cos(Pi/4 - (-2.5 + k(j))*Pi*t(i));
                        M4(i,j) = ans.real();}
                      M4(i,j) = M4(i,j)/sqrt(nn);
                    }
                  }
                  
#pragma omp parallel for
                  for (int i=0;i<nn; i++){
#pragma omp parallel for
                    for (int j=0;j<nn;j++){
                      if (j==0){ M5(i,j)=const0;}
                      else if (j==1){M5(i,j)= const1 + const2*t(i);}
                      else if (j==2){M5(i,j)= const3 + const4*t(i) + const5* pow(t(i),2);}
                      else if (j==3){M5(i,j)= const6 + const7*t(i) + const8* pow(t(i),2) + 
                        const9* pow(t(i),3);}
                      else if (j==4){M5(i,j)= const10 + const11*t(i) + const12* pow(t(i),2) + 
                        const13* pow(t(i),3) + const14* pow(t(i),4);}
                      else {dcomp ans =(sqrt(2)*(1 + sqrt(5)/2)- (comp5*(sqrt(10 - 2*sqrt(5)) + 
                                        sqrt(2*(5 + sqrt(5)))))/sqrt(2))*(pow(comp1,1 + k(j))/pow(E,pow(comp1,0.1)
                                                                                                    *(dcomp(-3 + k(j)))*Pi*(1 - t(i))) + pow(E,-(pow(comp1,0.1)*(dcomp(-3 + k(j)))
                                                                                                                                                   *Pi*t(i)))) + (-(1/sqrt(2)) - (comp5*(sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + 
                                                                                                                                                   sqrt(5)))))/sqrt(2))*(pow(comp1,1 + k(j))/pow(E,pow(comp1,0.3)*(dcomp (-3 + k(j)))
                                                                                                                                                                                                   *Pi*(1 - t(i))) + pow(E,-(pow(comp1,0.3)*(dcomp (-3 + k(j)))*Pi*t(i))))+ (sqrt(2)*
                                                                                                                                                                                                     (1 + sqrt(5)/2) + (comp5*(sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + sqrt(5)))))/
                                                                                                                                                                                                       sqrt(2)) *(pow(comp1,1 + k(j))/pow(E,(sqrt(0.625 + sqrt(5)/8) - comp6*
                                                                                                                                                                                                         (-1 + sqrt(5)))*(dcomp (-3 + k(j)))*Pi*(1 - t(i))) +pow(E,-((sqrt(0.625 + sqrt(5)/8.) 
                                                                                                                                                                                                                                                                        - comp6*(-1 + sqrt(5)))*(dcomp (-3 + k(j)))*Pi*t(i)))) +(-(1/sqrt(2)) + 
                                                                                                                                                                                                                                                                          (comp5*(sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + sqrt(5)))))/sqrt(2)) *(pow(comp1,1 + k(j))
                                                                                                                                                                                                                                                                                                                                              /pow(E,(dcomp(sqrt(0.625 - sqrt(5)/8)) - comp6*(1 + sqrt(5)))*(dcomp (-3 + k(j)))
                                                                                                                                                                                                                                                                                                                                                     *Pi*(1 - t(i))) +pow(E,-((sqrt(0.625 - sqrt(5)/8) - comp6*(1 + sqrt(5)))
                                                                                                                                                                                                                                                                                                                                                     *(dcomp (-3 + k(j)))*Pi*t(i))))- sqrt(2)*cos((-3 + k(j))*Pi*t(i));
                                                                                                                                                                                                                                                                                                                                                     M5(i,j) =ans.real();}
                      M5(i,j) =M5(i,j)/sqrt(nn);
                    }
                  }
                  
#pragma omp parallel for
                  for (int i=0;i<nn; i++){
#pragma omp parallel for
                    for (int j=0;j<nn;j++){
                      if (j==0){ M6(i,j)=const0;}
                      else if (j==1){M6(i,j)= const1 + const2*t(i);}
                      else if (j==2){M6(i,j)= const3 + const4*t(i) + const5* pow(t(i),2);}
                      else if (j==3){M6(i,j)= const6 + const7*t(i) + const8* pow(t(i),2) 
                        + const9* pow(t(i),3);}
                      else if (j==4){M6(i,j)= const10 + const11*t(i) + const12* pow(t(i),2) 
                        + const13* pow(t(i),3) + const14* pow(t(i),4);}
                      else if (j==5){M6(i,j)= const15 + const16*t(i) + const17* pow(t(i),2) 
                        + const18* pow(t(i),3) + const19* pow(t(i),4) + const20*pow(t(i),5);}
                      else {dcomp ans =(3 + 2*sqrt(3))*(pow(comp1,1 + k(j))/pow(E,(-3.5 + k(j))*
                      Pi*(1 - t(i))) + pow(E,-((-3.5 + k(j))*Pi*t(i)))) +((1 + sqrt(3))/2 - 
                      comp5*(5.0 + 3*sqrt(3)))*(pow(comp1,1 + k(j))/pow(E,pow(comp1,1.0/6.0)*
                      (-3.5 + k(j))*Pi*(1 - t(i))) +pow(E,-(pow(comp1,1.0/6.0)*(-3.5 + k(j))*
                      Pi*t(i)))) + ((-3 - sqrt(3))/2 - comp5*(1 + sqrt(3)))*(pow(comp1,1 + k(j))/
                        pow(E,pow(comp1,1.0/3.0)*(-3.5 + k(j))*Pi*(1 - t(i))) + 
                          pow(E,-(pow(comp1,1.0/3.0)*(-3.5 + k(j))*Pi*t(i)))) +((-3 - sqrt(3))/2 + 
                          comp5*(1 + sqrt(3)))* (pow(comp1,1 + k(j))/pow(E,(0.5 - comp5*sqrt(3))*
                          (-3.5 + k(j))*Pi*(1 - t(i)))+ pow(E,-((0.5 - comp5*sqrt(3))*(-3.5 + k(j))*Pi*t(i)))) 
                        + ((1 + sqrt(3))/2 + comp5*(5.0 + 3*sqrt(3)))*(pow(comp1,1 + k(j))/
                          pow(E,(comp2 + sqrt(3)/2)*(-3.5 + k(j))*Pi*(1 - t(i))) +pow(E,-((comp2 + sqrt(3)/2)*
                          (-3.5 + k(j))*Pi*t(i)))) - sqrt(2)*cos(Pi/4 + (-3.5 + k(j))*Pi*t(i));
                        M6(i,j) =ans.real();}
                      M6(i,j) =M6(i,j)/sqrt(nn);
                    }
                  }

                  arma::mat vectors1=M1; arma::mat vectorsQR1=Turn(M1); arma::vec values1=eigenvalues(nn,1);
                  arma::mat vectors2=M2; arma::mat vectorsQR2=Turn(M2); arma::vec values2=eigenvalues(nn,2);
                  arma::mat vectors3=M3; arma::mat vectorsQR3=Turn(M3); arma::vec values3=eigenvalues(nn,3);
                  arma::mat vectors4=M4; arma::mat vectorsQR4=Turn(M4); arma::vec values4=eigenvalues(nn,4);
                  arma::mat vectors5=M5; arma::mat vectorsQR5=Turn(M5); arma::vec values5=eigenvalues(nn,5);
                  arma::mat vectors6=M6; arma::mat vectorsQR6=Turn(M6); arma::vec values6=eigenvalues(nn,6);
                  
                  out=Rcpp::List::create(
                    Rcpp::List::create(
                      Rcpp::Named("eigenvectors") = vectors1,
                      Rcpp::Named("eigenvectorsQR")   = vectorsQR1,
                      Rcpp::Named("eigenvalues") = values1,
                      Rcpp::Named("x") = t),
                      Rcpp::List::create(
                        Rcpp::Named("eigenvectors") = vectors2,
                        Rcpp::Named("eigenvectorsQR")   = vectorsQR2,
                        Rcpp::Named("eigenvalues") = values2,
                        Rcpp::Named("x") = t),
                        Rcpp::List::create(
                          Rcpp::Named("eigenvectors") = vectors3,
                          Rcpp::Named("eigenvectorsQR")   = vectorsQR3,
                          Rcpp::Named("eigenvalues") = values3,
                          Rcpp::Named("x") = t),
                          Rcpp::List::create(
                            Rcpp::Named("eigenvectors") = vectors4,
                            Rcpp::Named("eigenvectorsQR")   = vectorsQR4,
                            Rcpp::Named("eigenvalues") = values4,
                            Rcpp::Named("x") = t),
                            Rcpp::List::create(
                              Rcpp::Named("eigenvectors") = vectors5,
                              Rcpp::Named("eigenvectorsQR")   = vectorsQR5,
                              Rcpp::Named("eigenvalues") = values5,
                              Rcpp::Named("x") = t),
                              Rcpp::List::create(
                                Rcpp::Named("eigenvectors") = vectors6,
                                Rcpp::Named("eigenvectorsQR")   = vectorsQR6,
                                Rcpp::Named("eigenvalues") = values6,
                                Rcpp::Named("x") = t));
                  break;
                }
  }
  return out;
  
}

