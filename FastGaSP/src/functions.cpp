
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*- 
 
// we only include RcppEigen.h which pulls Rcpp.h in for us 
#include <RcppEigen.h> 
#include <Rcpp.h> 
#include <cmath> 
#include "ctools.h"
// [[Rcpp::depends(RcppEigen)]] 

using namespace Rcpp;
using namespace std;
using namespace Eigen; 
 
 

//July 16, 2018
////Construct_W0_matern_5_2 
// [[Rcpp::export]] 
MatrixXd Construct_W0_matern_5_2(const double sigma2, const double lambda){ 
  //int num_dim=sigma2.size(); 
  
  Eigen::MatrixXd W0= Eigen::MatrixXd::Zero(3,3); 
  //Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  W0(0,0)=sigma2; 
  W0(0,2)=W0(2,0)=-sigma2*pow(lambda,2.0)/3.0; 
  W0(1,1)=sigma2*pow(lambda,2.0)/3.0; 
  W0(2,2)=sigma2*pow(lambda,4.0); 

  return W0; 
} 

//July 16, 2018
////Construct_W0_exp 
// [[Rcpp::export]] 
MatrixXd Construct_W0_exp(const double sigma2, const double lambda){ 
  //int num_dim=sigma2.size(); 
  
  Eigen::MatrixXd W0= Eigen::MatrixXd::Zero(1,1); 

  W0(0,0)=sigma2; 
  
  return W0; 
} 


 
 ////Construct_G_matern_5_2 
 // [[Rcpp::export]] 
 List Construct_G_matern_5_2(Eigen::VectorXd delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
   int num_obs=delta_x.size()+1; 
   //int num_dim=lambda.size();  
   List GG(num_obs);  
   GG[0]=Eigen::MatrixXd::Zero(3,3); 
   
   Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
   
   // num_dim list, each is 3(num_obs)\times 3 list 
  // for(int i_GG=0;i_GG<num_dim;i_GG++){ 
  //   Eigen::MatrixXd d= Eigen::MatrixXd::Zero(num_obs,9);  //the first row has all zeros  
     for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
       int j_GG_1=j_GG+1;    
       d(0,0)=pow(lambda,2.0)*pow(delta_x[j_GG],2.0)+2*lambda*delta_x[j_GG]+2; 
       d(1,0)=-pow(lambda,3.0)*pow(delta_x[j_GG],2.0); 
       d(2,0)=pow(lambda,4.0)*pow(delta_x[j_GG],2.0)-2*pow(lambda,3.0)*delta_x[j_GG]; 
       d(0,1)=2*(lambda*pow(delta_x[j_GG],2.0)+delta_x[j_GG]); 
       d(1,1)=-2*(pow(lambda,2.0)*pow(delta_x[j_GG],2.0)-lambda*delta_x[j_GG]-1); 
       d(2,1)=2*(pow(lambda,3.0)*pow(delta_x[j_GG],2.0)-3*pow(lambda,2.0)*delta_x[j_GG]); 
       d(0,2)=pow(delta_x[j_GG],2); 
       d(1,2)=2*delta_x[j_GG]-lambda*pow(delta_x[j_GG],2.0); 
       d(2,2)=pow(lambda,2.0)*pow(delta_x[j_GG],2.0)-4*lambda*delta_x[j_GG]+2;     
       d=exp(-lambda*delta_x[j_GG])/2.0*d;
       GG[j_GG_1]=d; 
     } 
    // GG[i_GG]=d; 
   //} 
   return GG; 
} 

////Construct_G_exp
// [[Rcpp::export]] 
List Construct_G_exp(Eigen::VectorXd delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=lambda.size();  
  List GG(num_obs);  
  GG[0]=Eigen::MatrixXd::Zero(1,1); 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(1,1); 
  
  for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
    d(0,0)=exp(-delta_x[j_GG]*lambda);
    GG[j_GG+1]=d; 
  }

  return GG;
}

////Construct_W_matern_5_2  
// [[Rcpp::export]] 
List Construct_W_matern_5_2(double sigma2,Eigen::VectorXd delta_x, double lambda, MatrixXd W0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=sigma2.size();  
  List Wi(num_obs);  
  Wi[0]=W0; 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  
 // List Wi(num_obs); 
  //for(int i_Wi=0;i_Wi<num_dim;i_Wi++){ 
    //Eigen::MatrixXd d= Eigen::MatrixXd::Zero(num_obs,9);   
    double  lambda_delta_x;
    double exp_neg_2_lambda_delta_x;
    int  j_Wi_1;
    for(int j_Wi=0;j_Wi<(num_obs-1);j_Wi++){ 
      j_Wi_1= j_Wi+1; 
      lambda_delta_x=lambda*delta_x[j_Wi];  //close and jump then it is... 
      exp_neg_2_lambda_delta_x=exp(-2*lambda_delta_x); 
      
      d(0,0)=(exp_neg_2_lambda_delta_x*(3+6*lambda_delta_x+6*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-3 )/(-4*pow(lambda,5.0)); 
      d(1, 0)=  d(0, 1)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Wi],4.0)/2.0; 
      d(2, 0)=  d(0, 2)=(exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))-1 )/(4*pow(lambda,3.0)); 
      d(1, 1)= (exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)-4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-1 )/(-4*pow(lambda,3.0)); 
      d(2, 1)=  d(1, 2)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Wi],2.0)*(4-4*lambda_delta_x+pow(lambda_delta_x,2.0) )/2.0; 
      d(2, 2)=(exp_neg_2_lambda_delta_x*(-3+10*lambda_delta_x-22*pow(lambda_delta_x,2.0)+12*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))+3 )/(4*lambda)     ;  
      d=d*(4*sigma2*pow(lambda,5.0)/3.0); 
      Wi[j_Wi_1]=d; 
      
    //} 
  } 
  return Wi; 
}


////Construct_W_matern_5_2  
// [[Rcpp::export]] 
List Construct_W_exp(double sigma2, Eigen::VectorXd delta_x, double lambda, MatrixXd W0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=sigma2.size();  
  List Wi(num_obs);  
  Wi[0]=W0; 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(1,1); 
  
  for(int j_Wi=0;j_Wi<(num_obs-1);j_Wi++){ 
    d(0,0)=1-exp(-2*delta_x[j_Wi]*lambda);
    Wi[j_Wi+1]=d;
  }
  
  return Wi;
}

////Get_Q_K  
// [[Rcpp::export]] 
List Get_Q_K(const List GG,const List  W,const Eigen::MatrixXd C0,const double VV){ 

   int n=GG.size();
   int k=C0.rows();
   
   Eigen::VectorXd Q=Eigen::VectorXd::Zero(n);
   Eigen::MatrixXd K=Eigen::MatrixXd::Zero(n,k);
   Eigen::MatrixXd C=C0;
    
   Eigen::MatrixXd GG_matrix;
   Eigen::MatrixXd W_matrix;
   
   Eigen::MatrixXd RR;
   
      
   // num_dim list, each is 3(num_obs)\times 3 list 
   for(int t=0;t<n;t++){ 
     GG_matrix=GG[t];
     W_matrix=W[t];
     RR=GG_matrix*C*GG_matrix.transpose()+W_matrix;
     //Q[t]=RR(0,0);
     Q[t]=RR(0,0)+VV;
     K.row(t)=RR.col(0).transpose()/Q[t];
     C=RR-RR.col(0)*RR.row(0)/Q[t];
   }

   List return_list;
   return_list.push_back(Q);
   return_list.push_back(K);
   
   return return_list;
}



////Get_log_det_S2
// [[Rcpp::export]] 
List Get_log_det_S2(const Eigen::VectorXd param,const bool have_noise,const Eigen::VectorXd delta_x,const Eigen::VectorXd output,
                    const  String kernel_type){
  
  //int n1=output_KF.rows();
  //int n2=output_KF.cols();
  
  int n=output.rows();
  
  double gamma=1.0/exp(param[0]);
  
  double VV=0;
  if(have_noise){
    VV=exp(param[1]);
  }
  
  Eigen::MatrixXd    W0;
  List    GG;
  
  List    W;
  List    Q_K;
  
  
  double lambda=0;
  if(kernel_type=="matern_5_2"){
    lambda=sqrt(5.0)/gamma;
  
  //param.tail(k).array().exp().matrix();

  W0=Construct_W0_matern_5_2(1.0,lambda);   
  GG=Construct_G_matern_5_2(delta_x,lambda);  
  W=Construct_W_matern_5_2(1.0,delta_x,lambda,W0);
  
  }else if(kernel_type=="exp"){
    
    lambda=1.0/gamma;
    W0=Construct_W0_exp(1.0,lambda);  
    GG=Construct_G_exp(delta_x,lambda);  
    W=Construct_W_exp(1.0,delta_x,lambda,W0);
    
  }

  Q_K=Get_Q_K(GG,W,W0,VV);
  
  
  Eigen::VectorXd Q=Q_K[0];
  Eigen::MatrixXd K=Q_K[1];
  
  double log_det_R=Q.array().log().sum();
  List return_vec;
  return_vec.push_back(log_det_R);
    
  //return_list.push_back(log_det_R);
  
  Eigen::MatrixXd GG_matrix;
  
  Eigen::MatrixXd m=Eigen::VectorXd::Zero(n);
  
  Eigen::VectorXd a;
  
  //Eigen::VectorXd Y_minus_a_1_scaled_vec=Eigen::VectorXd::Zero(n); 
  Eigen::VectorXd Y_minus_a_1=Eigen::VectorXd::Zero(n); 
  
  Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  for(int t=0;t<n;t++){
    GG_matrix=GG[t];
    a=GG_matrix*m;
    Y_minus_a_1[t]=(output[t]-a[0]);
    
   // Y_minus_a_1_scaled_vec[t]=(output[t]-a[0])/sqrt_Q[t];
    m=a+K.row(t).transpose()*(output[t]-a[0]);
  }
  

  ///sigma2_S2=(FRFt+sigma_2*eta)  ####using this will lead to |R|
    
    
  double S2=(Y_minus_a_1.array()*Y_minus_a_1.array()/Q.array()).sum(); 
    
  return_vec.push_back(S2);
  
  return return_vec;
}


//double log_lik_GP_fast(const vectorXd param, const VectorXd delta_x, const VectorXd output){
//  
//}


/////////the following code is for prediction 

//now consider about the prediction
////Get_C_R_K_pred, C and R is for smoothing, K is for filtering the mean  
// [[Rcpp::export]] 
List Get_C_R_K_pred(const VectorXi index, const List GG,const List  W,const Eigen::MatrixXd C0, double VV){ 
  
  
  //index is a sequence where 0 means NA and 1 means observations
  int n=GG.size();//total number
  int k=C0.rows();
  
  //Eigen::VectorXd Q=Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd K=Eigen::MatrixXd::Zero(n,k);
  List C(n+1);
  C[0]=C0;
  
  Eigen::MatrixXd GG_matrix;
  Eigen::MatrixXd W_matrix;
  Eigen::VectorXd Q=Eigen::VectorXd::Zero(n);
  
  Eigen::MatrixXd RR;
  
  List R(n);
  
  Eigen::MatrixXd C_cur=C[0];
  //int index_num=0;

  for(int t=0;t<n;t++){
    
    if(index[t]==1){
      //index_num=index_num+1;
      GG_matrix=GG[t];
      W_matrix=W[t];
      RR=GG_matrix*C_cur*GG_matrix.transpose()+W_matrix;
      R[t]=RR;
      
      //Q[t]=RR[t](0,0);
      Q[t]=RR(0,0)+VV;
      
      
      K.row(t)=RR.col(0).transpose()/Q[t];
      C[t+1]=RR-RR.col(0)*RR.row(0)/Q[t];
      C_cur=C[t+1];
      
    }else{
      GG_matrix=GG[t];
      W_matrix=W[t];
      
      R[t]=GG_matrix*C_cur*GG_matrix.transpose()+W_matrix;
      C[t+1]=C_cur=R[t];   
      
    }
  }
  
  List return_list;
  return_list.push_back(C);
  return_list.push_back(R);
  return_list.push_back(K);
  
  return return_list;
}




////Get_Y_minus_a_1_scaled_matrix_2d  
// [[Rcpp::export]] 
List Get_m_a_pred(const VectorXi index, const Eigen::VectorXd output_vec,const List GG,const Eigen::MatrixXd K){
  
  //output_KF
  int n=GG.size();//total number
  
  //int n1=output_KF.rows();
  //int n2=output_KF.cols();
  int k=K.cols();
  
  List m(n);  //m should have n+1 item but the first one is (0,0,0) so it's omitted.
  List a(n);
  Eigen::VectorXd a_vec;
  Eigen::VectorXd m_cur=Eigen::VectorXd::Zero(k); //3\times 1
  
  // =Eigen::MatrixXd::Zero(k,n2);
  
  //Eigen::MatrixXd Y_minus_a_1_scaled_matrix=Eigen::MatrixXd::Zero(n1,n2); 
  //Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  Eigen::MatrixXd GG_matrix;
  
  //int index_num=0;
  int index_count=0;
  for(int t=0;t<n;t++){
    if(index[t]==1){
      //index_num=index_num+1;
      GG_matrix=GG[t];
      a_vec=GG_matrix*m_cur;
      a[t]=a_vec;
      m[t]=a_vec+K.row(t).transpose()*(output_vec[index_count]-a_vec[0]);
      m_cur=m[t];
      index_count=index_count+1;
    }else{
      GG_matrix=GG[t];
      m[t]=a[t]=GG_matrix*m_cur;
      m_cur=m[t];
    }
  }
  
  List return_list;
  
  return_list.push_back(m);
  return_list.push_back(a);
  
  return return_list;
}


////Kalman_smoother, the index is 0 and 1 where 0 means missing value
// [[Rcpp::export]] 
List Get_S_KK(const VectorXi index,const List GG,const List C,const List R){
  
  int n=GG.size();
  //int k=3; //matern_5_2
  
  List S(n);
  List KK(n-1);
  
  S[n-1]= C[n];
  
  Eigen::MatrixXd GG_matrix;
  Eigen::MatrixXd C_matrix;
  Eigen::MatrixXd R_matrix;
  Eigen::MatrixXd KK_matrix;
  
  Eigen::MatrixXd S_cur=C[n];
  for(int t=n-2;t>=0;t--){
    GG_matrix=GG[t+1];
    R_matrix=R[t+1];
    C_matrix=C[t+1];  // C should have n+1 items
    
    LLT<MatrixXd> lltOfR(R_matrix);    // compute the cholesky decomposition of R called lltofR
    MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition
    
    KK_matrix= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(GG_matrix*C_matrix))); 
    KK[t]=KK_matrix;
    
    S[t]=C_matrix-KK_matrix.transpose()*(R_matrix-S_cur)*KK_matrix;
    S_cur=S[t];
  }
  
  List return_list;
  
  return_list.push_back(S);
  return_list.push_back(KK);
  
  return return_list;
  
}

// [[Rcpp::export]] 
MatrixXd Get_s_1st(const List m,const List a,const List C,const List KK){
  
  int n=C.size()-1;
  //int k=3; //matern_5_2
  
  Eigen::VectorXd s=m[n-1];
  //int n_2=s.cols();
  VectorXd s_1st=Eigen::VectorXd::Zero(n);
  MatrixXd KK_matrix;
  s_1st[n-1]=s[0];
  VectorXd a_vec;
  VectorXd m_vec;
  for(int t=n-2;t>=0;t--){
    KK_matrix=KK[t];
    a_vec=a[t+1];
    m_vec=m[t];
    s=m_vec+KK_matrix.transpose()*(s-a_vec);
    s_1st[t]=s[0];
  }
  
  return s_1st; 
}




// [[Rcpp::export]] 
List Kalman_smoother(const VectorXd param,const bool have_noise,const VectorXi index_obs, 
                     const VectorXd delta_x_all, const VectorXd output, const double sigma_2_hat,
                     const String kernel_type){
  //int n1=output_KF.rows();
  //int n2=output_KF.cols();
  
  double gamma=1.0/exp(param[0]);
  //double lambda=sqrt(5.0)/gamma;
  
  
  double VV=0;
  if(have_noise){
    VV=exp(param[1]); 
  }
  int n=delta_x_all.size()+1;
  //param.tail(k).array().exp().matrix();
  
  
  Eigen::MatrixXd W0;
  List GG;
  List W;
  
  double lambda=0;
  if(kernel_type=="matern_5_2"){
    lambda=sqrt(5.0)/gamma;
    
    //param.tail(k).array().exp().matrix();
    
    W0=Construct_W0_matern_5_2(1.0,lambda);   
    GG=Construct_G_matern_5_2(delta_x_all,lambda);  
    W=Construct_W_matern_5_2(1.0,delta_x_all,lambda,W0);
    
  }else if(kernel_type=="exp"){
    
    lambda=1.0/gamma;
    W0=Construct_W0_exp(1.0,lambda);  
    GG=Construct_G_exp(delta_x_all,lambda);  
    W=Construct_W_exp(1.0,delta_x_all,lambda,W0);
    
  }
  
  
  
  List C_R_K=Get_C_R_K_pred(index_obs, GG,W, W0, VV);
  
  
  List m_a=Get_m_a_pred(index_obs, output, GG,C_R_K[2]);
  
  List S_KK=Get_S_KK(index_obs, GG,C_R_K[0],C_R_K[1]);
  
  
  MatrixXd s_1st=Get_s_1st(m_a[0],m_a[1], C_R_K[0],S_KK[1]);
  
  List return_list;
  
  List S= S_KK[0];
  MatrixXd S_matrix;
  VectorXd Var=VectorXd::Zero(n);
  for(int t=0;t<n;t++){
    S_matrix=S[t];
    Var[t]=S_matrix(0,0);
  }
  Var=((Var.array()+VV)*sigma_2_hat).matrix();
  
  return_list.push_back(s_1st);
  return_list.push_back(Var);
  /*
  List return_list;
  
  return_list.push_back(W);
  */
  return return_list;
  
}




