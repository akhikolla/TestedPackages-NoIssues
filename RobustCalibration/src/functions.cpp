
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*- 
 
// we only include RcppEigen.h which pulls Rcpp.h in for us 
#include <RcppEigen.h> 
#include <Rcpp.h> 
#include <cmath>
#include "ctools.h"
//#include <chrono>
//#include <random>


// [[Rcpp::depends(RcppEigen)]] 

using namespace Rcpp;
using namespace std;
using namespace Eigen; 
 
 #define pi 3.141592653589793238462643383279502884197169399375105820974944592307816406L

//Cholesky decomposition in general
// [[Rcpp::export]]
MatrixXd Chol_Eigen(const Eigen::MatrixXd R){
     	    
      LLT<MatrixXd> lltOfR(R);            
      MatrixXd L =  lltOfR.matrixL();   //retrieve factor L  in the decomposition
      return L;
}

  
///////Kilauea computer model
// [[Rcpp::export]]
Eigen::VectorXd Mogihammer ( const MatrixXd obsCoords, const VectorXd m, int simul_type){
  VectorXd m_rev=m;
  Vector3d lookvec1(-0.616111273, -0.114189475, 0.779337965);
  Vector3d lookvec2(0.650447484, -0.119379453, 0.75011107);

  
  // m_rev(3)=m_rev(3)*17971200.0;
  //m_rev(3)=m_rev(3)*8600400.0;
  m_rev(3)=m_rev(3)*31557600.0;

  int nObs=obsCoords.rows(); 
  Eigen::VectorXd dx=obsCoords.col(0)-MatrixXd::Constant(nObs,1,m_rev(0));
  Eigen::VectorXd dy=obsCoords.col(1)-MatrixXd::Constant(nObs,1,m_rev(1));
  Eigen::VectorXd dd=MatrixXd::Constant(nObs,1,m_rev(2));
   
    Eigen::VectorXd R=(dx.array().pow(2)+dy.array().pow(2)+dd.array().pow(2)).matrix();

  Eigen::VectorXd K=MatrixXd::Constant(nObs,1,(1-m_rev(4))*m_rev(3)/pi);
  K=(K.array()/R.array().pow(1.5)).matrix();

  
  Eigen::MatrixXd A=MatrixXd::Zero(nObs,3);
  A.col(0)=K.cwiseProduct(dx);
  A.col(1)=K.cwiseProduct(dy);
  A.col(2)=K.cwiseProduct(dd);

    
  if(simul_type==2){
     return A*lookvec1;

  }else{
     return A*lookvec2;

  }

  // return A*lookvec;  
  //return K;

}

// [[Rcpp::export]]
bool Accept_proposal(double r){
  if(r>=1){
    return true;
  }else{

     //construct a trivial random generator engine from a time-based seed:
     //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
     //std::default_random_engine generator (seed);

     //std::uniform_real_distribution<> dis(0.0, 1.0);

     //double u=dis(generator);
     
     double u= R::runif(0.0,1.0);
     if(u<r){
       return true;
     }else{
       return false;
     }
  }

}

double Log_approx_ref_prior(const VectorXd param,double nugget, bool nugget_est, const Eigen::VectorXd CL,const double a,const double b ){

  Eigen::VectorXd beta;
  double nu=nugget;
  int param_size=param.size();
  if(!nugget_est){
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  double t=CL.cwiseProduct(beta).sum()+nu;
  double part_I=-b*t;
  double part_II= a*log(t);
  return part_I+part_II;
}


Eigen::MatrixXd Matern_5_2_funct (const MatrixXd d, double beta_i){
  //inline static Mat Matern_5_2_funct (const Eigen::Map<Eigen::MatrixXd> & d, double beta_i){
  const double cnst = sqrt(5.0);
  Eigen::MatrixXd matOnes = Eigen::MatrixXd::Ones(d.rows(),d.cols());
  Eigen::MatrixXd result = cnst*beta_i*d;
  return ((matOnes + result +
	   result.array().pow(2.0).matrix()/3.0).cwiseProduct((-result).array().exp().matrix()));
  
}

Eigen::MatrixXd Matern_3_2_funct (const MatrixXd d, double beta_i){
  const double cnst = sqrt(3.0);
  Eigen::MatrixXd matOnes = Eigen::MatrixXd::Ones(d.rows(),d.cols());
  Eigen::MatrixXd result = cnst*beta_i*d;
  return ((matOnes + result ).cwiseProduct((-result).array().exp().matrix()));
  
}


 Eigen::MatrixXd Pow_exp_funct (const MatrixXd d, double beta_i,double alpha_i){
  
  return (-(beta_i*d).array().pow(alpha_i)).exp().matrix();

}


Eigen::MatrixXd Separable_kernel (List R0, Eigen::VectorXd beta,String kernel_type, Eigen::VectorXd alpha ){
  Eigen::MatrixXd R0element = R0[0];
  int Rnrow = R0element.rows();
  int Rncol = R0element.cols();

  Eigen::MatrixXd R = R.Ones(Rnrow,Rncol);
  if(kernel_type=="matern_5_2"){
    for (int i_ker = 0; i_ker < beta.size(); i_ker++){
      R = (Matern_5_2_funct(R0[i_ker],beta[i_ker])).cwiseProduct(R);
    }
  }else if(kernel_type=="matern_3_2"){
    for (int i_ker = 0; i_ker < beta.size(); i_ker++){
      R = (Matern_3_2_funct(R0[i_ker],beta[i_ker])).cwiseProduct(R);
    }
  }
  else if(kernel_type=="pow_exp"){
    for (int i_ker = 0; i_ker < beta.size(); i_ker++){
      R = (Pow_exp_funct(R0[i_ker],beta[i_ker],alpha[i_ker])).cwiseProduct(R);
    }
  }
  return R;
}




////Get_R_z_new
// [[Rcpp::export]] 
MatrixXd Get_R_z_new(const Eigen::VectorXd beta_delta, const double  eta_delta, const double tilde_lambda, const  List R0, const  String kernel_type, const Eigen::VectorXd alpha,const  Eigen::VectorXd inv_output_weights){


    MatrixXd R_00= Separable_kernel(R0,beta_delta, kernel_type,alpha);
    
    int num_obs=R_00.cols();
    
    MatrixXd B=R_00+1.0/tilde_lambda*MatrixXd::Identity(num_obs,num_obs);
    LLT<MatrixXd> lltOfB(B);             // compute the cholesky decomposition of R called lltofB
  MatrixXd L_B = lltOfB.matrixL();   //retrieve factor L  in the decomposition      
  MatrixXd R_z= R_00- R_00*L_B.transpose().triangularView<Upper>().solve(L_B.triangularView<Lower>().solve(R_00)); //one forward and one backward to compute R.inv%*%X

  
  MatrixXd mat=  inv_output_weights.asDiagonal();

    
    // MatrixXd R_z_tilde=R_z+eta_delta*inv_output_weights_matrix;
  MatrixXd R_z_tilde=R_z+eta_delta*mat;

	    
  LLT<MatrixXd> lltOfR_z_tilde(R_z_tilde);            
  MatrixXd L_R_z =  lltOfR_z_tilde.matrixL();   //retrieve factor L  in the decomposition
  return L_R_z;

}

// Get_R_new
// [[Rcpp::export]] 
MatrixXd Get_R_new(const Eigen::VectorXd beta_delta, const double  eta_delta, const  List R0,  const  String kernel_type, const Eigen::VectorXd alpha,const  Eigen::VectorXd inv_output_weights){

      MatrixXd R= Separable_kernel(R0,beta_delta, kernel_type,alpha);
      MatrixXd mat=  inv_output_weights.asDiagonal();

      MatrixXd R_tilde=R+eta_delta*mat;

     	    
      LLT<MatrixXd> lltOfR_tilde(R_tilde);            
      MatrixXd L =  lltOfR_tilde.matrixL();   //retrieve factor L  in the decomposition
      return L;
}


//sample sigma_2 and theta_m
// [[Rcpp::export]]
Eigen::VectorXd Sample_sigma_2_theta_m(const Eigen::VectorXd  param, const Eigen::MatrixXd L_cur, const Eigen::VectorXd output,  const int p_theta, const int p_x, Eigen::MatrixXd X, bool have_mean,const VectorXd cm_obs){
  int num_obs=output.rows();
  //int p_x=input.cols();

  VectorXd theta=param.head(p_theta);
  VectorXd beta_delta=(param.segment(p_theta,p_x)).array().exp().matrix();
  // double sigma_2_delta =exp(param(p_theta+p_x));


  VectorXd theta_m=VectorXd::Zero(1);
  int p_theta_m=0;
  if(have_mean){
    p_theta_m=X.cols();
    theta_m=param.tail(p_theta_m);
  }


  
  // computer model outputs

  
  // VectorXd cm_obs=Mogihammer(theta,input,lookvec);


  VectorXd  output_tilde=output-cm_obs;
  VectorXd  output_tilde_normalized=output_tilde;

  if(have_mean){
       output_tilde_normalized=output_tilde-X*theta_m;
  }
  // return output_tilde_normalized;

  
    MatrixXd R_inv_y=L_cur.transpose().triangularView<Upper>().solve(L_cur.triangularView<Lower>().solve(output_tilde_normalized));

    VectorXd S_2_matrix=output_tilde_normalized.transpose()*R_inv_y;
    double S_2=S_2_matrix(0);
    //S_2_matrix(0)=S_2;
    //return S_2_matrix;
    
    //sample sigma_2
   //construct a trivial random generator engine from a time-based seed:
  //unsigned seed_1 = std::chrono::system_clock::now().time_since_epoch().count();
  //std::default_random_engine generator_1 (seed_1);

  //std::gamma_distribution<double> distribution_1 (num_obs/2.0,2.0/S_2); // this is k, theta so in the usual way we should have (alpha, 1/beta)
  //double phi_sample=distribution_1(generator_1);
     
  double phi_sample= R::rgamma(num_obs/2.0, 2.0/S_2);
     
  double  sigma_2_sample=1/phi_sample;

  VectorXd ans=VectorXd::Zero(1+p_theta_m);
  ans(0)=sigma_2_sample;

  //if needed, generate theta_m
  if(have_mean){
      MatrixXd R_inv_X=L_cur.transpose().triangularView<Upper>().solve(L_cur.triangularView<Lower>().solve(X));
      MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X;


      LLT<MatrixXd> lltOfXt_R_inv_X(Xt_R_inv_X);            
      MatrixXd LX =  lltOfXt_R_inv_X.matrixL();   

      MatrixXd Xt_R_inv_X_inv_Xt_R_inv=LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose()));

      VectorXd theta_m_hat=Xt_R_inv_X_inv_Xt_R_inv*output_tilde;

      //unsigned seed_2 = std::chrono::system_clock::now().time_since_epoch().count();
      //std::default_random_engine generator_2 (seed_2);

      //std::normal_distribution<double> distribution_2 (0.0,1.0);

      VectorXd random_norm=VectorXd::Zero(p_theta_m);
      for (int i=0; i<p_theta_m; ++i){
         random_norm(i)= R::rnorm(0,1.0);
           //distribution_2(generator_2); 
      }
      VectorXd theta_m_sample=theta_m_hat+((LX.inverse()*random_norm).array()*sqrt(sigma_2_sample)).matrix();
      
      ans.tail(p_theta_m)=theta_m_sample;
  }
   
  return ans;
  
 }

//posterior for gasp_z

// [[Rcpp::export]] 
double Log_marginal_post(const Eigen::VectorXd  param,Eigen::MatrixXd L_cur, const Eigen::VectorXd output,  const int p_theta,int p_x, Eigen::MatrixXd X, bool have_mean,const VectorXd CL,const double a,const double b, const VectorXd cm_obs){
  int num_obs=output.rows();
  //int p_x=input.cols();

  VectorXd theta=param.head(p_theta);
  // VectorXd beta_delta=(param.segment(p_theta,p_x)).array().exp().matrix();

  double sigma_2_delta =param(p_theta+p_x+1);


  
  VectorXd theta_m=VectorXd::Zero(1);
  int p_theta_m=0;
  if(have_mean){
    p_theta_m=X.cols();
    theta_m=param.tail(p_theta_m);
  }

    
  // computer model outputs
  //VectorXd theta_rev=theta;
  //theta_rev(3)=theta_rev(3)*17971200.0;

  //VectorXd cm_obs=Mogihammer(theta_rev,input,lookvec);


  VectorXd  output_tilde=output-cm_obs;
  VectorXd  output_tilde_normalized=output_tilde;

  if(have_mean){
       output_tilde_normalized=output_tilde-X*theta_m;
  }
  // return output_tilde_normalized;

  
    MatrixXd R_inv_y=L_cur.transpose().triangularView<Upper>().solve(L_cur.triangularView<Lower>().solve(output_tilde_normalized));

    VectorXd S_2_matrix=output_tilde_normalized.transpose()*R_inv_y;
    double S_2=S_2_matrix(0);

    double log_post=-num_obs/2.0*log(sigma_2_delta)-L_cur.diagonal().array().log().matrix().sum()-S_2/(2.0*sigma_2_delta)+log(sigma_2_delta)+Log_approx_ref_prior(param.segment(p_theta,p_x+1), 0,true,CL,a,b)+ (param.segment(p_theta,p_x+1)).sum();
    
    return log_post;

}





//sample sigma_2 and theta_m
// [[Rcpp::export]]
Eigen::VectorXd Sample_sigma_2_theta_m_no_discrepancy(const Eigen::VectorXd  param, const Eigen::VectorXd output,  const int p_theta, Eigen::MatrixXd X, bool have_mean, VectorXd inv_output_weights, const VectorXd cm_obs){
  int num_obs=output.rows();

  VectorXd theta=param.head(p_theta);
  // VectorXd beta_delta=(param.segment(p_theta,p_x)).array().exp().matrix();
  // double sigma_2_delta =exp(param(p_theta+p_x));


  VectorXd theta_m=VectorXd::Zero(1);
  int p_theta_m=0;
  if(have_mean){
    p_theta_m=X.cols();
    theta_m=param.tail(p_theta_m);
  }




  VectorXd  output_tilde=output-cm_obs;
  VectorXd  output_tilde_normalized=output_tilde;

  if(have_mean){
       output_tilde_normalized=output_tilde-X*theta_m;
  }
  // return output_tilde_normalized;

  double S_2=output_tilde_normalized.transpose()*((output_tilde_normalized.array()/inv_output_weights.array()).matrix());

    
    //sample sigma_2
   //construct a trivial random generator engine from a time-based seed:
  //unsigned seed_1 = std::chrono::system_clock::now().time_since_epoch().count();
  //std::default_random_engine generator_1 (seed_1);

  //std::gamma_distribution<double> distribution_1 (num_obs/2.0,2.0/S_2); // this is k, theta so in the usual way we should have (alpha, 1/beta)
  double phi_sample= R::rgamma(num_obs/2.0,2.0/S_2);
  //=distribution_1(generator_1);
     
  double  sigma_2_sample=1/phi_sample;

  VectorXd ans=VectorXd::Zero(1+p_theta_m);
  ans(0)=sigma_2_sample;

  //if needed, generate theta_m
  if(have_mean){
    //area size
    MatrixXd Ones=MatrixXd::Constant(1,num_obs,1.0);
    
    MatrixXd mat=  (Ones.array()/inv_output_weights.array()).matrix().asDiagonal();
    
    //MatrixXd mat=  inv_output_weights.asDiagonal();

    MatrixXd Xt_X=X.transpose()*mat*X;  
      

    LLT<MatrixXd> lltOfXt_X(Xt_X);            
    MatrixXd LX =  lltOfXt_X.matrixL();   
    
      MatrixXd Xt_X_inv_Xt=LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(X.transpose()*mat));

      VectorXd theta_m_hat=Xt_X_inv_Xt*output_tilde;

      //unsigned seed_2 = std::chrono::system_clock::now().time_since_epoch().count();
      //std::default_random_engine generator_2 (seed_2);

      //std::normal_distribution<double> distribution_2 (0.0,1.0);

      VectorXd random_norm=VectorXd::Zero(p_theta_m);
      for (int i=0; i<p_theta_m; ++i){
         random_norm(i)= R::rnorm(0,1.0);
           //distribution_2(generator_2); 
      }
      VectorXd theta_m_sample=theta_m_hat+((LX.inverse()*random_norm).array()*sqrt(sigma_2_sample)).matrix();
      
      ans.tail(p_theta_m)=theta_m_sample;
  }
   
  return ans;
  
 }

//posterior for gasp_z
// [[Rcpp::export]]
double Log_marginal_post_no_discrepancy(const Eigen::VectorXd  param, const Eigen::VectorXd output,  const int p_theta, Eigen::MatrixXd X, bool have_mean, VectorXd inv_output_weights,const VectorXd cm_obs){
  int num_obs=output.rows();
  //int p_x=input.cols();

  VectorXd theta=param.head(p_theta);
  // VectorXd beta_delta=(param.segment(p_theta,p_x)).array().exp().matrix();

  double sigma_2_delta =param(p_theta);


  
  VectorXd theta_m=VectorXd::Zero(1);
  int p_theta_m=0;
  if(have_mean){
    p_theta_m=X.cols();
    theta_m=param.tail(p_theta_m);
  }

    


  VectorXd  output_tilde=output-cm_obs;
  VectorXd  output_tilde_normalized=output_tilde;

  if(have_mean){
       output_tilde_normalized=output_tilde-X*theta_m;
  }

  double S_2=output_tilde_normalized.transpose()*((output_tilde_normalized.array()/inv_output_weights.array()).matrix());

  double log_post=-num_obs/2.0*log(sigma_2_delta)-S_2/(2.0*sigma_2_delta)+log(sigma_2_delta);
    
    return log_post;

}







//MCMC
// [[Rcpp::export]]
MatrixXd Update_R_inv_y(VectorXd R_inv_y, List R0, VectorXd beta_delta, String kernel_type, VectorXd alpha, double tilde_lambda,int num_obs){
    
     MatrixXd  R= Separable_kernel(R0,beta_delta, kernel_type,alpha);

     MatrixXd R_middle=MatrixXd::Identity(num_obs,num_obs)+ tilde_lambda*R;
              	    
      LLT<MatrixXd> lltOfR_middle(R_middle);            
      MatrixXd  L_middle =  lltOfR_middle.matrixL();

     VectorXd  R_inv_y_updated= L_middle.transpose().triangularView<Upper>().solve(L_middle.triangularView<Lower>().solve(R_inv_y));
     return R_inv_y_updated;
}

