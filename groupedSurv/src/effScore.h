/*****************************/
/* Jiaxing Lin               */
/* 03-13-2018                */
/* efficient score for single*/ 
/* SNP                       */
/*****************************/

#include "helper_EF.h"
#include <cmath>
#include <math.h> 
static double effScore_NF_S(double beta, 
  const VectorXd &params, const VectorXd &g_tmp,
  const MatrixXd &xmatrix_tmp, const VectorXd &kivec_tmp, 
  const VectorXd &deltavec_tmp, int ntps, bool reScore, VectorXd &USNPs) {

  
  int nanCount = 0;
  for(int non0 = 0; non0 < g_tmp.size(); non0++)
    if(isnan(g_tmp(non0))) 
       nanCount++;
  
  Eigen::MatrixXd xmatrix1 = xmatrix_tmp; //(g_tmp.size()-nanCount, xmatrix_tmp.cols());
  Eigen::VectorXd kivec1 = kivec_tmp;//(g_tmp.size()-nanCount);
  Eigen::VectorXd deltavec1 = deltavec_tmp;//(g_tmp.size()-nanCount);
  Eigen::VectorXd g1 = g_tmp;//(g_tmp.rows()-nanCount); 
  
  int nai = 0; 
  for(int non0 = 0; non0 < g_tmp.size(); non0++){
    if(!isnan(g_tmp(non0))) 
    {
       g1(nai) = g_tmp(non0);
       xmatrix1.row(nai) = xmatrix_tmp.row(non0);
       kivec1(nai) = kivec_tmp(non0);
       deltavec1(nai) = deltavec_tmp(non0);
       nai++;
    }
  }
 
  Eigen::MatrixXd xmatrix = xmatrix1.topRows(g_tmp.size()-nanCount);
  Eigen::VectorXd kivec = kivec1.head(g_tmp.size()-nanCount);
  Eigen::VectorXd deltavec = deltavec1.head(g_tmp.size()-nanCount);
  Eigen::VectorXd g = g1.head(g_tmp.size()-nanCount); 

  // parse data
  int n = g.size();
  int n_with_NA = g_tmp.size();
  Eigen::VectorXd gib = g*beta;
  Eigen::VectorXd gammavec = params.head(ntps);
  Eigen::VectorXd thetavec(params.size() - ntps);
  thetavec << params.tail(params.size() - ntps);
  Eigen::VectorXd xitheta = xmatrix * thetavec;
  
  // get repeated components of derivatives 
  Eigen::VectorXd dfrac = deltaifrac_EF(gammavec, gib, xitheta, kivec, deltavec);
  Eigen::MatrixXd summat = sumKim1mat_EF(gammavec, gib, xitheta, kivec, deltavec);
  
	// Partial wrt beta for each pt
  Eigen::VectorXd Sbetai = (dfrac - summat.rowwise().sum());
  for(int i = 0; i < n; i++)
    Sbetai(i) = Sbetai(i) *g(i);
 
  // Partial wrt gamma for each pt_i x gamma_j
  Eigen::MatrixXd Sgammai = summat*-1;
  for(int i = 0; i < n; i++)
    if(kivec(i) - 1 < ntps )
      Sgammai(i, kivec(i)-1) = dfrac(i);
 
  // Partial wrt theta for each pt_i x theta_l
  Eigen::MatrixXd Sthetai = xmatrix;
  Eigen::VectorXd dmsTmp  = dfrac - summat.rowwise().sum();
  for(int i = 0; i < n; i++)
    for(int j = 0; j < xmatrix.cols(); j++)
      Sthetai(i,j) = dmsTmp(i) * xmatrix(i,j); 

  // E[S_beta,gamma] = Sum over i of S_beta,i * S_gamma,i / n
  Eigen::VectorXd Ebeta_gamma = Sbetai.transpose()*Sgammai/n;
   
  // E[S_beta,theta] = Sum over i of S_beta,i * S_theta,i / n
  Eigen::VectorXd Ebeta_theta = Sbetai.transpose()*Sthetai/n; 
  //cout << "Ebeta_theta:\n" << Ebeta_theta << endl;
  
  Eigen::VectorXd Cov(Ebeta_theta.size()+ Ebeta_gamma.size());
  Cov << Ebeta_gamma, Ebeta_theta;

  // S_gamma,i x S_gamma,j (list of n ntps x ntps matrices)
  Eigen::MatrixXd Egamma_gamma(Sgammai.cols(), Sgammai.cols());
  Egamma_gamma.fill(0);
  for(int i = 0; i < Sgammai.rows(); i++)
	Egamma_gamma += Sgammai.row(i).transpose()*Sgammai.row(i);
  Egamma_gamma = Egamma_gamma/n;
 
  // S_theta,i x S_gamma,i (list of n l x ntps matrices)
  Eigen::MatrixXd Etheta_gamma(Sthetai.cols(), Sgammai.cols());
  Etheta_gamma.fill(0);
  for(int i = 0; i < Sgammai.rows(); i++)
     Etheta_gamma += Sthetai.row(i).transpose()*Sgammai.row(i);
  Etheta_gamma = Etheta_gamma/n;

  // S_theta,i x S_theta,i (list of n l x l matrices)
  Eigen::MatrixXd Etheta_theta(Sthetai.cols(), Sthetai.cols());
  Etheta_theta.fill(0);
  for(int i = 0; i < Sthetai.rows(); i++)
    Etheta_theta += Sthetai.row(i).transpose()*Sthetai.row(i);
  Etheta_theta = Etheta_theta/n;
  //cout << "Etheta_theta:\n" << Etheta_theta <<endl;
  // Here we are taking average before we make the inverse of the information matrix. 
  // put Efficent score together
  Eigen::MatrixXd Var(params.size(), params.size());
  Var.fill(0);
  Var.block(0, 0, Sgammai.cols(), Sgammai.cols()) = Egamma_gamma;
  Var.block(0, Sgammai.cols(), Sgammai.cols(), Sthetai.cols()) =  Etheta_gamma.transpose();
  Var.block(Sgammai.cols(), 0, Sthetai.cols(), Sgammai.cols()) = Etheta_gamma;
  Var.block(Sgammai.cols(), Sgammai.cols(), Sthetai.cols(), Sthetai.cols()) = Etheta_theta;
 
  // compute score
  Eigen::MatrixXd Sgamma_theta(n, params.size());
  Sgamma_theta.fill(0);
  Sgamma_theta.block(0, 0, n, Sgammai.cols()) = Sgammai;
  Sgamma_theta.block(0, Sgammai.cols(), n, Sthetai.cols()) = Sthetai;
  Eigen::VectorXd Ui(n);
  Ui.fill(0);

  Eigen::MatrixXd VarInverse = Var.inverse();
  
  for(int i = 0; i < n; i++)
  { 
 	Ui(i) = Sbetai(i)  - Cov.transpose()*VarInverse*(Sgamma_theta.row(i).transpose());
    if(isinf(Ui(i)) || isnan(Ui(i)))
				Ui(i) = 0;
	} 
	// return all the scores
	if(reScore)
  {  
    int i_Ui = 0; 
    for(int i = 0; i < n_with_NA; i++)
    {
       if(isnan(g_tmp(i))){
         USNPs(i) = NAN;
       }else{
         USNPs(i) = Ui(i_Ui);
         i_Ui++;
       }
    }
  }
  
  // return statistics 
  double stat = Ui.sum()*Ui.sum()/(Ui.cwiseProduct(Ui)).sum();
  if(isnan(stat))
    stat = 0;
  return stat;
}
