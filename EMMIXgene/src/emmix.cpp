//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

//'@importFrom Rcpp sourceCpp
//'@useDynLib EMMIXgene


#include <vector>
#include <string>
#include <algorithm>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/roots.hpp>
#include <RcppArmadillo.h>
#include <tkmeans.h>

//using namespace arma;
using namespace Rcpp;
//using namespace boost::math;


typedef std::vector<double> stdvec;

double df_eq_est (double v, double v_sum);
//function class with operator for root finder for df of t-dist
class df_eq_func{
public:
  df_eq_func(double v) :  v_sum(v){};

  double getSumV();

  double operator()(double v) {

    double temp = -boost::math::digamma(0.5*v) + std::log(0.5*v) + v_sum +
        boost::math::digamma(0.5*(v+1.)) - std::log(0.5*(v+1.))+1.;
    return temp;
    //return df_eq_est(v,v_sum);

  }

private:
  double v_sum;
};

double df_eq_func::getSumV()
{
  return v_sum;
}


Rcpp::NumericVector export_vec(arma::vec y)
{
  Rcpp::NumericVector tmp = Rcpp::wrap(y);
  tmp.attr("dim") = R_NilValue;
  return tmp;
}

Rcpp::NumericVector export_uvec(arma::uvec y)
{
  Rcpp::NumericVector tmp = Rcpp::wrap(y);
  tmp.attr("dim") = R_NilValue;
  return tmp;
}

//this is only 1D for now

// [[Rcpp::export]]
double mahalanobis_c(double y, double mu, double sigma)
{
  double delta = (y-mu)*(1./sigma)*(y-mu);
  return delta;
}

//this is only 1D for now

// [[Rcpp::export]]
double t_dist(double y, double mu, double sigma, double nu,  int p =1)
{
    double pdf =((tgamma(0.5*(nu+p))/sqrt(sigma))/
        (std::pow(arma::datum::pi*nu, 0.5*p ) *
        tgamma(0.5*nu)*std::pow(1+(mahalanobis_c(y,mu,sigma)/nu), 0.5*(nu+p))));
    return pdf;
}



// [[Rcpp::export]]
arma::mat estep(const arma::vec& dat, const arma::mat& params){
  int n = dat.size();
  int g = params.n_rows;
  int g2 =2*g;
  
  arma::mat new_params(g2, n, arma::fill::zeros);
  arma::mat tau(g,n, arma::fill::zeros);
  arma::rowvec tau_sum(n, arma::fill::zeros);
  arma::mat you(g,n, arma::fill::zeros);

  arma::vec pi = params.col(0);
  arma::vec mu = params.col(1);
  arma::vec nu = params.col(2);
  arma::vec sigma = params.col(3);
  

  for(int i=0;i<g;i++){
    stdvec dens(n);
    double mu_i = mu.at(i);
    double sigma_i = sigma.at(i);
    double nu_i = nu.at(i);
    std::transform(dat.begin(), dat.end(), dens.begin(), 
                   [mu_i, sigma_i, nu_i](double dat)
                       { return t_dist(dat, mu_i, sigma_i, nu_i); });
    arma::vec dens2 = arma::conv_to<arma::vec>::from(dens);
    tau.row(i)  = pi.at(i)*dens2.t();
    you.row(i) = ((nu.at(i)+1.0)/(nu.at(i) + 
        (dat-mu.at(i))%(dat-mu.at(i))*(1/sigma.at(i)))).t();

  }
  
  
  arma::mat tau2 = tau.each_row() / sum(tau,0);
  
  
  
  //tau.print();
  new_params.rows(0,g-1)=tau2;
  new_params.rows(g,g2-1)=you;

  return(new_params);
}

arma::mat start_kmeans(arma::vec dat, int g){
  arma::mat params(g,4, arma::fill::zeros);
  arma::vec weights(1, arma::fill::ones);
  arma::uvec alloc(dat.n_elem, arma::fill::zeros);
  arma::uvec cluster_sizes(g, arma::fill::zeros);
  
  params.col(1) = tkmeans(dat, g, 0.0, weights, 1);
  
  for(int i=0; i<((int)dat.n_elem);i++){
    arma::vec temp = abs(params.col(1)- dat(i)); 
    alloc.at(i) = arma::index_min(temp);
  }
  
  //cluster size and variance
  for(int i =0; i<g;i++){
    arma::uvec tmp1 = find(alloc == i);
    cluster_sizes.at(i) = tmp1.n_elem;
    params(i,3) = sum(pow(dat(tmp1)-params(i,1),2.0))/dat.n_elem;
    params(i,2) = 4.0+(sum(pow((dat(tmp1)-params(i,1))/sqrt(params(i,3)),4.0))/
        dat.n_elem)/6.0;
    if(params(i,2) < 0 ){params(i,2) = 2.5;}
    if(params(i,2) > 200 ){params(i,2) = 199.99;}
  }
  
  params.col(0) = arma::conv_to< arma::vec >::from(cluster_sizes)/dat.n_elem;
  return(params);
}


arma::mat start_random(arma::vec dat, int g){
  arma::mat params(g,4, arma::fill::zeros);
  arma::vec weights(1, arma::fill::ones);
  arma::uvec alloc(dat.n_elem, arma::fill::zeros);
  arma::uvec cluster_sizes(g, arma::fill::zeros);
  
  
  arma::ivec rstarts = arma::randi(g, arma::distr_param(0,dat.n_elem-1));
  for (int i=0; i<g; i++){
    params(i,1) = dat(rstarts(i));
   
  }
  
  
  for(int i=0; i<((int)dat.n_elem);i++){
    arma::vec temp = abs(params.col(1)- dat(i)); 
    alloc.at(i) = arma::index_min(temp);
  }
  
  //cluster size and variance
  for(int i =0; i<g;i++){
    arma::uvec tmp1 = find(alloc == i);
    cluster_sizes.at(i) = tmp1.n_elem;
    params(i,3) = sum(pow(dat(tmp1)-params(i,1),2.0))/dat.n_elem;
    params(i,2) = 4.0+(sum(pow((dat(tmp1)-params(i,1))/sqrt(params(i,3)),4.0))/
        dat.n_elem)/6.0;

  }
  
  params.col(0) = arma::conv_to< arma::vec >::from(cluster_sizes)/dat.n_elem;
  return(params);
}



// [[Rcpp::export]]
arma::mat mstep(arma::vec dat, arma::mat tau, arma::mat you, arma::mat params){
  int n = dat.n_elem;
  int g = tau.n_rows;

  arma::vec mu2 = arma::zeros<arma::vec>(g);
  arma::vec nu2 = arma::zeros<arma::vec>(g);
  arma::vec nu = params.col(2);
  arma::vec sigma2 = arma::zeros<arma::vec>(g);
  arma::vec s_tau = arma::zeros<arma::vec>(g);
  arma::vec pi2 = arma::zeros<arma::vec>(g);
  
  s_tau = sum(tau,1);
  pi2 = s_tau/n;
  

  
  double a = 0.01;
  double b = 199.99;
  boost::uintmax_t df_max_iter=500; // boost solver params, 
  //could be user params but relatively unimportant
  boost::math::tools::eps_tolerance<double> tol(30);
 
  for(int i=0;i<g;i++){
    mu2(i) = (sum(dat.t()%you.row(i)%tau.row(i)) / sum(you.row(i)%tau.row(i)));
    sigma2(i) = (sum(you.row(i)%tau.row(i)%((dat.t()-mu2(i))%(dat.t()-mu2(i)))) 
                     / s_tau(i));
    
  
    
    double v_sum = (1/s_tau(i))*sum(tau.row(i)%
    (arma::log(you.row(i))-you.row(i)));
    // Rcpp::Rcout << "(s_tau(i)) = " << (s_tau(i)) << std::endl;
    // Rcpp::Rcout << "(1/s_tau(i)) = " << (1/s_tau(i)) << std::endl;
    // Rcpp::Rcout << "v_sum = " << v_sum << std::endl;
    if(!std::isfinite(v_sum)){v_sum =0.0;}
    
    df_eq_func rootFun = df_eq_func(v_sum);
      
    try{
        std::pair<double, double>  r1= boost::math::tools::bisect(rootFun, a, b, tol, df_max_iter);
        nu2(i) = (r1.first + r1.second)/2.0;
    }catch(...){
        nu2(i) = 199.99;
        //nu2(i) = nu(i);
    }
    //nu2.t().print();
    
  }
  //nu2.t().print();
  
  //nu2.fill(20);  
    
  arma::mat new_params(g,4, arma::fill::zeros);
  new_params.col(0) = pi2;
  new_params.col(1) = mu2;
  new_params.col(2) = nu2;
  new_params.col(3) = sigma2;

  return(new_params);
}


double BIC_calc(double LL,  int k, int n){
  double BIC = std::log((double)n)* k - 2*LL ;
  return(BIC);
}

double AIC_calc(double LL,  int k){
  double AIC = 2*k - 2*LL;
  return(AIC);
}



// [[Rcpp::export]]
List emmix_t(arma::vec dat, int g=1, int random_starts=4, int max_it=100,
             double tol = 0.0001, std::string start_method = "both"){
  int n = dat.size();
  int k =0; int n_fixed = 0;
  double  LL =0;
  arma::mat params(g,4, arma::fill::zeros);
  List return_list;
  
  arma::vec pi = arma::zeros(g);
  arma::vec mu = arma::zeros(g);
  arma::vec nu = arma::zeros(g);
  arma::vec sigma = arma::zeros(g);
 
  arma::uvec clusters = arma::zeros<arma::uvec>(n);
  arma::uvec cluster_sizes = arma::zeros<arma::uvec>(g);
  arma::vec lik = arma::zeros(max_it);
  
  arma::mat tau = arma::zeros(g,n);
  arma::mat you = arma::zeros(g,n);
  arma::mat temp2(g,4, arma::fill::zeros);
  
  //double Q1 = 0.0;
  arma::vec Q2 =  arma::zeros(g);
  arma::mat Q3 =  arma::zeros(g,n);
  arma::mat Q2b =  arma::zeros(g,n);
  
  std::vector<List> tempOut;
  
  int starts =0; int safety_count = 0;
  while(starts < random_starts){
    
        if(start_method == "kmeans"){
          params = start_kmeans(dat, g);
          //params.print();
        }
        
       if(start_method == "random"){
         params = start_random(dat, g);
         //params.print();
       }
       
       if(start_method == "both"){
         if(starts<(random_starts/2)){
          params = start_random(dat, g);
         }else{
           params = start_kmeans(dat, g);
         }
         //params.print();
       }

      //primary loop
        int j =0; double diff=10.0; double old_LL=0;
        
        while((j<max_it) && (diff > tol)){
          
          arma::mat temp1 = estep(dat,params);
          tau = temp1.rows(0,g-1);
          you = temp1.rows(g,2*g-1);
          temp2 = mstep(dat, tau, you, params);
          params = temp2; 
          

          //log liklihood
          double accumQ1=0;
          double accumQ2=0;
          double accumQ3=0;
          
          mu=params.col(1);
          nu=params.col(2);
          sigma = params.col(3);

          for(int i=0; i<g; i++){
              Q2(i) = -std::log(boost::math::tgamma(0.5*(nu.at(i))))  + 0.5*nu.at(i)*std::log(0.5*nu.at(i)) - 0.5*nu.at(i)*(boost::math::digamma(0.5*(nu.at(i)+1.0))-std::log(0.5*(nu.at(i)+1.0)) + arma::sum(arma::log(you.row(i))-you.row(i)));
              Q3.row(i) = -0.5*1.0*std::log(2.0*arma::datum::pi) - 0.5*std::log(std::abs(sigma.at(i))) + 0.5*1.0*arma::log(you.row(i)) - 0.5*you.row(i)*(1.0/sigma.at(i))%((dat-mu.at(i))%(dat-mu.at(i))).t();
          }
          
          for(int i=0; i <g; i++){
              for(int k=0; k <n; k++){
                  accumQ1+=tau(i,k)*std::log(params(i,0));
                  accumQ2+=tau(i,k)*Q2(i);
                  accumQ3+=tau(i,k)*Q3(i,k);
              }
          }

          lik.at(j) = accumQ1 + accumQ2 + accumQ3 ;
          LL = lik.at(j);
          
          diff = std::abs(LL - old_LL);
          old_LL = LL;
          
          
          j++;
          
        }
        
   
    
    
    
    //cluster allocation
    for(int i =0; i<n;i++){
      clusters.at(i) = tau.col(i).index_max();
    }
    
    //cluster size
    for(int i =0; i<g;i++){
      arma::uvec tmp1 = find(clusters == i);
      cluster_sizes.at(i) = tmp1.n_elem;
    }
    //min cluster size
    int c_min = min(cluster_sizes);
    arma::vec lik2 = lik.rows(0,j-1);
    
    
    //n_params
    k = (g-1)+ g + g + g - n_fixed;//pi + mu + sigma + nu - fixed params
    if(std::isfinite(LL)){
      List temp_list = List::create(
        Named("mu")= export_vec(params.col(1)),
        Named("pi") = export_vec(params.col(0)),
        Named("nu") = export_vec(params.col(2)),
        Named("sigma") = export_vec(params.col(3)),
        Named("LL") = LL ,
        Named("BIC") = BIC_calc(LL, k, n),
        Named("AIC") = AIC_calc(LL, k),
        Named("tau") = tau,
        Named("lik") = lik2,
        Named("k") = k, //number of parameters
        Named("Clusters") = export_uvec(clusters),
        Named("c_min") = c_min, //smallest cluster,
        Named("n_iter") = j,//total iterations, if not equal max_iter then 
        //terminated due to meeting tolerance
        Named("components") = g,
        Named("Ratio") = 0.0
      );
      
      tempOut.push_back(temp_list);
      starts++;
    }else{
      if(safety_count>random_starts*10){
          //return_list = tempOut.back();
        List return_list = List::create(
          Named("LL") = - arma::datum::inf,
          Named("c_min") = 0
        );
        return return_list;
      }
    }
    safety_count++;
  }
    
  arma::vec temp_ll = arma::zeros(random_starts);
  
  for(int i=0; i<random_starts;i++){
    temp_ll(i)= as<double>(tempOut.at(i)["LL"]);
  }
  
  return_list = tempOut.at(temp_ll(find_finite(temp_ll)).index_max());
  
  return return_list;

}

// [[Rcpp::export]]
List each_gene(arma::vec dat, int random_starts=4, int max_it = 100,
               double ll_thresh = 8, int min_clust_size = 8,
               double tol = 0.0001, std::string start_method = "both",
               bool three=false){
  
  List g1 = emmix_t(dat, 1, random_starts, max_it, tol, start_method);
  List g2 = emmix_t(dat, 2, random_starts, max_it, tol, start_method);
  

  
  List best_g = g1;
 
  double lambda = 0;
  double lambda_s = 0;// ll_ratio(List g1, List g2)
  
  double ll_g1=as<double>(g1["LL"]);
  double ll_g2=as<double>(g2["LL"]);
  
  lambda = ll_g1-ll_g2;
  lambda_s = ll_g1-ll_g2;
  
  
  

  if(-2*(lambda) > ll_thresh){
    if(as<int>(g2["c_min"]) >= min_clust_size){
      best_g = g2;
      lambda_s = lambda;
    }else{
      
      if(three==true){
        List g3 = emmix_t(dat, 3, random_starts, max_it, tol, start_method);
        
        double ll_g3=as<double>(g3["LL"]);
        
        lambda = ll_g2-ll_g3;

        if(-2*(lambda) > ll_thresh){
          if(as<int>(g3["c_min"]) >= min_clust_size){
            best_g = g3;
            lambda_s = lambda;
          }
        }
      }
      
    }
  }
  
  
  best_g["Ratio"] = -2*(lambda_s);
 
  
  
  return(best_g);
}



// [[Rcpp::export]]
List emmix_gene(arma::mat& bigdat, int random_starts=4, int max_it = 100,
                double ll_thresh = 8, int min_clust_size = 8, 
                double tol = 0.0001, std::string start_method = "both",
                bool three = false){
 
 int n = bigdat.n_rows;
 Rcpp::List tmp;
 
 arma::vec stat =  arma::zeros(n);
 arma::vec comp =  arma::zeros(n);
 arma::vec it =  arma::zeros(n);
 
 for(int i=0; i<n;i++){

   tmp = (each_gene(bigdat.row(i).t(), random_starts, max_it, ll_thresh,
                    min_clust_size, tol, start_method));
   
   stat.at(i) = as<double>(tmp["Ratio"]);
   comp.at(i) = as<double>(tmp["components"]);
   it.at(i) = as<double>(tmp["n_iter"]);
   
   //Rcpp::Rcout << i<< std::endl;
  
   
 }
 
 
 Rcpp::List ret =List::create(
   Named("stat")= export_vec(stat),
   Named("g")= export_vec(comp),
   Named("it")= export_vec(it)
 );
 
 return(ret);
}
