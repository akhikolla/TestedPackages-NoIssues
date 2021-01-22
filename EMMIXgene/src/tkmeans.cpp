//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <tkmeans.h>
// constructing priority queues
#include <iostream>
#include <queue>
#include <vector>
#include <functional>
#include <tuple>
// maths
#include <math.h>


const double log2pi = std::log(2.0 * M_PI);

typedef std::tuple<double, int, int>  di_pair;

class CompareDist
{
public:
  bool operator()(di_pair n1, di_pair n2) {
    return std::get<0>(n1)<std::get<0>(n2);
  }
};


arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov) {
  unsigned int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (unsigned int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return arma::sum((x_cen * cov.i()) % x_cen, 1);
}



arma::vec dmvnorm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool log_flag = false) {
  arma::vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = arma::sum(arma::log(arma::eig_sym(sigma)));
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;

  if (log_flag) {
    return(logretval);
  } else {
    return(exp(logretval));
  }
}



arma::vec distCentre(int k, arma::mat ui, arma::mat centres, double lambda, int d){
  //set up sigma
  arma::vec sigmaDiag(d);  sigmaDiag.fill(lambda);
  arma::mat sigma = diagmat(sigmaDiag);

   //store distances
  arma::vec dist =  arma::zeros(k);

  for (int j=0; j<k; j++)
  {
    dist.at(j) = 1.0/k * dmvnorm_arma(ui, centres.row(j), sigma, false)[0];
  }
  //Rcpp::Rcout << "size: " << size(dist) << std::endl;
  return(dist);
}


arma::vec distCentre2(int k, arma::mat ui, arma::mat centres, double lambda, int d){

  //store distances
  arma::vec dist =  arma::zeros(k);

  for (int j=0; j<k; j++)
  {
    dist.at(j) = std::sqrt(arma::sum(arma::pow(centres.row(j) - ui, 2.0)));
  }


  return(dist);
}


arma::vec distCentre3(int k, arma::mat ui, arma::mat centres, double lambda, int d){
  arma::vec dist = arma::sqrt(arma::sum(arma::pow(centres - repmat(ui,k,1), 2.0),1));
  return(dist);
}



int max_index(arma::vec x){

  //returns index of biggest value in x.
  Rcpp::NumericVector y = NumericVector(x.begin(),x.end());
  Rcpp::NumericVector::iterator it =  std::max_element(y.begin(), y.end());

  //Rcpp::Rcout << "it  " << it <<  " y " << y.begin() <<  " diff " << std::distance(y.begin(), it) << std::endl;
  return (std::distance(y.begin(), it));
}



int min_index(arma::vec x){

  //returns index of biggest value in x.
  Rcpp::NumericVector y = NumericVector(x.begin(),x.end());
  Rcpp::NumericVector::iterator it =  std::min_element(y.begin(), y.end());

  //Rcpp::Rcout << "it  " << it <<  " y " << y.begin() <<  " diff " << std::distance(y.begin(), it) << std::endl;
  return (std::distance(y.begin(), it));
}



arma::mat init_centres(arma::mat M, int k, bool verbose){
  int n = size(M)[0];
  int d = size(M)[1];

  // if(verbose){
  //   Rcpp::Rcout << "Randomly initialising centres..." << std::endl;
  // }

  //randomly pick centres
  arma::mat B = arma::randi<arma::mat>(k,1,arma::distr_param(0,n-1));

  arma::mat cent =  arma::zeros(k,d);

  for(int i=0; i<k; i++){
    //Rcpp::Rcout << "i = " << i << " n = " << B.at(i,0) <<  std::endl;

    cent.row(i) = M.row(B.at(i,0));
  }

  return cent;
}

// [[Rcpp::export]]
arma::mat tkmeans(arma::mat& M, int k , double alpha, arma::vec weights,  int nstart, int iter , double tol, bool verbose){
  unsigned int n = size(M)[0];
  unsigned int d = size(M)[1];

  if(weights.n_elem!=d){
    stop("weights vector must be same length as number of columns in M.");
  }
// 
//   if((alpha>=1) | (alpha <=0)){
//     stop("alpha must be in (0,1).");
//   }

  if(nstart<1){
    stop("Number of starts must be at least 1.");
  }

  if(iter<5){
    stop("Maximum number of iterations is too low.");
  }

  if(verbose){
    Rcpp::Rcout << "n = " << n << " d = " << d <<  " k = " << k << std::endl;
    Rcpp::Rcout << "n starts = " << nstart << " max iterations " << iter <<  " alpha = " << alpha << std::endl;
  }

  //apply weights
  M.each_row() /= weights.t();

  arma::mat best_means = init_centres(M, k, verbose);
  //double temp_BIC = 0.0;
  //double best_BIC = 0.0;
  //int best_j = 0;

  for( int j=0; j<nstart;j++){

    if(verbose){
      Rcpp::Rcout << j+1 << " of " << nstart<< " starts." << std::endl;
    }

    arma::mat centres = init_centres(M, k, verbose);
    arma::mat means = centres;

    double lambda = 1;
    //double tol = 0.001;
    double diff = 1;
    int m=0;

    // cluster membership record
    //std::vector<int> cluster_membership;
    //cluster_membership.reserve(n);

    //keep list of smallest using priority queue
    std::priority_queue< di_pair, std::vector<di_pair>, CompareDist  >  smallest;

    unsigned int queue_len;
    if(alpha<=0.5){




      if(alpha != 0.0){
        queue_len = floor((n-1)*(alpha))+1;
      }else{
        queue_len = 0;
      }

      if(verbose){
        Rcpp::Rcout << "Removing " << alpha << " of points. Discard que length " <<   queue_len << std::endl;
      }


      //Rcpp::Rcout << "number of outliers = " <<  queue_len   << std::endl;
      //Rcpp::Rcout << "Begin main loop..."  << std::endl;



      while((m<iter) & (diff > tol)){


        //Rcpp::Rcout << m << " of " << iter  << " iterations" << std::endl;

        arma::mat new_centres = arma::zeros(k, d);
        arma::mat centre_members = arma::zeros(k, 1);

        for(unsigned int i=0; i<n; i++)
        {
          arma::mat Dc = distCentre2(k, M.row(i), centres, lambda, d);
          //distCentre(k, M.row(i), centres, lambda, d).print();

          unsigned int cluster = min_index(Dc);

          di_pair temp = di_pair(Dc.at(cluster), cluster, i);

          new_centres.row(cluster) =  new_centres.row(cluster)+ M.row(i);
          centre_members(cluster,0) =  centre_members(cluster,0) + 1;

          if(queue_len > 0){
            if(smallest.size() < queue_len){
              smallest.push(temp);
            }else{
              if(std::get<0>(smallest.top()) > std::get<0>(temp)){
                //Rcpp::Rcout << "out: "   <<  std::get<0>(smallest.top()) <<  " in: " << std::get<0>(temp) << std::endl;
                smallest.pop();
                smallest.push(temp);
              }
            }
          }
          //Rcpp::Rcout << "size of heap:"  <<  smallest.size() << std::endl;
        }


        //remove all the outliers
        if(queue_len > 0){

          while (!smallest.empty())
          {
            new_centres.row(std::get<1>(smallest.top())) = new_centres.row(std::get<1>(smallest.top())) - M.row(std::get<2>(smallest.top()));

            centre_members(std::get<1>(smallest.top()),0) = centre_members(std::get<1>(smallest.top()),0) - 1;

            //Rcpp::Rcout << " " << std::get<0>(smallest.top()) <<  " " << std::get<1>(smallest.top()) << " " << std::get<2>(smallest.top()) << std::endl;

            smallest.pop();
          }
        }

        if(centre_members.min() > 0) {
          means = new_centres.each_col() / centre_members;
        }else{
          if(verbose){
              Rcpp::Rcout << "Empty cluster, resetting centres with " << iter-m << " iterations remaining." << std::endl;
            }
          means = init_centres(M, k, verbose);
        }

        //centre_members.print();
        //Rcpp::Rcout << diff << std::endl;

        diff = arma::accu(arma::abs(centres-means));
        centres = means;
        m++;
        if(verbose){
          Rcpp::Rcout << "iteration  " << m << " of " << iter <<  ". diff > tol :  " << diff << " > " << tol << std::endl;
        }
      }
    }
    else{
      if(alpha != 0.0){
        queue_len = floor((n-1)*(1-alpha))+1;
      }else{
        queue_len = 0;
      }


      //Rcpp::Rcout << "number of outliers = " <<  queue_len   << std::endl;
      //Rcpp::Rcout << "Begin main loop..."  << std::endl;


      while((m<iter) & (diff > tol)){

        //Rcpp::Rcout << m << " of " << iter  << " iterations" << std::endl;

        arma::mat new_centres = arma::zeros(k, d);
        arma::mat centre_members = arma::zeros(k, 1);

        for(unsigned int i=0; i<n; i++)
        {
            arma::mat Dc = distCentre2(k, M.row(i), centres, lambda, d);
            //distCentre(k, M.row(i), centres, lambda, d).print();

            unsigned int cluster = min_index(Dc);

            di_pair temp = di_pair(Dc.at(cluster), cluster, i);

            //new_centres.row(cluster) =  new_centres.row(cluster)+ M.row(i);
            //centre_members(cluster,0) =  centre_members(cluster,0) + 1;

            if(queue_len > 0){
              if(smallest.size() < queue_len){
                 smallest.push(temp);
              }else{
                if(std::get<0>(smallest.top()) > std::get<0>(temp)){
                  //Rcpp::Rcout << "out: "   <<  std::get<0>(smallest.top()) <<  " in: " << std::get<0>(temp) << std::endl;
                  smallest.pop();
                  smallest.push(temp);
                }
              }
            }
            //Rcpp::Rcout << "size of heap:"  <<  smallest.size() << std::endl;
        }


      //remove all the outliers
        if(queue_len > 0){
            while (!smallest.empty())
            {
              new_centres.row(std::get<1>(smallest.top())) = new_centres.row(std::get<1>(smallest.top())) + M.row(std::get<2>(smallest.top()));

              centre_members(std::get<1>(smallest.top()),0) = centre_members(std::get<1>(smallest.top()),0) + 1;

              //Rcpp::Rcout << " " << std::get<0>(smallest.top()) <<  " " << std::get<1>(smallest.top()) << " " << std::get<2>(smallest.top()) << std::endl;

              smallest.pop();
            }
        }

        if(centre_members.min() > 0) {
            means = new_centres.each_col() / centre_members;
        }else{
          if(verbose){
            Rcpp::Rcout << "Empty cluster, resetting centres with " << iter-m << " iterations remaining." << std::endl;
          }
          means = init_centres(M, k, verbose);
        }

        //centre_members.print();
        //Rcpp::Rcout << diff << std::endl;

        diff = arma::accu(arma::abs(centres-means));
        centres = means;
        m++;

        if(verbose){
          Rcpp::Rcout << "iteration  " << m << " of " << iter <<  ". diff > tol :  " << diff << " > " << tol << std::endl;
        }
      }
    }

    

  }


// return to original because in place edited
  M.each_row() %= weights.t();
  return best_means;
}

// [[Rcpp::export]]
arma::uvec nearest_cluster(arma::mat& data, arma::mat& centres){
  unsigned int k =  centres.n_rows;
  double lambda = 1.0;
  unsigned int d = size(data)[0];
  
  if(data.n_cols!=centres.n_cols){
    stop("Cluster centre dimensionality does not match data.");
  }
  
  arma::uvec clusters = arma::zeros<arma::uvec>(d);
  
  for(unsigned int i=0;i<d;i++){
    arma::mat Dc = distCentre2(k, data.row(i), centres, lambda, d);
    unsigned int cluster = min_index(Dc);
    clusters.at(i) = cluster+1;//going back to 1 indexed R
  }
  
  return clusters;
}
