//
//  RMSD.cpp
//  
//
//  Created by Dylan Shi on 2018-09-28.
//

// #include "RMSD.hpp"
#include <iostream>
#include <stdio.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
using namespace Rcpp;
using namespace RcppEigen;
// using Eigen::MatrixXd;
// using Eigen::VectorXd;

// [[Rcpp::depends(RcppEigen)]]

// Calculates squared euclidean distance between two vectors 
double dis(Eigen::VectorXd x, Eigen::VectorXd y){
    const Eigen::VectorXd diff = x - y;
    return diff.dot(diff);
}

// Returns Matrix x after subtracting its centroid
Eigen::MatrixXd center(Eigen::MatrixXd x){
    int num_cols = x.cols();
    Eigen::ArrayXd row_x = x.row(0);
    Eigen::ArrayXd row_y = x.row(1);
    Eigen::ArrayXd row_z = x.row(2);
    double x_c = row_x.mean();
    double y_c = row_y.mean();
    double z_c = row_z.mean();
    Eigen::MatrixXd centered(3, num_cols);
    centered.row(0) = row_x-x_c;
    centered.row(1) = row_y-y_c;
    centered.row(2) = row_z-z_c;
    return centered;
}

// Calculates the RMSD between two matricies x, y
// [[Rcpp::export]]
double RMSD(Eigen::MatrixXd x, Eigen::MatrixXd y){
    // Eigen::MatrixXd x = center(xr);
    // Eigen::MatrixXd y = center(yr);
    double size = x.cols();
    double sum = 0;
    for (int i=0; i<size; i++){
        sum += dis(x.col(i), y.col(i));
    }
    return sqrt(sum/size);
}

// [[Rcpp::export]]
List LRMSD(Eigen::MatrixXd xr, Eigen::MatrixXd yr){
  Eigen::MatrixXd x = center(xr);
  Eigen::MatrixXd y = center(yr);
  Eigen::MatrixXd C = x*y.adjoint();
  Eigen::MatrixXd T = Eigen::MatrixXd::Zero(3, 3);
  T(0,0) = 1;
  T(1,1) = 1;
  T(2,2) = (C.determinant()>=0) ? 1 : -1;
  Eigen::BDCSVD<Eigen::MatrixXd> const svd(C , Eigen::ComputeFullU | Eigen::ComputeFullV);
  const Eigen::MatrixXd U = svd.matrixU();
  const Eigen::MatrixXd V = svd.matrixV();
  const Eigen::MatrixXd s = svd.singularValues();
  const Eigen::MatrixXd rotation = V*T*U.adjoint();
  // std::cout << "Rotation: " << rotation << std::endl;
  List L = List::create(_["RMSD"]=RMSD(rotation * x, y), 
                        _["rotation"]=rotation);
return L;
}

