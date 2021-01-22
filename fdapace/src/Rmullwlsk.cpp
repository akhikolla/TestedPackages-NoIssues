#include <RcppEigen.h>
#include <map>          // to map kernels to integers for the switch
#include <string>       // to read in the kernel name
#include <vector>       // to use vectors
#include <algorithm>    // to get the intersect and sort

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]


Eigen::MatrixXd Rmullwlsk( const Eigen::Map<Eigen::VectorXd> & bw, const std::string kernel_type, const Eigen::Map<Eigen::MatrixXd> & tPairs, const Eigen::Map<Eigen::MatrixXd> & cxxn, const Eigen::Map<Eigen::VectorXd> & win,  const Eigen::Map<Eigen::VectorXd> & xgrid, const Eigen::Map<Eigen::VectorXd> & ygrid, const bool & bwCheck, const bool & transp = true){ 

  // tPairs : xin (in MATLAB code)
  // cxxn : yin (in MATLAB code)
  // xgrid: out1 (in MATLAB code)
  // ygrid: out2 (in MATLAB code)
  // bwCheck : boolean/ cause the function to simply run the bandwidth check.

  const double invSqrt2pi=  1./(sqrt(2.*M_PI));

  // Map the kernel name so we can use switches  
  std::map<std::string,int> possibleKernels; 
  possibleKernels["epan"]    = 1;   possibleKernels["rect"]    = 2;
  possibleKernels["gauss"]   = 3;   possibleKernels["gausvar"] = 4; 
  possibleKernels["quar"]    = 5; 
   
  // The following test is here for completeness, we mightwant to move it up a 
  // level (in the wrapper) in the future. 

  // If the kernel_type key exists set KernelName appropriately
  int KernelName = 0;
  if ( possibleKernels.count( kernel_type ) != 0){ 
    KernelName = possibleKernels.find( kernel_type )->second; //Set kernel choice
  } else {
  // otherwise use "epan"as the kernel_type 
    //Rcpp::Rcout << "Kernel_type argument was not set correctly; Epanechnikov kernel used." << std::endl;
    Rcpp::warning("Kernel_type argument was not set correctly; Epanechnikov kernel used.");
    KernelName = possibleKernels.find( "epan" )->second;;
  }

  // Check that we do not have zero weights // Should do a try-catch here
  // Again this might be best moved a level-up. 
  if ( !(win.all()) ){  // 
    Rcpp::Rcout << "Cases with zero-valued windows are not yet implemented" << std::endl;
    return (tPairs);
  } 

  // Start the actual smoother here  
  unsigned int xgridN = xgrid.size();  
  unsigned int ygridN = ygrid.size();  
  
  Eigen::MatrixXd mu(ygrid.size(),xgrid.size());
  mu.setZero();    
  const double bufSmall = 1.0e-6; // pow(double(10),-6);

  for (unsigned int i = 0; i != ygridN; ++i){  
    for (unsigned int j = i; j < xgridN ; ++j){ 

      //locating local window (LOL) (bad joke)
      std::vector <unsigned int> indx; 
      //if the kernel is not Gaussian
      if ( KernelName != 3) { 
        //construct listX as vectors / size is unknown originally
        for (unsigned int y = 0; y != tPairs.cols(); y++){ 
          if ( std::abs( tPairs(0,y) - xgrid(j) ) <= (bw(0)+ bufSmall) && std::abs( tPairs(1,y) - ygrid(i) ) <= (bw(1)+ bufSmall) ) {
          // legacy MATLAB equivalent form :  
          // if ( (tPairs(0,y) >= (xgrid(j) -(bw(0)+ pow(double(10),-6)))) & (tPairs(0,y) <= (xgrid(j) + (bw(0)+  pow(double(10),-6))))) {
            indx.push_back(y);
          }         
        }
      } else{ // just get the whole deal
        for (unsigned int y = 0; y != tPairs.cols(); ++y){
          indx.push_back(y);
        }
      }  

      unsigned int indxSize = indx.size();
      Eigen::VectorXd lw(indxSize);  
      Eigen::VectorXd ly(indxSize);
      Eigen::MatrixXd lx(2,indxSize);

      for (unsigned int u = 0; u !=indxSize; ++u){ 
        lx.col(u) = tPairs.col(indx[u]); 
        lw(u) = win(indx[u]); 
        ly(u) = cxxn(indx[u]); 
      }

      // check enough points are in the local window 
      unsigned int meter=1;  
      for (unsigned int u =0; u< indxSize; ++u) { 
        for (unsigned int t = u + 1; t < indxSize; ++t) {
          if ( (lx(0,u) !=  lx(0,t) ) || (lx(1,u) != lx(1,t) ) ) {
            meter++;
          }
        }
        if (meter >= 3) { 
          break; 
        }
      }
   
      //computing weight matrix 
      if (meter >=  3 && !bwCheck) { 
        Eigen::VectorXd temp(indxSize);
        Eigen::MatrixXd llx(2, indxSize );  
        llx.row(0) = (lx.row(0).array() - xgrid(j))/bw(0);  
        llx.row(1) = (lx.row(1).array() - ygrid(i))/bw(1); 

        //define the kernel used 

        switch (KernelName){
          case 1: // Epan
            temp=  ((1-llx.row(0).array().pow(2))*(1- llx.row(1).array().pow(2))).array() * 
                   ((9./16)*lw).transpose().array(); 
            break;  
          case 2 : // Rect
            temp=(lw.array())*.25 ; 
            break;
          case 3 : // Gauss
            temp = ((-.5*(llx.row(1).array().pow(2))).exp()) * invSqrt2pi  *   
                   ((-.5*(llx.row(0).array().pow(2))).exp()) * invSqrt2pi  *
                   (lw.transpose().array()); 
            break;
          case 4 : // GausVar
            temp = (lw.transpose().array()) * 
		   ((-0.5 * llx.row(0).array().pow(2)).array().exp() * invSqrt2pi).array() *
                   ((-0.5 * llx.row(1).array().pow(2)).array().exp() * invSqrt2pi).array() * 
                   (1.25 - (0.25 * (llx.row(0).array().pow(2))).array())  * 
                   (1.50 - (0.50 * (llx.row(1).array().pow(2))).array()); 
            break;
          case 5 :  // Quar
              temp = (lw.transpose().array()) * 
                     ((1.-llx.row(0).array().pow(2)).array().pow(2)).array() *
                     ((1.-llx.row(1).array().pow(2)).array().pow(2)).array() * (225./256.);
            break;
        } 

        // make the design matrix
        Eigen::MatrixXd X(indxSize ,3);
        X.setOnes();    
        X.col(1) = lx.row(0).array() - xgrid(j);
        X.col(2) = lx.row(1).array() - ygrid(i); 
        Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWX(X.transpose() * temp.asDiagonal() *X);
        // The solver should stop if the value is NaN. See the HOLE example in gcvlwls2dV2.
        Eigen::VectorXd beta = ldlt_XTWX.solve(X.transpose() * temp.asDiagonal() * ly);
        mu(i,j)=beta(0);  
      } else if(meter < 3){
        // Rcpp::Rcout <<"The meter value is:" << meter << std::endl;  
        if (bwCheck) {
            Eigen::MatrixXd checker(1,1);
            checker(0,0) = 0.;
            return(checker);
        } else {
            Rcpp::stop("No enough points in local window, please increase bandwidth.");
        }
      }
    }
  }

  if (bwCheck){
     Eigen::MatrixXd checker(1,1); 
     checker(0,0) = 1.; 
     return(checker);
  }
    
  if (transp) {
  //  Eigen::MatrixXd result(ygrid.size(),xgrid.size());
  //  result = mu + mu.transpose();
  //  result.diagonal() = 0.5 * result.diagonal();
  //  return ( result ); 
    return ( Eigen::MatrixXd(mu.triangularView<Eigen::StrictlyUpper>().transpose()) + Eigen::MatrixXd(mu.triangularView<Eigen::Upper>()));
  } else {
    return( mu );
  }
}

 
