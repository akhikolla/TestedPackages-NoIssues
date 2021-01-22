#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "ParallelFunctions.h"
using namespace Rcpp;
using namespace RcppArmadillo;
using namespace RcppParallel;

// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate k-th order moment of normal distribution
//' @description This function recursively calculates k-th order moment of 
//' normal distribution.
//' @param k non-negative integer moment order.
//' @param mean numeric expected value.
//' @param sd positive numeric standard deviation.
//' @param return_all_moments logical; if \code{TRUE}, function returns 
//' (k+1)-dimensional numeric vector of moments of normally distributed random 
//' variable with mean = \code{mean} and standard deviation = \code{sd}. 
//' Note that i-th vector's component value corresponds to the (i-1)-th moment.
//' @template is_validation_Template
//' @param is_central logical; if \code{TRUE}, then central moments 
//' will be calculated.
//' @details This function estimates \code{k}-th order moment of normal 
//' distribution which mean equals to \code{mean} and standard deviation 
//' equals to \code{sd}.\cr
//' @template k_integer_Template
//' @template diff_type_Template
//' @return This function returns \code{k}-th order moment of
//' normal distribution which mean equals to \code{mean} and standard deviation 
//' is \code{sd}. If \code{return_all_moments} is \code{TRUE} then see this 
//' argument description above for output details.
//' @examples
//' ## Calculate 5-th order moment of normal random variable which
//' ## mean equals to 3 and standard deviation is 5.
//'
//' # 5-th moment
//' normalMoment(k = 5, mean = 3, sd = 5)
//' 
//' # (0-5)-th moments
//' normalMoment(k = 5, mean = 3, sd = 5, return_all_moments = TRUE)
//' 
//' # 5-th moment derivative respect to mean
//' normalMoment(k = 5, mean = 3, sd = 5, diff_type = "mean")
//' 
//' # 5-th moment derivative respect to sd
//' normalMoment(k = 5, mean = 3, sd = 5, diff_type = "sd")
//'
//' @export
// [[Rcpp::export]]
NumericVector normalMoment(int k = 0,
						   double mean = 0, 
						   double sd = 1,
						   bool return_all_moments = false,
						   bool is_validation = true,
						   bool is_central = false,
						   String diff_type = "NO")
{
	// Validation
	// ------------------------------------------------------------
	if (is_validation)
	{
		if (k < 0)
		{
			stop("parameter k should be non-negative integer");
		}
		if (sd <= 0)
		{
			stop("parameter sd should be positive integer");
		}
		if((diff_type != "NO") & (diff_type != "mean") & (diff_type != "sd"))
		{
		  stop("diff_type argument should take value 'NO', 'mean' or 'sd'");
		}
	}
	// ------------------------------------------------------------
	
	// Estimate unconditional (on truncation) variance
	double sd_squared = pow(sd, 2);
	
	// Initialize matrix to store the moments
	NumericVector moments(k + 1, 1.0);
	
	// Initialize matrix to store the moments derivatives if need
	NumericVector moments_diff(k + 1, 1.0);
	
	// Check weather central moment should be calculated
	if(is_central)
	{
	  mean = 0;
	}
	
	// The zero moment always equals 1 and its derivatives respect to
	// mean and sd parameters are 0
	moments[0] = 1;
	moments_diff[0] = 0;

	if (k == 0)
	{
	  if(diff_type != "NO")
	  {
	    return(moments_diff);
	  }
		return(moments);
	}
	
	// If the moment is 1 it equals to mean, its derivative respect to
	// mean equals 1 and respect to sd equals zero
	moments[1] = mean;
	if(diff_type == "mean")
  {
    moments_diff[1] = 1;
  }
  if(diff_type == "sd")
  {
    moments_diff[1] = 0;
  }

	if (k == 1)
	{
		if (!return_all_moments)
		{
		  if(diff_type != "NO")
		  {
		    NumericVector moments_diff_return = NumericVector::create(moments_diff[k]);
		    return(moments_diff_return);
		  }
		  NumericVector moments_return = NumericVector::create(moments[k]);
		  return(moments_return);
		}
		return(moments);
	}

	// Recursively calculate other moments or derivatives

	  // calculate moments itself
	for (int i = 2; i <= k; i++)
	{
		moments[i] = (i - 1) * sd_squared * moments[i - 2] + mean * moments[i - 1];
	}
	
	  // calculate derivative respect to mean if need
	if(diff_type == "mean")
	{
	  for (int i = 2; i <= k; i++)
	  {
	    moments_diff[i] = (i - 1) * sd_squared * moments_diff[i - 2] + 
	    mean * moments_diff[i - 1] + moments[i - 1];
	  }
	 }
	  
	  // calculate derivative respect to mean if need
	 if(diff_type == "sd")
	 {
	   for (int i = 2; i <= k; i++)
	   {
	     moments_diff[i] = (i - 1) * (2 * sd) * moments[i - 2] + 
	     (i - 1) * sd_squared * moments_diff[i - 2] +
	     mean * moments_diff[i - 1];
	   }
	 }

	// Return depends on return_all_moments value
	
	  // if all moments should be returned
	if (!return_all_moments)
	{
	  if(diff_type != "NO")
	  {
	    NumericVector moments_diff_return = NumericVector::create(moments_diff[k]);
	    return(moments_diff_return);
	  }
		NumericVector moments_return = NumericVector::create(moments[k]);
		return(moments_return);
	}

	  // if only k-th moment should be returned
	if(diff_type != "NO")
	{
	  return(moments_diff);
	}
	return(moments);
}

//' Calculate k-th order moment of truncated normal distribution
//' @description This function recursively calculates k-th order moment of 
//' truncated normal distribution.
//' @param k non-negative integer moment order.
//' @param x_lower numeric vector of lower truncation points.
//' @param x_upper numeric vector of upper truncation points.
//' @param mean numeric expected value.
//' @param sd positive numeric standard deviation.
//' @template pdf_lower_Template
//' @template cdf_lower_Template
//' @template pdf_upper_Template
//' @template cdf_upper_Template
//' @template cdf_difference_Template
//' @template is_validation_Template
//' @template is_parallel_Template
//' @template diff_type_Template
//' @param return_all_moments logical; if \code{TRUE}, function returns the 
//' matrix of moments of normally distributed random variable with 
//' mean = \code{mean} and standard deviation = \code{sd} under lower and upper 
//' truncation points \code{x_lower} and \code{x_upper} correspondingly. 
//' Note that element in i-th row and j-th column of this matrix corresponds to 
//' the i-th observation (j-1)-th order moment.
//' @details This function estimates \code{k}-th order moment of
//' normal distribution which mean equals to \code{mean} and standard deviation 
//' equals to \code{sd} truncated at points given by \code{x_lower} and 
//' \code{x_upper}. Note that the function is vectorized so you can provide
//' \code{x_lower} and \code{x_upper} as vectors of equal size. If vectors values 
//' for \code{x_lower} and \code{x_upper} are not provided then their default 
//' values will be set to \code{-(.Machine$double.xmin * 0.99)} and 
//' \code{(.Machine$double.xmax * 0.99)} correspondingly.
//' @template k_integer_Template
//' @template pdf_cdf_precalculated_Template
//' @return This function returns vector of k-th order moments for normally 
//' distributed random variable with mean = \code{mean} and standard 
//' deviation = \code{sd} under \code{x_lower} and \code{x_upper} truncation 
//' points \code{x_lower} and \code{x_upper} correspondingly. 
//' If \code{return_all_moments} is \code{TRUE} then see this argument 
//' description above for output details.
//' @examples
//' ## Calculate 5-th order moment of three truncated normal random  
//' ## variables (x1, x2, x3) which mean is 5 and standard deviation is 3. 
//' ## These random variables truncation points are given 
//' ## as follows:-1<x1<1, 0<x2<2, 1<x3<3.
//' k <- 3
//' x_lower <- c(-1, 0, 1, -Inf, -Inf)
//' x_upper <- c(1, 2 , 3, 2, Inf)
//' mean <- 3
//' sd <- 5
//' 
//' # get the moments
//' truncatedNormalMoment(k, x_lower, x_upper, mean, sd)
//'
//' # get matrix of (0-5)-th moments (columns) for each variable (rows)
//' truncatedNormalMoment(k, x_lower, x_upper, 
//'                       mean, sd, 
//'                       return_all_moments = TRUE)
//'
//' # get the moments derivatives respect to mean
//' truncatedNormalMoment(k, x_lower, x_upper, 
//'                       mean, sd, 
//'                       diff_type = "mean")
//' 
//' # get the moments derivatives respect to standard deviation
//' truncatedNormalMoment(k, x_lower, x_upper, 
//'                       mean, sd, 
//'                       diff_type = "sd")
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix truncatedNormalMoment(int k = 1,
	NumericVector x_lower = NumericVector(0),
	NumericVector x_upper = NumericVector(0),
	double mean = 0, double sd = 1,
	NumericVector pdf_lower = NumericVector(0),
	NumericVector cdf_lower = NumericVector(0),
	NumericVector pdf_upper = NumericVector(0),
	NumericVector cdf_upper = NumericVector(0),
	NumericVector cdf_difference = NumericVector(0),
	bool return_all_moments = false, 
	bool is_validation = true,
	bool is_parallel = false,
	String diff_type = "NO") 
{
	// Assign default truncation values
	double max_value = 0.99*std::numeric_limits<double>::max();
  double min_value = -max_value;

	if (x_lower.size() == 0)
	{
		x_lower = NumericVector::create(min_value);
	}

	if (x_upper.size() == 0)
	{
		x_lower = NumericVector::create(max_value);
	}

	// Get number of observations
	int n = x_lower.size();

	// Validation
	//------------------------------------------------------------
	if (is_validation)
	{
		if (k < 0)
		{
			stop("parameter k should be non-negative integer");
		}
		if (sd <= 0)
		{
			stop("parameter sd should be positive integer");
		}
		if ((x_upper.size() != n) & (x_upper[0] != max_value))
		{
			stop("vectors x_lower and x_upper should have the same length");
		}
		if ((x_lower.size() != n) & (x_lower[0] != min_value))
		{
			stop("vectors x_lower and x_upper should have the same length");
		}
		if ((pdf_lower.size() != n) & (pdf_lower.size() != 0))
		{
			stop("vectors x_lower and pdf_lower should have the same length");
		}
		if ((cdf_lower.size() != n) & (cdf_lower.size() != 0))
		{
			stop("vectors x_lower and cdf_lower should have the same length");
		}
		if ((pdf_upper.size() != n) & (pdf_upper.size() != 0))
		{
			stop("vectors x_lower and pdf_upper should have the same length");
		}
		if ((cdf_upper.size() != n) & (cdf_upper.size() != 0))
		{
			stop("vectors x_lower and cdf_upper should have the same length");
		}
		if ((cdf_difference.size() != n) & (cdf_difference.size() != 0))
		{
			stop("vectors x_lower and cdf_difference should have the same length");
		}
		if((diff_type != "NO") & (diff_type != "mean") & (diff_type != "sd") &
       (diff_type != "x_lower") & (diff_type != "x_upper"))
		{
		  stop("diff_type argument should take value 'NO', 'mean', 'sd', 'x_lower' or 'x_upper'");
		}
		if(sum(x_lower >= x_upper) > 0)
		{
		  stop("x_lower values should not be greater then x_upper values");
		}
	}
	//------------------------------------------------------------
	
	// Store some information during calculations in order to
	// use if for the differentiation if need
	NumericMatrix x_pow_prod_pdf_upper = NumericMatrix(n, k + 1);
	NumericMatrix x_pow_prod_pdf_lower = NumericMatrix(n, k + 1);
	NumericMatrix x_pow_prod_pdf_difference = NumericMatrix(n, k + 1);
	NumericVector cdf_difference_squared;
	NumericVector x_upper_adj;
	NumericVector x_lower_adj;
	NumericVector x_upper_adj_2;
	NumericVector x_lower_adj_2;

	// Initialize matrix to store the moments 
	// and their derivatives if need
	NumericMatrix tr_moments(n, k + 1);
	NumericMatrix tr_moments_diff(n, k + 1);

	std::fill(tr_moments.begin(), tr_moments.end(), 1);

	// If order is 0 just return 1 for  
	// the moment or 0 for its derivatives
	if (k == 0)
	{
	  if(diff_type != "NO")
	  {
	    return(tr_moments_diff);
	  }
		return(tr_moments);
	}

	// PDF calculation (if not provided)
	if (pdf_lower.size() == 0)
	{
	  pdf_lower = dnorm_parallel(x_lower, mean, sd, is_parallel);
	}
	if (pdf_upper.size() == 0)
	{
	  pdf_upper = dnorm_parallel(x_upper, mean, sd, is_parallel);
	}

	// CDF calculation (if not provided)
	if (cdf_difference.size() == 0)
	{
		if (cdf_lower.size() == 0)
		{
		  cdf_lower = pnorm_parallel(x_lower, mean, sd, is_parallel);
		}
		if (cdf_upper.size() == 0)
		{
		  cdf_upper = pnorm_parallel(x_upper, mean, sd, is_parallel);
		}
		cdf_difference = cdf_upper - cdf_lower;
	}
	
	// Set infinity to zero in order to nullify power * pdf
	LogicalVector lower_cond = is_infinite(x_lower);
	LogicalVector upper_cond = is_infinite(x_upper);
	NumericVector x_lower_new = clone(x_lower);
	NumericVector x_upper_new = clone(x_upper);
	x_lower_new[lower_cond] = 0;
	x_upper_new[upper_cond] = 0;
	
	
	// Prepare some values
	double sd_squared = pow(sd, 2);
	
	if(diff_type != "NO")
	{
	  cdf_difference_squared = cdf_difference * cdf_difference;
	  
	  x_upper_adj = (x_upper_new - mean) / sd_squared;
	  x_lower_adj = (x_lower_new - mean) / sd_squared;
	  x_upper_adj_2 = (pow(x_upper_adj, 2) - 1 / sd_squared) * sd;
	  x_lower_adj_2 = (pow(x_lower_adj, 2) - 1 / sd_squared) * sd;
	}

	// The first moment
	NumericVector pdf_difference = pdf_upper - pdf_lower;
	tr_moments(_, 1) = mean - sd_squared * 
	                   pdf_difference / cdf_difference;
	
	// The first moment derivative
	
	  // for mean
	if(diff_type == "mean")
	{
	  tr_moments_diff(_, 1) = 1 - sd_squared * 
	                          ((pdf_upper * x_upper_adj - 
	                            pdf_lower * x_lower_adj) * cdf_difference +
	                            pdf_difference * pdf_difference) / 
	                          cdf_difference_squared;
	}

	  // for sd
	if(diff_type == "sd")
	{
	  tr_moments_diff(_, 1) = -2 * sd * pdf_difference / cdf_difference - 
	                          sd_squared * 
                      	    ((pdf_upper * x_upper_adj_2 - 
                      	      pdf_lower * x_lower_adj_2) *
                      	     cdf_difference + sd *
                      	     (pdf_upper * x_upper_adj -
                      	      pdf_lower * x_lower_adj) * 
                      	     pdf_difference) / 
                      	    cdf_difference_squared;
	}
	
	  // for x_upper
	if(diff_type == "x_upper")
	{
	  tr_moments_diff(_, 1) = sd_squared * pdf_upper *
	                          (x_upper_adj * cdf_difference +
	                           pdf_difference) / cdf_difference_squared;
	}
	
	  // for x_lower
	if(diff_type == "x_lower")
	{
	  tr_moments_diff(_, 1) = -sd_squared * pdf_lower *
	    (x_lower_adj * cdf_difference +
	    pdf_difference) / cdf_difference_squared;
	}

	// Recursively calculate other moments
	for (int i = 2; i <= k; i++)
	{
	  x_pow_prod_pdf_upper(_, i) = pow(x_upper_new, i - 1) * pdf_upper;
	  x_pow_prod_pdf_lower(_, i) = pow(x_lower_new, i - 1) * pdf_lower;
	  
	  x_pow_prod_pdf_difference(_, i) = (x_pow_prod_pdf_upper(_, i) - 
	                                     x_pow_prod_pdf_lower(_, i));
	  
		tr_moments(_, i) = (i - 1) * sd_squared * tr_moments(_, i - 2) +
			                 mean * tr_moments(_, i - 1) - sd_squared *
			                 (x_pow_prod_pdf_difference(_, i) / cdf_difference);
	}
	
	if(diff_type == "mean")
	{
	  for (int i = 2; i <= k; i++)
	  {
  	  tr_moments_diff(_, i) = (i - 1) * sd_squared * tr_moments_diff(_, i - 2) +
  	                          mean * tr_moments_diff(_, i - 1) + 
  	                          tr_moments(_, i - 1) - sd_squared *
  	                          (((x_pow_prod_pdf_upper(_, i) * x_upper_adj -
  	                             x_pow_prod_pdf_lower(_, i) * x_lower_adj) * 
  	                            cdf_difference) +
  	                            x_pow_prod_pdf_difference(_, i) *
  	                            pdf_difference) / cdf_difference_squared;
	  }
	}

	if(diff_type == "sd")
	{
	  for (int i = 2; i <= k; i++)
	  {
	    tr_moments_diff(_, i) = (i - 1) * 
	                            (2 * sd * tr_moments(_, i - 2) 
                               + sd_squared * tr_moments_diff(_, i - 2)) +
                               mean * tr_moments_diff(_, i - 1) - 
                               2 * sd * x_pow_prod_pdf_difference(_, i) / 
                               cdf_difference - sd_squared *
                               ((x_pow_prod_pdf_upper(_, i) * x_upper_adj_2 -
                                 x_pow_prod_pdf_lower(_, i) * x_lower_adj_2) *
                                cdf_difference + sd *
                                (pdf_upper * x_upper_adj -
                                 pdf_lower * x_lower_adj) * 
                                 x_pow_prod_pdf_difference(_, i)) / 
                               cdf_difference_squared;
	  }
	}
	
	if(diff_type == "x_upper")
	{
	  for (int i = 2; i <= k; i++)
	  {
	    NumericVector x_upper_adj_3 = pow(x_upper, i - 2) * 
                            	      (sd_squared * (i - 1) - 
                            	      x_upper * (x_upper - mean)) / 
                            	      sd_squared;
	    tr_moments_diff(_, i) = (i - 1) * sd_squared * tr_moments_diff(_, i - 2) +
	                            mean * tr_moments_diff(_, i - 1) - sd_squared * 
	                            (pdf_upper * x_upper_adj_3 * cdf_difference - 
	                             pdf_upper * x_pow_prod_pdf_difference(_, i)) / 
	                             cdf_difference_squared;
	                            
	  }
	}
	
	if(diff_type == "x_lower")
	{
	  for (int i = 2; i <= k; i++)
	  {
	    NumericVector x_lower_adj_3 = pow(x_lower, i - 2) * 
                            	      (sd_squared * (i - 1) - 
                            	      x_lower * (x_lower - mean)) / 
                            	      sd_squared;
	    tr_moments_diff(_, i) = (i - 1) * sd_squared * tr_moments_diff(_, i - 2) +
                      	      mean * tr_moments_diff(_, i - 1) + sd_squared * 
                      	      (pdf_lower * x_lower_adj_3 * cdf_difference - 
                      	      pdf_lower * x_pow_prod_pdf_difference(_, i)) / 
                      	      cdf_difference_squared;
	    
	  }
	}
	
	// If return_all_moments is TRUE then 
	// return matrix of all moments from 0 to k
	if (return_all_moments)
	{
	  if(diff_type != "NO")
	  {
	    return(tr_moments_diff);
	  }
		return(tr_moments);
	}

	// If return_all_moments is FALSE then return k-th moment only
	if(diff_type != "NO")
	{
	  NumericMatrix tr_moments_diff_new(n, 1);
	  tr_moments_diff_new(_, 0) = tr_moments_diff(_, k);
	  
	  return(tr_moments_diff_new);
	}
	
	NumericMatrix tr_moments_new(n, 1);
	tr_moments_new(_, 0) = tr_moments(_, k);

	return(tr_moments_new);
}
