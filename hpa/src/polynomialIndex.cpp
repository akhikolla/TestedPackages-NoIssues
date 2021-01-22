#include "polynomialIndex.h"
#include "hpaValidation.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]

//' Multivariate Polynomial Representation
//' @name polynomialIndex
//' @description Function \code{\link[hpa]{polynomialIndex}} 
//' provides matrix which allows to iterate through the elements 
//' of multivariate polynomial being aware of these elements powers. 
//' So (i, j)-th element of the matrix is power of j-th variable in i-th 
//' multivariate polynomial element.
//' 
//' Function \code{\link[hpa]{printPolynomial}} prints multivariate polynomial
//' given its degrees (\code{pol_degrees}) and coefficients 
//' (\code{pol_coefficients}) vectors.
//' @template pol_degrees_Template
//' @template pol_coefficients_Template
//' @details Multivariate polynomial of degrees   
//' \eqn{(K_{1},...,K_{m})} (\code{pol_degrees}) has the form:
//' \deqn{a_{(0,...,0)}x_{1}^{0}*...*x_{m}^{0}+ ... + 
//' a_{(K_{1},...,K_{m})}x_{1}^{K_{1}}*...*x_{m}^{K_{m}},}
//' where \eqn{a_{(i_{1},...,i_{m})}} are polynomial coefficients, while
//' polynomial elements are:
//' \deqn{a_{(i_{1},...,i_{m})}x_{1}^{i_{1}}*...*x_{m}^{i_{m}},}
//' where \eqn{(i_{1},...,i_{m})} are polynomial element's powers corresponding
//' to variables \eqn{(x_{1},...,x_{m})} respectively. Note that 
//' \eqn{i_{j}\in \{0,...,K_{j}\}}. 
//' 
//' Function \code{\link[hpa]{printPolynomial}} removes polynomial elements 
//' which coefficients are zero and variables which powers are zero. Output may 
//' contain long coefficients representation as they are not rounded.
//' @return Function \code{\link[hpa]{polynomialIndex}} 
//' returns matrix which rows are 
//' responsible for variables while columns are related to powers. 
//' So \eqn{(i, j)}-th element of this matrix corresponds to the 
//' power \eqn{i_{j}} of the \eqn{x_{j}} variable in \eqn{i}-th polynomial 
//' element. Therefore \eqn{i}-th column of this matrix contains vector of
//' powers \eqn{(i_{1},...,i_{m})} for the \eqn{i}-th polynomial element.
//' So the function transforms \eqn{m}-dimensional elements indexing
//' to one-dimensional.
//' 
//' Function \code{\link[hpa]{printPolynomial}} returns the string which 
//' contains polynomial symbolic representation.
//' @template polynomialIndex_examples_Template
//' @template printPol_examples_Template
//' @template is_validation_Template
//' @export
// [[Rcpp::export]]
NumericMatrix polynomialIndex(NumericVector pol_degrees = NumericVector(0),
                              bool is_validation = true) 
{
  // Validation stuff
  if(is_validation)
  {
    pol_Validate(pol_degrees, NumericVector(0));
  }
  
  // Convert pol_degrees to std vector of integer values
	std::vector<int> degrees = as<std::vector<int> >(pol_degrees);

	// Initiale degrees related variables
	int degrees_size = degrees.size();
	std::vector<int> degrees_products(degrees_size);

	// Calculate number of coefficients and degrees products
	int coefficients_size = 1;

	for (int i = 0; i < degrees_size; i++)
	{
		coefficients_size *= (degrees[i] + 1);      // +1 because degrees order 
		degrees_products[i] = 1;                    // starts from zero

		for (int j = i + 1; j < degrees_size; j++)
		{
			degrees_products[i] *= (degrees[j] + 1);
		}
	}

	// Assign vector index to coefficients
	std::vector<std::vector<int>> ind_pattern_full = std::vector<std::vector<int>>(degrees_size);
	std::vector<std::vector<int>> coefficients_ind(coefficients_size);
	NumericMatrix coefficients_vec(degrees_size, coefficients_size);

	for (int i = 0; i < degrees_size; i++)
	{
		// Calculate pattern for i-th variable
		std::vector<int> ind_pattern = std::vector<int>(degrees_products[i] * 
		                               (degrees[i] + 1));
		int counter = 0;

		for (int j = 0; j <= degrees[i]; j++)
		{
			for (int k = 0; k < degrees_products[i]; k++)
			{
				ind_pattern[counter] = j;
				counter++;
			}
		}

		int ind_pattern_times = coefficients_size / ind_pattern.size(); //times pattern repeats
		ind_pattern_full[i].reserve(coefficients_size); //preallocate memorry to increase insertation speed

		for (int j = 0; j < ind_pattern_times; j++)
		{
			// Repeat pattern untill the end of the pattern matrix row
			ind_pattern_full[i].insert(ind_pattern_full[i].end(), 
                                 ind_pattern.begin(), ind_pattern.end());
		}

		// Pattern defined for rows while coefficients indexes are located in columns
		// of pattern matrix. Lets copy values from patterns to coefficients indexes
		for (int j = 0; j < coefficients_size; j++)
		{
			coefficients_vec(i, j) = (double)(ind_pattern_full[i][j]);
		}
	}

	return(coefficients_vec);
}

//' @name polynomialIndex
//' @export
// [[Rcpp::export]]
Rcpp::String printPolynomial(NumericVector pol_degrees, 
                             NumericVector pol_coefficients,
                             bool is_validation = true)
{
  // Validation stuff
  if(is_validation)
  {
    pol_Validate(pol_degrees, pol_coefficients);
  }
  
	// Load R environments
	Rcpp::Environment base_env("package:base");
	Rcpp::Function paste0 = base_env["paste0"];

	// Get polynomial matrix from polynomialIndex function
	NumericVector pol_ind_mat = polynomialIndex(pol_degrees);

	// Create dimensions related variables
	int pol_coefficients_n = pol_coefficients.size();
	int pol_degrees_n = pol_degrees.size();

	// Initialize variable to store the polynomial symbolic representation
	std::string pol_string = "";

	// Iteratite throught polynomial coefficients and variables
	for (int i = 0; i < pol_coefficients_n; i++)
	{
		if ((pol_coefficients[i] != 0) | (i == 0))
		{
			if ((pol_coefficients[i] != 1) | (i == 0))
			{
				String pol_string_R = paste0((double)pol_coefficients[i]);
				pol_string += pol_string_R;
			}
			for (int j = 0; j < pol_degrees_n; j++)
			{
				int pol_pow = pol_ind_mat(j, i);
				if (pol_pow != 0)
				{
					pol_string += "x" + std::to_string(j + 1);
					if (pol_pow != 1)
					{
						pol_string += "^" + std::to_string(pol_pow);
					}
				}
			}
		}
			if (i < (pol_coefficients_n - 1))
			{
				if (pol_coefficients[i + 1] > 0)
				{
					pol_string += " + ";
				}
				if (pol_coefficients[i + 1] < 0)
				{
					pol_coefficients[i + 1] = -pol_coefficients[i + 1];
					pol_string += " - ";
				}
			}
	}

	return(pol_string);
}
