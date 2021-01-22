#include "polynomialIndex.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]

void pol_Validate(NumericVector pol_degrees = NumericVector(0),
                  NumericVector pol_coefficients = NumericVector(0))
{
  int n = pol_degrees.size();
  int m = pol_coefficients.size();
  
  // Validation for polynomial degrees
  
  bool pol_degrees_is_na = any(is_na(pol_degrees));
  bool pol_degrees_is_nan = any(is_nan(pol_degrees));
  bool pol_degrees_is_epmpy = n == 0;
  
  bool pol_degrees_is_notposint = false;
  for(int i = 0; i < n; i ++)
  {
    if((floor(pol_degrees[i]) != pol_degrees[i]) |
       (pol_degrees[i] < 0))
    {
      pol_degrees_is_notposint = true;
      break;
    }
  }
  
  if(pol_degrees_is_na | pol_degrees_is_nan | 
     pol_degrees_is_epmpy | pol_degrees_is_notposint)
  {
      stop("pol_degrees should be not empty vector of non-negative integer values.");
  }
  
  // Validation for polynomial coefficients

  if(m > 0)
  {
    int pol_degrees_prod = 1;
    for(int i = 0; i < n; i++)
    {
      pol_degrees_prod *= (pol_degrees[i] + 1.0);
    }
    
    if(pol_degrees_prod != m)
    {
      stop("pol_coefficients length do not much pol_degrees elements. Please insure that: length(pol_degrees) == prod(pol_coefficients + 1).");
    }
    
    bool pol_coefficients_is_na = any(is_na(pol_coefficients));
    bool pol_coefficients_is_nan = any(is_nan(pol_coefficients));
    
    if(pol_coefficients_is_na | pol_coefficients_is_nan)
    {
      warning("pol_coefficients contains NA and (or) NaN values.");
    }
  }
}

void ind_Validate(LogicalVector given_ind,
                  LogicalVector omit_ind)
{
  // Get sizes of given_ind and omit_ind
  int n_given = given_ind.size();
  int n_omit = omit_ind.size();
  
  // Check that at least one element is FALSE
  if(n_given != 0)
  {
    if(sum(given_ind) == n_given)
    {
      stop("At least one given_ind component should be FALSE.");
    }
  }
  
  if(n_omit != 0)
  {
    if(sum(omit_ind) == n_omit)
    {
      stop("At least one omit_ind component should be FALSE.");
    }
  }
  
  // Check that given_ind and omit_ind are match
  if((n_given != 0) & (n_omit != 0))
  {
    LogicalVector is_both_true = given_ind & omit_ind;
    for(int i = 0; i < n_given; i ++)
    {
      bool is_both_true_bool = is_both_true[i];
      if(is_both_true_bool)
      {
        stop("Ambiguity since for some 'i' both given_ind[i] and omit_ind[i] are TRUE.");
      }
      if(sum(given_ind + omit_ind) == n_given)
      {
        stop("At least one omit_ind or given_ind component should be FALSE.");
      }
    }
    
    if(n_given != n_omit)
    {
      stop("given_ind and omit_ind should be of the same size.");
    }
  }
}

void mean_Validate(NumericVector mean)
{
  if(mean.size() != 0)
  {
    bool mean_is_na = any(is_na(mean));
    bool mean_is_nan = any(is_nan(mean));
    
    if(mean_is_na | mean_is_nan)
    {
      warning("mean contains NA or NaN values.");
    }
  }
}

void sd_Validate(NumericVector sd)
{
  int n = sd.size();
  
  if(n != 0)
  {
    bool sd_is_na = any(is_na(sd));
    bool sd_is_nan = any(is_nan(sd));
    
    if(sd_is_na | sd_is_nan)
    {
      warning("sd contains NA or NaN values.");
    }
    
    for(int i = 0; i < n; i ++)
    {
      if(sd[i] <= 0)
      {
        stop("sd should not contain zero or negative values.");
      }
    }
  }
}
  
void expectation_powers_Validate(NumericVector expectation_powers)
{
  int n = expectation_powers.size();
    
  // Validation for polynomial degrees
    
  bool expectation_powers_is_na = any(is_na(expectation_powers));
  bool expectation_powers_is_nan = any(is_nan(expectation_powers));
    
  bool expectation_powers_is_notint = false;
  for(int i = 0; i < n; i ++)
  {
    if(floor(expectation_powers[i]) != expectation_powers[i])
    {
      expectation_powers_is_notint = true;
      break;
    }
  }
    
  if(expectation_powers_is_na | expectation_powers_is_nan | 
     expectation_powers_is_notint)
  {
    stop("expectation_powers should be a vector of non-negative integer values.");
  }
}
