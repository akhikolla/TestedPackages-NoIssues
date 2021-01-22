#include "normalMoments.h"
#include "polynomialIndex.h"
#include "hpaMain.h"
#include "ParallelFunctions.h"
#include "hpaValidation.h"
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace RcppParallel;

// Hermite polynomial density,
// cumulative distribution function and moments approximations.
List hpaMain(
	NumericMatrix x_lower = NumericMatrix(1,1),
	NumericMatrix x_upper = NumericMatrix(1,1),
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	String type = "pdf",
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	NumericVector expectation_powers = NumericVector(0),
	String grad_type = "NO",
	bool is_parallel = false,
	bool is_cdf = false,
  bool log = false,
  bool is_validation = true)
{
	// Get number of observations
	int n = x_upper.nrow();
  
  // Validation Stuff
  if (is_validation)
  {
    int target_dim = x_upper.ncol();
    
    // Validate polynomial coefficients
    pol_Validate(pol_degrees, pol_coefficients);
    
    if((type != "expectation") & (type != "expectation truncated"))
    {
      if (pol_degrees.size() != target_dim)
      {
        if (type == "pdf")
        {
          stop("pol_degrees length should be the same as the number of x columns.");
        } else {
          stop("pol_degrees length should be the same as the number of x_upper columns.");
        }
      }
    } else{
      target_dim = pol_degrees.size();
    }
    
    // Validate conditional and omitted values
    ind_Validate(given_ind, omit_ind);
    
    int n_given_ind = given_ind.size();
    if((n_given_ind != 0) & (n_given_ind != target_dim))
    {
      stop("given_ind length should be the same as the length of pol_degrees.");
    }
    
    int n_omit_ind = omit_ind.size();
    if((n_omit_ind != 0) & (n_omit_ind != target_dim))
    {
      stop("omit_ind length should be the same as the length of pol_degrees.");
    }
    
    // Validate mean and sd vectors
    int n_mean = mean.size();
    if((n_mean != 0) & (n_mean != target_dim))
    {
      stop("mean length should be the same as the length of pol_degrees.");
    }
    mean_Validate(mean);
    
    int n_sd = sd.size();
    if((n_sd != 0) & (n_sd != target_dim))
    {
      stop("sd length should be the same as the length of pol_degrees.");
    }
    sd_Validate(sd);
    
    // Validate expectations powers vectors
    int n_expectation_powers = expectation_powers.size();
    if((n_expectation_powers != 0) & 
       (n_expectation_powers != target_dim))
    {
      stop("expectation_powers length should be the same as the length of pol_degrees.");
    }
    expectation_powers_Validate(expectation_powers);
    
    // Validate x_lower and x_upper
    if((type != "expectation") & (type != "expectation truncated"))
    {
      if (((x_lower.nrow() > 1) | (x_lower.ncol() > 1)) & !is_cdf)
      {
        if((x_lower.ncol() != x_upper.ncol()))
        {
          stop("x_lower should have the same number of columns as x_upper.");
        }
        if((x_lower.nrow() != x_upper.nrow()))
        {
          stop("x_lower should have the same number of rows as x_upper.");
        }
      }
    } else {
      if ((x_upper.nrow() > 1) | (x_upper.ncol() > 1))
      {
        if(x_upper.ncol() != target_dim)
        {
          stop("pol_degrees length should be the same as the number of x columns.");
        }
      }
    }
  }

	// Initialize polynomial structure related values
	int pol_degrees_n = pol_degrees.size();
	int pol_coefficients_n = pol_coefficients.size();

	// Fill x_lower with (-INF) if need
	if (!((type == "interval") | (type == "expectation truncated")) | is_cdf)
	{
		x_lower = NumericMatrix(n, pol_degrees_n);
		std::fill(x_lower.begin(), x_lower.end(), R_NegInf);
	}

	// Fill given_ind and omit_ind with defaults if need
	if (given_ind.size() == 0)
	{
		given_ind = LogicalVector(pol_degrees_n);
	}
	
	if (omit_ind.size() == 0)
	{
		omit_ind = LogicalVector(pol_degrees_n);
	}
	
	// Fill mean and sd with defaults if need
	if (mean.size() == 0)
	{
		mean = NumericVector(pol_degrees_n);
		std::fill(mean.begin(), mean.end(), 0);
	}

	if (sd.size() == 0)
	{
		sd = NumericVector(pol_degrees_n);
		std::fill(sd.begin(), sd.end(), 1);
	}

	// Control for the expected powered product powers values
	if ((expectation_powers.size() == 0) | ((type != "expectation") & 
                                          (type != "expectation truncated")))
	{
		expectation_powers = NumericVector(pol_degrees_n);
		std::fill(expectation_powers.begin(), expectation_powers.end(), 0);
	}

	// Initialize indexes for observable unconditioned components
	LogicalVector d_cond = ((!given_ind) & (!omit_ind));

	// Define vectors related to normal distribution
	
	  // pdf associated values
	NumericMatrix pdf_upper = NumericMatrix(n, pol_degrees_n);
	NumericMatrix pdf_lower = NumericMatrix(n, pol_degrees_n);
	NumericVector pdf_product(n);
	  
	  // cdf associated values
	NumericMatrix cdf_upper = NumericMatrix(n, pol_degrees_n);
	NumericMatrix cdf_lower = NumericMatrix(n, pol_degrees_n);
	
	NumericMatrix cdf_difference = NumericMatrix(n, pol_degrees_n);
	NumericMatrix pdf_difference;
	
	NumericVector cdf_difference_product = NumericVector(n);
	std::fill(cdf_difference_product.begin(), 
            cdf_difference_product.end(), 
            1);
	
	// control for zero moments during numeric differentiation if need
	if((grad_type == "mean") | (grad_type == "sd") | (grad_type == "all"))
	{
	  mean[mean == 0] = std::numeric_limits<double>::epsilon();
	}
	
	// control for zero x_lower during numeric differentiation if need
	if((grad_type == "x_lower") | (grad_type == "all"))
	{
	  for(int i = 0; i < pol_degrees_n; i++)
	  {
	    NumericVector x_i_tmp = x_lower(_, i);
	    x_i_tmp[x_i_tmp == 0] = std::numeric_limits<double>::epsilon();
	    x_lower(_, i) = x_i_tmp;
	  }
	}
	
	// control for zero x_upper during numeric differentiation if need
	// and substitute x for x_upper
	if(grad_type == "x")
	{
	  grad_type = "x_upper";
	}
	
	if((grad_type == "x_upper") | (grad_type == "all"))
	{
	  for(int i = 0; i < pol_degrees_n; i++)
	  {
	    NumericVector x_i_tmp = x_upper(_, i);
	    x_i_tmp[x_i_tmp == 0] = std::numeric_limits<double>::epsilon();
	    x_upper(_, i) = x_i_tmp;
	  }
	}

	if (type != "expectation")
	{
		// Initialize densities

  		// Upper densities
  			pdf_upper = NumericMatrix(n, pol_degrees_n);
  			for (int i = 0; i < pol_degrees_n; i++)
  			{
  				if (!omit_ind[i])
  				{
  				  pdf_upper(_, i) = dnorm_parallel(x_upper(_, i), 
                                             mean[i], sd[i], 
                                             is_parallel);
  				}
  			}

  		// Lower densities
  		if ((type == "interval") | (type == "expectation truncated"))
  		{
  			for (int i = 0; i < pol_degrees_n; i++)
  			{
  				if (!omit_ind[i])
  				{
  				  pdf_lower(_, i) = dnorm_parallel(x_lower(_, i), 
                                             mean[i], sd[i], 
                                             is_parallel);
  				}
  			}
  		}

  		// Product of densities
  		if (type == "pdf")
  		{
  			std::fill(pdf_product.begin(), 
                  pdf_product.end(), 
                  1);
  			for (int i = 0; i < pol_degrees_n; i++)
  			{
  				if (d_cond[i])
  				{
  					pdf_product = pdf_product * pdf_upper(_, i);
  				}
  			}
  		}

		// Initialize cumulative distribution functions (cdfs)

  		// Upper cdf
  		if (type != "pdf")
  		{
  			if (cdf_upper(0, 0) == 0)
  			{
  				for (int i = 0; i < pol_degrees_n; i++)
  				{
  					if (d_cond[i])
  					{
  					  cdf_upper(_, i) = pnorm_parallel(x_upper(_, i), 
                                               mean[i], sd[i], 
                                               is_parallel);
  					}
  				}
  			}
  		}

  		// Lower cdf
  		if (((type == "interval") | (type == "expectation truncated")))
  		{
  			cdf_lower = NumericMatrix(n, pol_degrees_n);
  			for (int i = 0; i < pol_degrees_n; i++)
  			{
  				if (d_cond[i])
  				{
  					cdf_lower(_, i) = pnorm_parallel(x_lower(_, i), 
                                             mean[i], sd[i], 
                                             is_parallel);
  				}
  			}
  		}

  		// Calculate cdf_difference and pdf_difference if need
  		// for gradient calculations
  		if((grad_type == "mean") | (grad_type == "x_lower") |
         (grad_type == "x_upper") | (grad_type == "all"))
  		{
  		  pdf_difference = NumericMatrix(n, pol_degrees_n);
  		}
  		for (int i = 0; i < pol_degrees_n; i++)
  		{
  			if (d_cond[i])
  			{
  				cdf_difference(_, i) = cdf_upper(_, i) - cdf_lower(_, i);
  			  if((grad_type == "mean") | (grad_type == "x_lower") |
             (grad_type == "x_upper") | (grad_type == "all"))
  			  {
  			    pdf_difference(_, i) = pdf_upper(_, i) - pdf_lower(_, i);
  			  }
  			}
  		}

  		// Estimate cdf_difference product
  		if (type != "pdf")
  		{
  			for (int i = 0; i < pol_degrees_n; i++)
  			{
  				if (d_cond[i])
  				{
  					cdf_difference_product = cdf_difference_product * 
  					                         cdf_difference(_, i);
  				}
  			}
  		}
	}

	// Define vector indexing system for polynomial
	NumericMatrix polynomial_index = polynomialIndex(pol_degrees, false);

	// Calculate moments
	List moments(pol_degrees_n);            // for function value
	List moments_diff_mean(pol_degrees_n);  // for derivative w.r.t. mean
	List moments_diff_sd(pol_degrees_n);    // for derivative w.r.t. sd
	
	int max_degree;

	for (int i = 0; i < pol_degrees_n; i++)
	{
		if (!given_ind[i])
		{
			max_degree = 2 * pol_degrees[i] + expectation_powers[i];
			moments[i] = normalMoment(max_degree,
									              mean[i], sd[i], 
									              true, false, 
									              false, "NO");
			if((grad_type == "mean") | (grad_type == "all"))
			{
			  moments_diff_mean[i] = normalMoment(max_degree,
                                            mean[i], sd[i], 
                                            true, false, 
                                            false, "mean");
			}
			if((grad_type == "sd") | (grad_type == "all"))
			{
			  moments_diff_sd[i] = normalMoment(max_degree,
                                          mean[i], sd[i], 
                                          true, false, 
                                          false, "sd");
			}
		}
	}

	// Calculate truncated moments
	List tr_moments(pol_degrees_n);                // for function value
	List tr_moments_diff_mean(pol_degrees_n);      // for derivative w.r.t. mean
	List tr_moments_diff_sd(pol_degrees_n);        // for derivative w.r.t. sd
	List tr_moments_diff_x_upper(pol_degrees_n);   // for derivative w.r.t. x_upper
	List tr_moments_diff_x_lower(pol_degrees_n);   // for derivative w.r.t. x_lower

	if ((type != "pdf") & (type != "expectation"))
	{
		for (int i = 0; i < pol_degrees_n; i++)
		{
			if (d_cond[i])
			{
				max_degree = 2 * pol_degrees[i] + expectation_powers[i];
				tr_moments[i] = truncatedNormalMoment(max_degree,
					x_lower(_, i), x_upper(_, i),
					mean[i], sd[i],
					pdf_lower(_, i), cdf_lower(_, i),
					pdf_upper(_, i), cdf_upper(_, i),
					cdf_difference(_, i), true, false, is_parallel, "NO");
				if((grad_type == "mean") | (grad_type == "all"))
				{
				  tr_moments_diff_mean[i] = truncatedNormalMoment(max_degree,
            x_lower(_, i), x_upper(_, i),
            mean[i], sd[i],
            pdf_lower(_, i), cdf_lower(_, i),
            pdf_upper(_, i), cdf_upper(_, i),
            cdf_difference(_, i), true, false, is_parallel, "mean");
				}
				if((grad_type == "sd") | (grad_type == "all"))
				{
				  tr_moments_diff_sd[i] = truncatedNormalMoment(max_degree,
            x_lower(_, i), x_upper(_, i),
            mean[i], sd[i],
            pdf_lower(_, i), cdf_lower(_, i),
            pdf_upper(_, i), cdf_upper(_, i),
            cdf_difference(_, i), true, false, is_parallel, "sd");
				}
				if((grad_type == "x_upper") | (grad_type == "all"))
				{
				  tr_moments_diff_x_upper[i] = truncatedNormalMoment(max_degree,
            x_lower(_, i), x_upper(_, i),
            mean[i], sd[i],
            pdf_lower(_, i), cdf_lower(_, i),
            pdf_upper(_, i), cdf_upper(_, i),
            cdf_difference(_, i), true, false, is_parallel, "x_upper");
				}
				if(((grad_type == "x_lower") | (grad_type == "all")) & 
           (type == "interval"))
				{
				  tr_moments_diff_x_lower[i] = truncatedNormalMoment(max_degree,
            x_lower(_, i), x_upper(_, i),
            mean[i], sd[i],
            pdf_lower(_, i), cdf_lower(_, i),
            pdf_upper(_, i), cdf_upper(_, i),
            cdf_difference(_, i), true, false, is_parallel, "x_lower");
				}
			}
		}
	}

	// Calculate truncated moments derivatives if need

	// Calculate x powers (x ^ polynomial_degree)
	LogicalVector x_cond(pol_degrees_n);

	if (type == "pdf")
	{
		x_cond = !omit_ind;
	} else {
		x_cond = given_ind;
	}

	List x_pow(pol_degrees_n);

	int k = 0;

	for (int i = 0; i < pol_degrees_n; i++)
	{
		if (x_cond[i])
		{
			k = 2 * pol_degrees[i] + 1;
			x_pow[i] = NumericMatrix(n, k);
			NumericMatrix x_pow_i = x_pow[i]; // it is reference
			NumericVector x_upper_i = x_upper(_, i);
			for (int j = 0; j < k; j++)
			{
			  if(is_parallel)
			  {
			    x_pow_i(_, j) = ParallelVectorPow(x_upper_i, j);
			  } else {
			    x_pow_i(_, j) = pow(x_upper_i, j);
			  }
			}
		}
	}

	// Calculate main expression

	// Initialize values to store temporal results
	NumericVector value_pgn(n);                          // nominator
	std::fill(value_pgn.begin(), value_pgn.end(), 0);

	NumericVector psi(n);                                // denominator
	std::fill(psi.begin(), psi.end(), 0);

	NumericVector value_sum_element(n);
	NumericVector psi_sum_element(n);
	int polynomial_sum = 0;
	
	// Initialize values to store gradient 
	// information if need
	  // for polynomial coefficients
	NumericMatrix pc_grad;
	NumericMatrix pc_grad_value;
	NumericMatrix pc_grad_psi;
	  // for mean
	NumericMatrix mean_grad;
	NumericMatrix mean_grad_value;
	NumericMatrix mean_grad_psi;
	  // for sd
	NumericMatrix sd_grad;
	NumericMatrix sd_grad_value;
	NumericMatrix sd_grad_psi;
	  // for x_upper
	NumericMatrix x_upper_grad;
	NumericMatrix x_upper_grad_value;
	NumericMatrix x_upper_grad_psi;
	  // for x_lower
	NumericMatrix x_lower_grad;
	NumericMatrix x_lower_grad_value;
	NumericMatrix x_lower_grad_psi;
	
	  // preallocate memory to store gradient specific information if need
	if((grad_type == "pol_coefficients") | (grad_type == "all"))
	{
	  pc_grad = NumericMatrix(n, pol_coefficients_n);
	  pc_grad_value = NumericMatrix(n, pol_coefficients_n);
	  pc_grad_psi = NumericMatrix(n, pol_coefficients_n);
	}
	if((grad_type == "mean") | (grad_type == "all"))
	{
	  mean_grad = NumericMatrix(n, pol_degrees_n);
	  mean_grad_value = NumericMatrix(n, pol_degrees_n);
	  mean_grad_psi = NumericMatrix(n, pol_degrees_n);
	}
	if((grad_type == "sd") | (grad_type == "all"))
	{
	  sd_grad = NumericMatrix(n, pol_degrees_n);
	  sd_grad_value = NumericMatrix(n, pol_degrees_n);
	  sd_grad_psi = NumericMatrix(n, pol_degrees_n);
	}
	if((grad_type == "x_upper") | (grad_type == "all"))
	{
	  x_upper_grad = NumericMatrix(n, pol_degrees_n);
	  x_upper_grad_value = NumericMatrix(n, pol_degrees_n);
	  x_upper_grad_psi = NumericMatrix(n, pol_degrees_n);
	}
	if(((grad_type == "x_lower") | (grad_type == "all")) & (type == "interval"))
	{
	  x_lower_grad = NumericMatrix(n, pol_degrees_n);
	  x_lower_grad_value = NumericMatrix(n, pol_degrees_n);
	  x_lower_grad_psi = NumericMatrix(n, pol_degrees_n);
	}

	// Perform main calculations
	for (int i = 0; i < pol_coefficients_n; i++)
	{
		for (int j = i; j < pol_coefficients_n; j++)
		{
		  double pol_coefficients_prod = pol_coefficients[i] * 
		                                 pol_coefficients[j];
		  
			// Initialize temporal value
			std::fill(value_sum_element.begin(), value_sum_element.end(), 1);
			std::fill(psi_sum_element.begin(), psi_sum_element.end(), 1);

			// Main calculations for each element of sum
			for (int r = 0; r < pol_degrees_n; r++)
			{
				polynomial_sum = polynomial_index(r, i) + 
				                 polynomial_index(r, j);
				if (!omit_ind[r])
				{
					if ((type == "pdf") | (given_ind[r]))
					{
						NumericMatrix x_pow_r = x_pow[r];
						value_sum_element = value_sum_element * 
						                    x_pow_r(_, polynomial_sum);
					} else {
						if (type != "expectation")
						{
							NumericMatrix tr_moments_r = tr_moments[r];
							value_sum_element = value_sum_element *
								tr_moments_r(_, polynomial_sum + expectation_powers[r]);
						} else {
							NumericVector moments_r_e = moments[r];
							value_sum_element = value_sum_element *
								moments_r_e[polynomial_sum + expectation_powers[r]];
						}
					}
				} else {
					NumericVector moments_r = moments[r];
					value_sum_element = value_sum_element * moments_r[polynomial_sum];
				}
				// psi
				if (given_ind[r])
				{
					NumericMatrix x_pow_r = x_pow[r];
					psi_sum_element = psi_sum_element * x_pow_r(_, polynomial_sum);
				} else {
					if (type != "expectation truncated")
					{
						NumericVector moments_r = moments[r];
						psi_sum_element = psi_sum_element * 
						                  moments_r[polynomial_sum];
					} else {
						NumericMatrix tr_moments_r = tr_moments[r];
						psi_sum_element = psi_sum_element * 
						                  tr_moments_r(_, polynomial_sum);
					}
				}
			}
			// Each iteration perform results storage
			int mult_for_unequal_i_j = (1 + (i != j));
			
			NumericVector value_sum_element_adj = mult_for_unequal_i_j * 
			                                      value_sum_element * 
			                                      pol_coefficients_prod;
			value_pgn = value_pgn + value_sum_element_adj;
			
			NumericVector psi_sum_element_adj = mult_for_unequal_i_j * 
			                                    psi_sum_element * 
			                                    pol_coefficients_prod;
			psi = psi + psi_sum_element_adj;

			// gradient specific storage respect to
			
			  // mean
			if((grad_type == "mean") | (grad_type == "all"))
			{
			  for (int r = 0; r < pol_degrees_n; r++)
			  {
			    if(!given_ind[r])
			    {
  			    polynomial_sum = polynomial_index(r, i) + 
  			                     polynomial_index(r, j);
  			    NumericVector moments_r = moments[r];
  			    NumericVector moments_r_diff_mean = moments_diff_mean[r];
  			    double moments_ratio = moments_r_diff_mean[polynomial_sum] /
  			                           moments_r[polynomial_sum];
    			  mean_grad_psi(_, r) = mean_grad_psi(_, r) + 
    			                        psi_sum_element_adj * moments_ratio;
  			    if(omit_ind[r])
  			    {
  			      mean_grad_value(_, r) = mean_grad_value(_, r) + 
  			                              value_sum_element_adj * moments_ratio;
  			    }
  			    if((type == "interval") & d_cond[r])
  			    {
  			      NumericMatrix tr_moments_r = tr_moments[r];
  			      NumericMatrix tr_moments_diff_mean_r = tr_moments_diff_mean[r];
  			      NumericVector tr_moments_ratio = tr_moments_diff_mean_r(_, polynomial_sum) /
  			                                       tr_moments_r(_, polynomial_sum);
  			      mean_grad_value(_, r) = mean_grad_value(_, r) +
  			                              value_sum_element_adj * tr_moments_ratio;
  			    }
			    }
			  }
			}

			  // sd
			if((grad_type == "sd") | (grad_type == "all"))
			{
			  for (int r = 0; r < pol_degrees_n; r++)
			  {
			    if(!given_ind[r])
			    {
			      polynomial_sum = polynomial_index(r, i) + 
			                       polynomial_index(r, j);
			      NumericVector moments_r = moments[r];
			      NumericVector moments_r_diff_sd = moments_diff_sd[r];
			      double moments_ratio = moments_r_diff_sd[polynomial_sum] /
			                             moments_r[polynomial_sum];
			      sd_grad_psi(_, r) = sd_grad_psi(_, r) + 
			                          psi_sum_element_adj * moments_ratio;
			      if(omit_ind[r])
			      {
			        sd_grad_value(_, r) = sd_grad_value(_, r) + 
			                              value_sum_element_adj * moments_ratio;
			      }
			      if((type == "interval") & d_cond[r])
			      {
			        NumericMatrix tr_moments_r = tr_moments[r];
			        NumericMatrix tr_moments_diff_sd_r = tr_moments_diff_sd[r];
			        NumericVector tr_moments_ratio = tr_moments_diff_sd_r(_, polynomial_sum) /
			                                         tr_moments_r(_, polynomial_sum);
			        sd_grad_value(_, r) = sd_grad_value(_, r) +
			                              value_sum_element_adj * tr_moments_ratio;
			      }
			    }
			  }
			}
			
			  // x_upper
			if((grad_type == "x_upper") | (grad_type == "all"))
			{
			  for (int r = 0; r < pol_degrees_n; r++)
			  {
			    polynomial_sum = polynomial_index(r, i) + 
              			       polynomial_index(r, j);
			    
			    if(given_ind[r])
			    {
			      x_upper_grad_psi(_, r) = x_upper_grad_psi(_, r) + 
			                               psi_sum_element_adj * polynomial_sum /
			                               x_upper(_, r);
			      x_upper_grad_value(_, r) = x_upper_grad_value(_, r) + 
                        			         value_sum_element_adj * polynomial_sum /
                        			         x_upper(_, r);
			    }
			    if(d_cond[r])
			    {
			      if(type == "interval")
			      {
			        NumericMatrix tr_moments_r = tr_moments[r];
			        NumericMatrix tr_moments_diff_x_upper_r = tr_moments_diff_x_upper[r];
			        NumericVector tr_moments_ratio = tr_moments_diff_x_upper_r(_, polynomial_sum) /
                              			           tr_moments_r(_, polynomial_sum);
			        x_upper_grad_value(_, r) = x_upper_grad_value(_, r) +
                        			           value_sum_element_adj * tr_moments_ratio;
			      } else {
			        x_upper_grad_value(_, r) = x_upper_grad_value(_, r) + 
                          			         value_sum_element_adj * polynomial_sum /
                          			         x_upper(_, r);
			      }
			    }
			  }
			}
			
			// x_lower
			if(((grad_type == "x_lower") | (grad_type == "all")) & 
         (type == ("interval")))
			{
			  for (int r = 0; r < pol_degrees_n; r++)
			  {
			    polynomial_sum = polynomial_index(r, i) + 
              			       polynomial_index(r, j);
			    
			    if(d_cond[r] & (type == "interval"))
			    {
			      NumericMatrix tr_moments_r = tr_moments[r];
			      NumericMatrix tr_moments_diff_x_lower_r = tr_moments_diff_x_lower[r];
			      NumericVector tr_moments_ratio = tr_moments_diff_x_lower_r(_, polynomial_sum) /
                                			       tr_moments_r(_, polynomial_sum);
			      x_lower_grad_value(_, r) = x_lower_grad_value(_, r) +
                      			           value_sum_element_adj * tr_moments_ratio;
			    }
			  }
			}
			
			// polynomial coefficients
			if((grad_type == "pol_coefficients") | (grad_type == "all"))
			{
			  pc_grad_value(_, i) = pc_grad_value(_, i) + mult_for_unequal_i_j * 
			                        value_sum_element * pol_coefficients[j];
			  pc_grad_value(_, j) = pc_grad_value(_, j) + mult_for_unequal_i_j * 
			                        value_sum_element * pol_coefficients[i];
			  pc_grad_psi(_, i) = pc_grad_psi(_, i) + mult_for_unequal_i_j * 
			                      psi_sum_element * pol_coefficients[j];
			  pc_grad_psi(_, j) = pc_grad_psi(_, j) + mult_for_unequal_i_j * 
			                      psi_sum_element * pol_coefficients[i];
			}
		}
	}
	
	// Return the probabilities depending on the type of calculations
	
	NumericVector return_value;
	
	if(!log | (grad_type == "NO"))
	{
  	if ((type == "expectation") | (type == "expectation truncated"))
  	{
  	  return_value = value_pgn / psi;
  	}
  	
  	if (type != "pdf")
  	{
  	  return_value = value_pgn * cdf_difference_product / psi;
  	} else {
  	  return_value = value_pgn * pdf_product / psi;
  	}
  	// Take the logarithm if need
  	if(log)
  	{
  	  return_value = Rcpp::log(return_value);
  	}
	}

	// Return gradient specific values if need

	  // for polynomial coefficients
	if ((grad_type == "pol_coefficients") | (grad_type == "all"))
	{
	  for (int i = 0; i < pol_coefficients_n; i++)
	  {
	    pc_grad(_, i) = pc_grad_value(_, i) / value_pgn -
              	      pc_grad_psi(_, i) / psi;
	    
	    if(!log)
	    {
	      pc_grad(_, i) = pc_grad(_, i) * return_value;
	    }
	  }
	  
	  if(grad_type != "all")
	  {
    	return(List::create(Named("pc_grad") = pc_grad));
	  }
	}
	
	  // for mean
	if ((grad_type == "mean") | (grad_type == "all"))
	{
	  for (int r = 0; r < pol_degrees_n; r++)
	  {
	    if(omit_ind[r] | (d_cond[r] & (type == "interval")))
	    {
	      mean_grad(_, r) = mean_grad(_, r) + 
	                        mean_grad_value(_, r) / value_pgn;
	    } 
	    
	    if(!omit_ind[r] & !given_ind[r])
	    {
	      if(type == "pdf")
	      {
  	      mean_grad(_, r) = mean_grad(_, r) +
  	                        (x_upper(_, r) - mean[r]) / 
  	                        (sd[r] * sd[r]);
	      } else {
	        mean_grad(_, r) = mean_grad(_, r) -
	                          pdf_difference(_, r) / cdf_difference(_, r);
	      }
	    }
	    
	    if(!given_ind[r])
	    {
	      mean_grad(_, r) = mean_grad(_, r) - 
	                        mean_grad_psi(_, r) / psi;
	    }
	    
	    if(!log)
	    {
  	    mean_grad(_, r) = mean_grad(_, r) * return_value;
	    }
	  }
	  
	  if(grad_type != "all")
	  {
  	  return(List::create(Named("mean_grad") = mean_grad));
	  }
	}
	
	  // for sd
	if ((grad_type == "sd") | (grad_type == "all"))
	{
	  for (int r = 0; r < pol_degrees_n; r++)
	  {
	    if(omit_ind[r] | (d_cond[r] & (type == "interval")))
	    {
	      sd_grad(_, r) = sd_grad(_, r) + 
	                      sd_grad_value(_, r) / value_pgn;
	    } 
	    
	    if(!omit_ind[r] & !given_ind[r])
	    {
	      if(type == "pdf")
	      {
	        sd_grad(_, r) = sd_grad(_, r) + 
	                        (pow((x_upper(_, r) - mean[r]) / (sd[r] * sd[r]), 2) - 
	                        1 / (sd[r] * sd[r])) * sd[r];
	      } else {
	        // Deal with infinite values
	        NumericVector x_upper_r = x_upper(_, r) * 1;
	        NumericVector x_lower_r = x_lower(_, r) * 1;
	        x_upper_r[is_infinite(x_upper_r)] = 0;
	        x_lower_r[is_infinite(x_lower_r)] = 0;
	        sd_grad(_, r) = sd_grad(_, r) - 
                	        (pdf_upper(_, r) * 
                	        (x_upper_r - mean[r]) / sd[r] -
                	         pdf_lower(_, r) * 
                	        (x_lower_r - mean[r]) / sd[r]) /
                	        cdf_difference(_, r);
	      }
	    }
	    
	    if(!given_ind[r])
	    {
	      sd_grad(_, r) = sd_grad(_, r) - 
                	      sd_grad_psi(_, r) / psi;
	    }
	    
	    if(!log)
	    {
  	    sd_grad(_, r) = sd_grad(_, r) * return_value;
	    }
	  }
	  
	  if(grad_type != "all")
	  {
	    return(List::create(Named("sd_grad") = sd_grad));
	  }
	}
	
	  // for x_upper
	if ((grad_type == "x_upper") | (grad_type == "all"))
	{
	  for (int r = 0; r < pol_degrees_n; r++)
	  {
	    if(!omit_ind[r])
	    {
	      x_upper_grad(_, r) = x_upper_grad(_, r) + 
                    	       x_upper_grad_value(_, r) / value_pgn;
	    }
	    
	    if(d_cond[r])
	    {
	      if(type == "pdf")
	      {
	        x_upper_grad(_, r) = x_upper_grad(_, r) -
                  	          (x_upper(_, r) - mean[r]) / 
                  	          (sd[r] * sd[r]);
	      } else {
	        x_upper_grad(_, r) = x_upper_grad(_, r) +
                	             pdf_upper(_, r) / cdf_difference(_, r);
	      }
	    }
	    
	    if(given_ind[r])
	    {
	      x_upper_grad(_, r) = x_upper_grad(_, r) - 
                  	         x_upper_grad_psi(_, r) / psi;
	    }
	    
	    if(!log)
	    {
  	    x_upper_grad(_, r) = x_upper_grad(_, r) * return_value;
	    }
	  }
	  
	  if(grad_type != "all")
	  {
	    return(List::create(Named("x_upper_grad") = x_upper_grad));
	  }
	}
	
  	// for x_lower
	if (((grad_type == "x_lower") | (grad_type == "all")) & 
      (type == ("interval")))
	{
	  for (int r = 0; r < pol_degrees_n; r++)
	  {
	    if(d_cond[r] & (type == "interval"))
	    {
	      x_lower_grad(_, r) = x_lower_grad(_, r) + 
	                           x_lower_grad_value(_, r) / value_pgn;
	      x_lower_grad(_, r) = x_lower_grad(_, r) -
	                           pdf_lower(_, r) / cdf_difference(_, r);
	      if(!log)
	      {
  	      x_lower_grad(_, r) = x_lower_grad(_, r) * 
                    	         cdf_difference_product * value_pgn / psi;
	      }
	    }
	  }
	  
	  if(grad_type != "all")
	  {
	    return(List::create(Named("x_lower_grad") = x_lower_grad));
	  }
	}
	
	if(grad_type == "all")
	{
	  List return_List = List::create(Named("pc_grad") = pc_grad,
            	                      Named("mean_grad") = mean_grad,
                                    Named("sd_grad") = sd_grad,
                                    Named("x_lower_grad") = x_lower_grad,
                                    Named("x_upper_grad") = x_upper_grad);
	  
	  return(return_List);
	}
	
	return(List::create(Named("values") = return_value));
}

//' Probabilities and Moments Hermite Polynomial Approximation
//' @name hpaDist
//' @template dhpa_formula_Template
//' @template x_pdf_Template
//' @template x_lower_Template
//' @template x_upper_Template
//' @template pol_coefficients_Template
//' @template pol_degrees_Template
//' @template given_ind_Template
//' @template omit_ind_Template
//' @template mean_Template
//' @template sd_Template
//' @template is_parallel_Template
//' @template log_Template
//' @template expectation_powers_Template
//' @template tr_left_Template
//' @template tr_right_Template
//' @template type_diff_Template
//' @template is_validation_Template
//' @template GN_details_Template
//' @template dhpa_examples_Template
//' @export
// [[Rcpp::export]]
NumericVector dhpa(
	NumericMatrix x,
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	bool is_parallel = false,
	bool log = false,
	bool is_validation = true)
{

  List return_List = hpaMain(
    NumericMatrix(1, 1),            // x_lower
    x,                              // x_upper
    pol_coefficients, pol_degrees,
    "pdf",                          // type
    given_ind, omit_ind,
    mean, sd,
    NumericVector(0), "NO", 
    is_parallel,
    false, log,
    is_validation);
		                   
		NumericVector return_value = return_List["values"];
		                   
		return(return_value);
}

//' @name hpaDist
//' @template phpa_examples_Template
//' @export
// [[Rcpp::export]]
NumericVector phpa(
	NumericMatrix x,
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	bool is_parallel = false,
	bool log = false,
	bool is_validation = true) 
{
  List return_List = hpaMain(
    NumericMatrix(1, 1),                     // x_lower
    x,                                       // x_upper
    pol_coefficients, pol_degrees,
    "cdf",                                   // type
    given_ind, omit_ind,
    mean, sd,
    NumericVector(0), "NO", 
    is_parallel,
    true, log,
    is_validation);
  
  NumericVector return_value = return_List["values"];
  
  return(return_value);
}

//' @name hpaDist
//' @template ihpa_examples_Template
//' @export
// [[Rcpp::export]]
NumericVector ihpa(
	NumericMatrix x_lower = NumericMatrix(1, 1),
	NumericMatrix x_upper = NumericMatrix(1, 1),
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	bool is_parallel = false,
	bool log = false,
	bool is_validation = true) 
{
  List return_List = hpaMain(
    x_lower,                                 // x_lower
    x_upper,                                 // x_upper
    pol_coefficients, pol_degrees,
    "interval",                              // type
    given_ind, omit_ind,
    mean, sd,
    NumericVector(0),"NO", 
    is_parallel,
    false, log,
    is_validation);
  
  NumericVector return_value = return_List["values"];
  
  return(return_value);
}

//' @name hpaDist
//' @template ehpa_examples_Template
//' @export
// [[Rcpp::export]]
NumericVector ehpa(NumericMatrix x = NumericMatrix(1, 1), //for given
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	NumericVector expectation_powers = NumericVector(0),
	bool is_parallel = false,
  bool is_validation = true) 
{
  List return_List = hpaMain(
    NumericMatrix(1, 1),                     // x_lower
    x,                                       // x_upper
    pol_coefficients, pol_degrees,
    "expectation",                           // type
    given_ind, omit_ind,
    mean, sd,
    expectation_powers,
    "NO", 
    is_parallel,
    false, false, is_validation);
  
  NumericVector return_value = return_List["values"];
  
  return(return_value);
}

//' @name hpaDist
//' @template etrhpa_examples_Template
//' @export
// [[Rcpp::export]]
NumericVector etrhpa(
	NumericMatrix tr_left = NumericMatrix(1, 1),
	NumericMatrix tr_right = NumericMatrix(1, 1),
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	NumericVector expectation_powers = NumericVector(0),
	bool is_parallel = false,
  bool is_validation = true) 
{
  List return_List = hpaMain(
    tr_left,                                 // x_lower
    tr_right,                                // x_upper
    pol_coefficients, pol_degrees,
    "expectation truncated",                 // type
    LogicalVector(0), LogicalVector(0),      // given_ind, omit_ind
    mean, sd,
    expectation_powers,
    "NO", 
    is_parallel,
    false, false, is_validation);
  
  NumericVector return_value = return_List["values"];
  
  return(return_value);
}

//' @name hpaDist
//' @template dtrhpa_examples_Template
//' @export
// [[Rcpp::export]]
NumericVector dtrhpa(
	NumericMatrix x,
	NumericMatrix tr_left = NumericMatrix(),
	NumericMatrix tr_right = NumericMatrix(),
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	bool is_parallel = false,
	bool log = false,
	bool is_validation = true)
{  
  int n = x.nrow();
  int m = x.ncol();
    
  if (is_validation)
  {
    if((tr_left.nrow() != tr_right.nrow()) | 
       (tr_left.ncol() != tr_right.ncol()))
    {
      stop("tr_left and tr_right should be matrices of the same dimensions.");
    }
    
    // Insure that all values are between
    // lower and upper truncation points
    if ((tr_left.nrow() == 1) | (tr_right.nrow() == 1))
    {
      for(int i = 0; i < m; i++)
      {
        double tr_left_value = tr_left[i];
        double tr_right_value = tr_right[i];
        if(tr_left_value >= tr_right_value)
        {
          stop("tr_right element's should greater than corresponding tr_left elements");
        }
        for(int j = 0; j < n; j++)
        {
          if((x(j, i) < tr_left_value) |
             (x(j, i) > tr_right_value))
          {
            NumericVector tr_return = rep(0.0, n);
            if(log)
            {
              std::fill(tr_return.begin(), tr_return.end(), R_NegInf);
            }
            return(tr_return);
          }
        }
      }
    } else {
      for(int i = 0; i < m; i++)
      {
        for(int j = 0; j < n; j++)
        {
          if((x(j, i) < tr_left(j, i)) | 
             (x(j, i) > tr_right(j, i)))
          {
            NumericVector tr_return = rep(0.0, n);
            if(log)
            {
              std::fill(tr_return.begin(), tr_return.end(), R_NegInf);
            }
            return(tr_return);
          }
        }
      }
    }
  }

	// Calculate the nominator
	NumericVector density_main = dhpa(
		x,
		pol_coefficients, pol_degrees,
		given_ind, omit_ind,
		mean, sd,
		is_parallel,
		log, is_validation);

	// Calculate the denominator
	NumericVector cdf_tr = ihpa( 
	  tr_left, tr_right,
		pol_coefficients, pol_degrees,
		given_ind, omit_ind, 
		mean, sd, 
	  is_parallel,
	  log, is_validation);
	
	NumericVector return_value;

	if ((tr_left.nrow() == 1) | (tr_right.nrow() == 1))
	{
	  if(log)
	  {
	    return_value = density_main - cdf_tr[0];
	  } else {
	    return_value = density_main / cdf_tr[0];
	  }
	  
	  return(return_value);
	} 

	if(log)
	{
	  return_value = density_main - cdf_tr;
	} else {
	  return_value = density_main / cdf_tr;
	}
	
	return(return_value);
}

//' @name hpaDist
//' @template itrhpa_examples_Template
//' @export
// [[Rcpp::export]]
NumericVector itrhpa(
	NumericMatrix x_lower = NumericMatrix(1, 1),
	NumericMatrix x_upper = NumericMatrix(1, 1),
	NumericMatrix tr_left = NumericMatrix(1, 1),
	NumericMatrix tr_right = NumericMatrix(1, 1),
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	bool is_parallel = false,
	bool log = false,
	bool is_validation = true)
{
  int n = x_upper.nrow();
  int m = x_upper.ncol();
  
  if (is_validation)
  {
    for(int i = 0; i < m; i++)
    {
      for(int j = 0; j < n; j++)
      {
        if(x_lower(j, i) >= x_upper(j, i))
        {
          stop("x_lower elements should be less than corresponding x_upper elements");
        }
      }
    }
    
    if((tr_left.nrow() != tr_right.nrow()) | 
       (tr_left.ncol() != tr_right.ncol()))
    {
      stop("tr_left and tr_right should be matrices of the same dimensions");
    }
    
    // Set values under and above truncation points
    // equal to truncation points
    if ((tr_left.nrow() == 1) | (tr_right.nrow() == 1))
    {
      for(int i = 0; i < m; i++)
      {
        double tr_left_value = tr_left[i];
        double tr_right_value = tr_right[i];
        if(tr_left_value >= tr_right_value)
        {
          stop("tr_right element's should greater than corresponding tr_left elements");
        }
        for (int j = 0; j < n; j++)
        {
          if (x_lower(j, i) < tr_left_value)
          {
            x_lower(j, i) = tr_left_value;
          }
          if (x_upper(j, i) > tr_right_value)
          {
            x_upper(j, i) = tr_right_value;
          }
        }
      }
    } else {
      for(int i = 0; i < m; i++)
      {
        for(int j = 0; j < n; j++)
        {
          if (x_lower(j, i) < tr_left(j, i))
          {
            x_lower(j, i) = tr_left(j, i);
          }
          if (x_upper(j, i) > tr_right(j, i))
          {
            x_upper(j, i) = tr_right(j, i);
          }
        }
      }
    }
  }
  
	// Calculate the nominator
	NumericVector interval_main = ihpa(
		x_lower, x_upper,
		pol_coefficients, pol_degrees,
		given_ind, omit_ind,
		mean, sd,
    is_parallel,
    log, is_validation);

	// Calculate the denominator
	NumericVector interval_tr = ihpa(
		tr_left, tr_right,
		pol_coefficients, pol_degrees,
		given_ind, omit_ind,
		mean, sd,
		is_parallel,
		log, is_validation);

	NumericVector return_value;
	
	if ((tr_left.nrow() == 1) | (tr_right.nrow() == 1))
	{
	  if(log)
	  {
	    return_value = interval_main - interval_tr[0];
	  } else {
	    return_value = interval_main / interval_tr[0];
	  }
		return(return_value);
	}
	
	if(log)
	{
	  return_value = interval_main - interval_tr;
	} else {
	  return_value = interval_main / interval_tr;
	}
	
	return(return_value);
}

//' @name hpaDist
//' @template dhpaDiff_examples_Template
//' @export
// [[Rcpp::export]]
NumericMatrix dhpaDiff(
    NumericMatrix x,
    NumericVector pol_coefficients = NumericVector(0),
    NumericVector pol_degrees = NumericVector(0),
    LogicalVector given_ind = LogicalVector(0),
    LogicalVector omit_ind = LogicalVector(0),
    NumericVector mean = NumericVector(0),
    NumericVector sd = NumericVector(0),
    String type = "pol_coefficients",
    bool is_parallel = false,
    bool log = false,
    bool is_validation = true)
{
  List return_List;
  
  NumericMatrix pc_grad;
  NumericMatrix mean_grad;
  NumericMatrix sd_grad;
  NumericMatrix x_upper_grad;
  NumericMatrix all_grad;
  
  StringVector pc_names;
  StringVector mean_names;
  StringVector sd_names;
  StringVector x_upper_names;
  StringVector all_names;
  
  int pol_degrees_n = pol_degrees.size();
  int pol_coefficients_n = pol_coefficients.size();

  if ((type == "pol_coefficients") | (type == "mean") |
      (type == "sd") | (type == "x") |
      (type == "all"))
  {
  return_List = hpaMain(
    NumericMatrix(1, 1),            // x_lower
    x,                              // x_upper
    pol_coefficients, pol_degrees,
    "pdf",                          // type
    given_ind, omit_ind,
    mean, sd,
    NumericVector(0),
    type,                           // grad type
    is_parallel,
    false, log,
    is_validation);                          
  } else {
    stop("Argument type should be 'pol_coefficients', 'mean', 'sd', 'x' or 'all'.");
  }

  // Prepare output and assign the names
  
    // to polynomial coefficients
  if((type == "pol_coefficients") | (type == "all"))
  {
    NumericMatrix pol_ind = polynomialIndex(pol_degrees, false);
    
    pc_names = StringVector(pol_coefficients_n);
    
    for (int i = 0; i < pol_coefficients_n; i++)
    {
      pc_names[i] = "a";
      for (int j = 0; j < pol_degrees_n; j++)
      {
        int my_int = pol_ind(j, i);
        pc_names[i] = as<std::string>(pc_names[i]) + 
                      "_" + std::to_string(my_int);
      }
    }
    
    NumericMatrix tmp_grad = return_List["pc_grad"];
    pc_grad = tmp_grad;
    colnames(pc_grad) = pc_names;
    if(type != "all")
    {
      return(pc_grad);
    }
  }
  
    // to mean
  if((type == "mean") | (type == "all"))
  {
    mean_names = StringVector(pol_degrees_n);
    
    for (int i = 0; i < pol_degrees_n; i++)
    {
      mean_names[i] = "mean_" + std::to_string(i + 1);
    }
    
    NumericMatrix tmp_grad = return_List["mean_grad"];
    mean_grad = tmp_grad;
    colnames(mean_grad) = mean_names;
    if(type != "all")
    {
      return(mean_grad);
    }
  }
  
    // to sd
  if((type == "sd") | (type == "all"))
  {
    sd_names = StringVector(pol_degrees_n);
    
    for (int i = 0; i < pol_degrees_n; i++)
    {
      sd_names[i] = "sd_" + std::to_string(i + 1);
    }
    
    NumericMatrix tmp_grad = return_List["sd_grad"];
    sd_grad = tmp_grad;
    colnames(sd_grad) = sd_names;
    if(type != "all")
    {
      return(sd_grad);
    }
  }
  
    // to x
  if((type == "x") | (type == "all"))
  {
    x_upper_names = StringVector(pol_degrees_n);
    
    for (int i = 0; i < pol_degrees_n; i++)
    {
      x_upper_names[i] = "x_" + std::to_string(i + 1);
    }
    
    NumericMatrix tmp_grad = return_List["x_upper_grad"];
    x_upper_grad = tmp_grad;
    colnames(x_upper_grad) = x_upper_names;
    if(type != "all")
    {
      return(x_upper_grad);
    }
  }
  
  if(type == "all")
  {
    int n = pc_grad.nrow();
    int all_grad_n = pol_coefficients_n + 3 * pol_degrees_n;
    
    all_names = StringVector(all_grad_n);
    all_grad = NumericMatrix(n, all_grad_n);

    for(int i = 0; i < pol_coefficients_n; i++)
    {
      all_grad(_, i) = pc_grad(_, i);
      all_names[i] = pc_names[i];
    }

    for(int i = pol_coefficients_n; 
        i < pol_coefficients_n + pol_degrees_n; i++)
    {
      all_grad(_, i) = mean_grad(_, i - pol_coefficients_n);
      all_names[i] = mean_names[i - pol_coefficients_n];
    }

    for(int i = pol_coefficients_n + pol_degrees_n; 
        i < pol_coefficients_n + 2 * pol_degrees_n; i++)
    {
      all_grad(_, i) = sd_grad(_, i - pol_coefficients_n - pol_degrees_n);
      all_names[i] = sd_names[i - pol_coefficients_n - pol_degrees_n];
    }

    for(int i = pol_coefficients_n + 2 * pol_degrees_n; 
        i < pol_coefficients_n + 3 * pol_degrees_n; i++)
    {
      all_grad(_, i) = x_upper_grad(_, i - pol_coefficients_n - 2 * pol_degrees_n);
      all_names[i] = x_upper_names[i - pol_coefficients_n - 2 * pol_degrees_n];
    }
    
    colnames(all_grad) = all_names;
  }
  
  return(all_grad);
}

//' @name hpaDist
//' @template ihpaDiff_examples_Template
//' @export
// [[Rcpp::export]]
NumericMatrix ihpaDiff(
    NumericMatrix x_lower = NumericMatrix(1, 1),
    NumericMatrix x_upper = NumericMatrix(1, 1),
    NumericVector pol_coefficients = NumericVector(0),
    NumericVector pol_degrees = NumericVector(0),
    LogicalVector given_ind = LogicalVector(0),
    LogicalVector omit_ind = LogicalVector(0),
    NumericVector mean = NumericVector(0),
    NumericVector sd = NumericVector(0),
    String type = "pol_coefficients",
    bool is_parallel = false,
    bool log = false,
    bool is_validation = true)
{
  List return_List;
  
  NumericMatrix pc_grad;
  NumericMatrix mean_grad;
  NumericMatrix sd_grad;
  NumericMatrix x_lower_grad;
  NumericMatrix x_upper_grad;
  NumericMatrix all_grad;
  
  StringVector pc_names;
  StringVector mean_names;
  StringVector sd_names;
  StringVector x_lower_names;
  StringVector x_upper_names;
  StringVector all_names;
  
  int pol_degrees_n = pol_degrees.size();
  int pol_coefficients_n = pol_coefficients.size();
  
  if ((type == "pol_coefficients") | (type == "mean") |
      (type == "sd") | (type == "x_lower") | (type == "x_upper") |
      (type == "all"))
  {
    return_List = hpaMain(
      x_lower,                        // x_lower
      x_upper,                        // x_upper
      pol_coefficients, pol_degrees,
      "interval",                     // type
      given_ind, omit_ind,
      mean, sd,
      NumericVector(0),
      type,                           // grad type
      is_parallel,
      false, log,
      is_validation);                          
  } else {
    stop("Argument type should be 'pol_coefficients', 'mean', 'sd', 'x_lower', 'x_upper' or 'all'.");
  }
  
  // Prepare output and assign the names
  
    // to polynomial coefficients
  if((type == "pol_coefficients") | (type == "all"))
  {
    NumericMatrix pol_ind = polynomialIndex(pol_degrees, false);
    
    pc_names = StringVector(pol_coefficients_n);
    
    for (int i = 0; i < pol_coefficients_n; i++)
    {
      pc_names[i] = "a";
      for (int j = 0; j < pol_degrees_n; j++)
      {
        int my_int = pol_ind(j, i);
        pc_names[i] = as<std::string>(pc_names[i]) + 
          "_" + std::to_string(my_int);
      }
    }
    
    NumericMatrix tmp_grad = return_List["pc_grad"];
    pc_grad = tmp_grad;
    colnames(pc_grad) = pc_names;
    if(type != "all")
    {
      return(pc_grad);
    }
  }
  
    // to mean
  if((type == "mean") | (type == "all"))
  {
    mean_names = StringVector(pol_degrees_n);
    
    for (int i = 0; i < pol_degrees_n; i++)
    {
      mean_names[i] = "mean_" + std::to_string(i + 1);
    }
    
    NumericMatrix tmp_grad = return_List["mean_grad"];
    mean_grad = tmp_grad;
    colnames(mean_grad) = mean_names;
    if(type != "all")
    {
      return(mean_grad);
    }
  }
  
    // to sd
  if((type == "sd") | (type == "all"))
  {
    sd_names = StringVector(pol_degrees_n);
    
    for (int i = 0; i < pol_degrees_n; i++)
    {
      sd_names[i] = "sd_" + std::to_string(i + 1);
    }
    
    NumericMatrix tmp_grad = return_List["sd_grad"];
    sd_grad = tmp_grad;
    colnames(sd_grad) = sd_names;
    if(type != "all")
    {
      return(sd_grad);
    }
  }
  
    // to x_lower
  if((type == "x_lower") | (type == "all"))
  {
    x_lower_names = StringVector(pol_degrees_n);
    
    for (int i = 0; i < pol_degrees_n; i++)
    {
      x_lower_names[i] = "x_lower_" + std::to_string(i + 1);
    }
    
    NumericMatrix tmp_grad = return_List["x_lower_grad"];
    x_lower_grad = tmp_grad;
    colnames(x_lower_grad) = x_lower_names;
    if(type != "all")
    {
      return(x_lower_grad);
    }
  }
  
    // to x_upper
  if((type == "x_upper") | (type == "all"))
  {
    x_upper_names = StringVector(pol_degrees_n);
    
    for (int i = 0; i < pol_degrees_n; i++)
    {
      x_upper_names[i] = "x_upper_" + std::to_string(i + 1);
    }
    
    NumericMatrix tmp_grad = return_List["x_upper_grad"];
    x_upper_grad = tmp_grad;
    colnames(x_upper_grad) = x_upper_names;
    if(type != "all")
    {
      return(x_upper_grad);
    }
  }
  
  if(type == "all")
  {
    int n = pc_grad.nrow();
    int all_grad_n = pol_coefficients_n + 4 * pol_degrees_n;
    
    all_names = StringVector(all_grad_n);
    all_grad = NumericMatrix(n, all_grad_n);

    for(int i = 0; i < pol_coefficients_n; i++)
    {
      all_grad(_, i) = pc_grad(_, i);
      all_names[i] = pc_names[i];
    }

    for(int i = pol_coefficients_n; 
        i < pol_coefficients_n + pol_degrees_n; i++)
    {
      all_grad(_, i) = mean_grad(_, i - pol_coefficients_n);
      all_names[i] = mean_names[i - pol_coefficients_n];
    }

    for(int i = pol_coefficients_n + pol_degrees_n; 
        i < pol_coefficients_n + 2 * pol_degrees_n; i++)
    {
      all_grad(_, i) = sd_grad(_, i - pol_coefficients_n - pol_degrees_n);
      all_names[i] = sd_names[i - pol_coefficients_n - pol_degrees_n];
    }

    for(int i = pol_coefficients_n + 2 * pol_degrees_n; 
        i < pol_coefficients_n + 3 * pol_degrees_n; i++)
    {
      all_grad(_, i) = x_lower_grad(_, i - pol_coefficients_n - 
                                       2 * pol_degrees_n);
      all_names[i] = x_lower_names[i - pol_coefficients_n - 2 * pol_degrees_n];
    }

    for(int i = pol_coefficients_n + 3 * pol_degrees_n; 
        i < pol_coefficients_n + 4 * pol_degrees_n; i++)
    {
      all_grad(_, i) = x_upper_grad(_, i - pol_coefficients_n - 
                                       3 * pol_degrees_n);
      all_names[i] = x_upper_names[i - pol_coefficients_n - 3 * pol_degrees_n];
    }

    colnames(all_grad) = all_names;
  }
  
  return(all_grad);
}
