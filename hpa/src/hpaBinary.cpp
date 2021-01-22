#include "hpaMain.h"
#include "hpaML.h"
#include "hpaBinary.h"
#include "polynomialIndex.h"
#include "hpaValidation.h"
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace RcppParallel;

// [[Rcpp::depends(RcppArmadillo)]]

//' Semi-nonparametric single index binary choice model estimation
//' @description This function performs semi-nonparametric (SNP) maximum 
//' likelihood estimation of single index binary choice model 
//' using Hermite polynomial based approximating function proposed by Gallant 
//' and Nychka in 1987. Please, see \code{\link[hpa]{dhpa}} 'Details' section to 
//' get more information concerning this approximating function.
//' @template formula_Template
//' @template data_Template
//' @template K_Template
//' @template z_mean_fixed_Template
//' @template z_sd_fixed_Template
//' @template z_constant_fixed_Template
//' @template z_coef_first_fixed_Template
//' @template x0_binary_Template
//' @template cov_type_Template
//' @template boot_iter_Template
//' @template is_parallel_Template
//' @template opt_type_Template
//' @template opt_control_Template
//' @template is_validation_Template
//' @param is_x0_probit logical; if \code{TRUE} (default) then initial
//' points for optimization routine will be
//' obtained by probit model estimated via \link[stats]{glm} function.
//' @template is_sequence_Template
//' @template GN_details_Template
//' @template hpaBinary_formula_Template
//' @template is_numeric_Template
//' @template parametric_paradigm_Template
//' @template optim_details_Template
//' @template opt_control_details_Template
//' @template opt_control_details_hpaBinary_Template
//' @return This function returns an object of class "hpaBinary".\cr \cr
//' An object of class "hpaBinary" is a list containing the 
//' following components:
//' \itemize{
//' \item \code{optim} - \code{\link[stats]{optim}} function output. 
//' If \code{opt_type = "GA"} then it is the list containing 
//' \code{\link[stats]{optim}} and \code{\link[GA]{ga}} functions outputs.
//' \item \code{x1} - numeric vector of distribution parameters estimates.
//' \item \code{mean} - mean (mu) parameter of density function estimate.
//' \item \code{sd} - sd (sigma) parameter of density function estimate.
//' \item \code{pol_coefficients} - polynomial coefficients estimates.
//' \item \code{pol_degrees} - the same as \code{K} input parameter.
//' \item \code{coefficients} - regression (single index) 
//' coefficients estimates.
//' \item \code{cov_mat} - covariance matrix estimate.
//' \item \code{marginal_effects} - marginal effects matrix where columns are
//' variables and rows are observations.
//' \item \code{results} - numeric matrix representing estimation results.
//' \item \code{log-likelihood} - value of Log-Likelihood function.
//' \item \code{AIC} - AIC value.
//' \item \code{errors_exp} - random error expectation estimate.
//' \item \code{errors_var} - random error variance estimate.
//' \item \code{dataframe} - data frame containing variables mentioned in 
//' \code{formula} without \code{NA} values.
//' \item \code{model_Lists} - lists containing information about 
//' fixed parameters and parameters indexes in \code{x1}.
//' \item \code{n_obs} - number of observations.
//' \item \code{z_latent} - latent variable (single index) estimates.
//' \item \code{z_prob} - probabilities of positive 
//' outcome (i.e. 1) estimates.}
//' @seealso \link[hpa]{summary.hpaBinary}, \link[hpa]{predict.hpaBinary}, 
//' \link[hpa]{plot.hpaBinary},
//' \link[hpa]{logLik.hpaBinary}
//' @template hpaBinary_examples_Template
//' @export	
// [[Rcpp::export]]
List hpaBinary(Rcpp::Formula formula,
	DataFrame data,
	int K = 1,
	double z_mean_fixed = NA_REAL,
	double z_sd_fixed = NA_REAL,
	double z_constant_fixed = 0,
	bool is_z_coef_first_fixed = true,
	bool is_x0_probit = true,
	bool is_sequence = false,
	NumericVector x0 = NumericVector(0),
	String cov_type = "sandwich",
	int boot_iter = 100,
	bool is_parallel = false,
	String opt_type = "optim",
	List opt_control = R_NilValue,
	bool is_validation = true)
{
  // Validation
  
  if (is_validation)
  {
      // Check covariance matrix type
    if ((cov_type != "sandwich") & (cov_type != "sandwichFD") &
        (cov_type != "bootstrap") & (cov_type != "gop") & 
        (cov_type != "hessian") & (cov_type != "hessianFD"))
    {
      stop("Incorrect cov_type argument value.");
    }
    
      // Check opt_type
      if ((opt_type != "optim") & (opt_type != "GA"))
      {
        stop("Incorrect opt_type argument value.");
      }
      
      // Warning concerning large number of bootstrap iterations
      if (boot_iter > 1000)
      {
        warning("Since boot_iter is large estimation may take lots of time.");
      }
      
      // Validate polynomial degree
      pol_Validate(NumericVector::create(K), NumericVector(0));
  }

	// Load additional environments
	
	  // stats environment
	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function optim = stats_env["optim"];
	Rcpp::Function model_frame = stats_env["model.frame"];
	Rcpp::Function na_omit_R = stats_env["na.omit"];
	Rcpp::Function glm = stats_env["glm"];
	Rcpp::Function binomial = stats_env["binomial"];
	Rcpp::Function na_pass = stats_env["na.pass"];
	Rcpp::Function cov_R = stats_env["cov"];
	
	  // base environment
	Rcpp::Environment base_env("package:base");
	Rcpp::Function class_R = base_env["class"];
	Rcpp::Function c_R = base_env["c"];
	Rcpp::Function diag_R = base_env["diag"];
	Rcpp::Function requireNamespace_R = base_env["requireNamespace"];
	Rcpp::Function cat_R = base_env["cat"];
	
	  // GA environment
	Rcpp::Function ga_R = stats_env["optim"];
	Rcpp::Function ga_summary_R = stats_env["optim"];
	
	    // Should be LogicalVector not bool otherwise not working
	LogicalVector is_ga_installed = requireNamespace_R(Rcpp::_["package"] = "GA", 
                                                     Rcpp::_["quietly"] = true);
	
	    // Check weather GA package has been installed
	if(opt_type == "GA")
	{
	  if(is_ga_installed[0])
	  {
	    Rcpp::Environment ga_env = Rcpp::Environment::namespace_env("GA");
	    ga_R = ga_env["ga"];
	    ga_summary_R = ga_env["summary.ga"];
	  } else {
	    stop("Package GA needed if (opt_type = 'GA'). Please install it.");
	  }
	}
	
	// Initialize polynomial structure related values

	int pol_coefficients_n = K;

	// Initialize bool values related to fixed parameters
	bool is_z_mean_fixed = !R_IsNA(z_mean_fixed);
	bool is_z_sd_fixed = !R_IsNA(z_sd_fixed);
	bool is_z_constant_fixed = !R_IsNA(z_constant_fixed);

	// Working with Data

		// Extract data frame from formula
	DataFrame z_df = model_frame(Rcpp::_["formula"] = formula, 
                               Rcpp::_["data"] = data,
                               Rcpp::_["na.action"] = na_pass);

	z_df = na_omit_R(z_df);
	int z_df_n = z_df.size();

		// Extract binary dependent variable
	NumericVector z = z_df[0];
	
	  // Get observations number
	int n_obs = z.size();

		// Extract independent variables (regressors)
	NumericMatrix z_d(n_obs, (z_df_n - 1) + 
	                  !is_z_constant_fixed); // -1 because of dependent variable
	int z_d_col = z_d.ncol();                // number of independent variables

	    // The constant located in last column of regressors matrix
	if (!is_z_constant_fixed)
	{
		z_d(_, z_d_col -  1) = (NumericVector(n_obs) + 1); // add constant
	}

	    // Provide variables from data frame to regressors matrix
	for (int i = 0; i < (z_d_col - !is_z_constant_fixed); i++)
	{
		z_d(_, i) = as<NumericVector>(z_df[i + 1]); // +1 because of 
	}                                             // dependent variable

	// Sequential estimation
	if (is_sequence)
	{
		// for K=0
		List results(K + 1);
		results[0] = hpaBinary(formula, data, 
                         0, z_mean_fixed, z_sd_fixed, z_constant_fixed, 
                         is_z_coef_first_fixed, true, false, NumericVector(0), 
                         "sandwich", 100, is_parallel, 
                         opt_type, opt_control, false);
		List results_current = results[0];
		NumericVector x1 = results_current["x1"];
		int x0_n = x1.size() + 1; // add one more alpha parameter 
		                          // for the next estimation
		NumericVector x0 = NumericVector(x0_n);
		for (int i = 1; i < x0_n; i++)
		{
			x0[i] = x1[i - 1];
		}
		
		// for other K
		for (int i = 1; i <= K; i++)
		{
			if (is_z_sd_fixed)
			{
				z_sd_fixed = results_current["sd"];
			}
			results[i] = hpaBinary(formula, data, i, 
                          z_mean_fixed, z_sd_fixed, z_constant_fixed, 
                          is_z_coef_first_fixed, 
                          false, false, x0, "sandwich", 100, is_parallel, 
                          opt_type, opt_control);
			results_current = results[i];
			x1 = results_current["x1"];
			x0_n++;
			x0 = NumericVector(x0_n);
			for (int j = 0; j < x0_n; j++)
			{
				if (i > j)
				{
					x0[j] = x1[j];
				}
				if (i < j)
				{
					x0[j] = x1[j - 1];
				}
			}
		}
		return(results);
	}

	// Create initial values vector
	bool x0_given = true;

	if (x0.size() == 0)
	{
		x0_given = false;
		// x0 dimensions are equal to the estimated (nonfixed) parameters number
		// note thet constant is already in z_d_col or not
		x0 = NumericVector(pol_coefficients_n +
							         !is_z_mean_fixed + !is_z_sd_fixed +
							         z_d_col - is_z_coef_first_fixed);
	}

	// Assign indexes

		// Initialize additional index and upper value for some loops
	int k = 0; // to account for fixed parameters

	int lower_ind = 0;
	int upper_ind = pol_coefficients_n;

		// for polynomial coefficients
	NumericVector pol_coefficients_ind(pol_coefficients_n);

	if (K != 0)
	{
		for (int i = lower_ind; i < upper_ind; i++)
		{
			pol_coefficients_ind[i] = i;
		}
	} else {
		pol_coefficients_ind = NumericVector::create(0);
	}
	
		// for mean vector
	int z_mean_ind = 0;

	if (!is_z_mean_fixed)
	{
		z_mean_ind = pol_coefficients_n + k;
		k++;
	}

		// for sd vector
	int z_sd_ind = 0;

	if (!is_z_sd_fixed)
	{
		z_sd_ind = pol_coefficients_n + k;
		k++;
	}
	
		// for z coefficients
		// Note that z_d_col may contain or not the constant term
	int n_coef = z_d_col - is_z_coef_first_fixed; // number of 
	                                              // estimated coefficients
	NumericVector z_coef_ind(n_coef);

	lower_ind = pol_coefficients_n + k;
	upper_ind = pol_coefficients_n + z_coef_ind.size() + k;

	for (int i = lower_ind; i < upper_ind; i++)
	{
		z_coef_ind[i - lower_ind] = i;
	}

	// Convert to arma
	arma::vec z_arma = as<arma::vec>(z);
	arma::mat z_d_arma = as<arma::mat>(z_d);

	// Divide into 0 and 1 samples
	arma::vec z_1 = as<arma::vec>(z[z == 1]);
	arma::mat z_d_1 = (as<arma::mat>(z_d)).rows(arma::find(z_arma == 1));

	arma::vec z_0 = as<arma::vec>(z[z == 0]);
	arma::mat z_d_0 = (as<arma::mat>(z_d)).rows(arma::find(z_arma == 0));

	// Set initial sd value to 1
	if (!is_z_sd_fixed & !x0_given)
	{
		x0[z_sd_ind] = 1;
	}

	// Estimate intial values via probit model using glm function
	if (is_x0_probit & !x0_given)
	{
	  // Calculate probit model via glm function
		List model_probit = glm(Rcpp::_["formula"] = formula, 
                            Rcpp::_["data"] = data,
		                        Rcpp::_["family"] = binomial(
		                          Rcpp::_["link"] = "probit"));
	  
	  // Extract probit model coefficients estimates
		NumericVector glm_coef = model_probit["coefficients"];
		
		// Coefficient for the first regressor which under some
		// input parameters should be used for adjust purposes
		double coef_adjust = std::abs(glm_coef[1]);
		
		if (is_z_coef_first_fixed)
		{
		  // Adjust all coefficients
		  glm_coef = glm_coef / coef_adjust;
		  
			// Adjust sd to the first coefficient value
			if (is_z_sd_fixed)
			{
				z_sd_fixed /= coef_adjust;
			} else {
				x0[z_sd_ind] = x0[z_sd_ind] / coef_adjust;
			}
		} else {
		  // Adjust coefficients to sd parameter
			if (is_z_sd_fixed)
			{
				glm_coef = glm_coef / z_sd_fixed;
			}
		}
		if (!is_z_constant_fixed)
		{
			x0[z_coef_ind[n_coef - 1]] = glm_coef[0]; // already adjusted because 
		} else {                                    // glm_coef = glm_coef / coef_adjust 
			if (!is_z_mean_fixed)
			{
				x0[z_mean_ind] = glm_coef[0];           // already adjusted because 
			} else {                                  // glm_coef = glm_coef / coef_adjust 
				z_mean_fixed = glm_coef[0];
			}
		}
		for (int i = 0; i < (n_coef - !is_z_constant_fixed); i++)
		{
			x0[z_coef_ind[i]] = glm_coef[i + 1 +                // + 1 to omit constant 
			                            is_z_coef_first_fixed]; // assigned previously
		}
	}

	// Create list for some variables because unfortunately in Rcpp optim function 
	// has limitation for the input arguments number
	
	    // Collect some values to lists since there are limited 
	    // number of objects could be stored in the list in Rcpp
	List is_List = List::create(
	  Named("is_z_coef_first_fixed") = is_z_coef_first_fixed, 
		Named("is_z_mean_fixed") = is_z_mean_fixed,
		Named("is_z_sd_fixed") = is_z_sd_fixed,
		Named("is_z_constant_fixed") = is_z_constant_fixed);
	
	List ind_List = List::create(
	  Named("z_mean_ind") = z_mean_ind, 
    Named("z_sd_ind") = z_sd_ind,
    Named("z_coef_ind") = z_coef_ind,
    Named("pol_coefficients_ind") = pol_coefficients_ind);

	List fixed_List = List::create(
	  Named("z_mean_fixed") = z_mean_fixed,
	  Named("z_sd_fixed") = z_sd_fixed,
	  Named("z_constant_fixed") = z_constant_fixed);
	
	
	  // Store all the values into the list
	List hpaBinary_args = List::create(Named("is_List") = is_List, 
                                     Named("ind_List") = ind_List,
                                     Named("fixed_List") = fixed_List,
                                     Named("z_1") = z_1,
                                     Named("z_0") = z_0,
                                     Named("z_d_1") = z_d_1,
                                     Named("z_d_0") = z_d_0,
                                     Named("K") = K,
                                     Named("is_parallel") = is_parallel);

	// Apply optimization routine
	
	  // Set optim control parameters
	List PGN_control = List::create(
	     Named("maxit") = 10000000, 
       Named("fnscale") = -1.0,
       Named("abstol") = std::sqrt(
         std::numeric_limits<double>::epsilon()) * 0.01,
       Named("reltol") = std::sqrt(
         std::numeric_limits<double>::epsilon()) * 0.01);
	
	    // Perform the optimization
	List optim_results = optim(
		Rcpp::_["par"] = x0,
		Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaBinaryLnLOptim),
		Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaBinaryLnLOptim_grad),
		Rcpp::_["control"] = PGN_control,
		Rcpp::_["method"] = "BFGS",
		Rcpp::_["hessian"] = true,
		Rcpp::_["hpaBinary_args"] = hpaBinary_args);

	// Extract coefficients and function value
	NumericVector x1 = optim_results["par"];
	
	int x1_n = x0.size();
	
	  // Evolutionary algorithm
	List ga_List;
	List ga_summary;
	
	    // lower and upper bounds for parameters space
	NumericVector ga_lower = NumericVector(x1_n);
	NumericVector ga_upper = NumericVector(x1_n);
	
	if(opt_type == "GA")
	{
	  //suggest initial solution based on local optimization
	  NumericMatrix ga_suggestions;

	  if(opt_control.containsElementNamed("suggestions"))
	  {
	    ga_suggestions = Rcpp::as<Rcpp::NumericMatrix>(opt_control["suggestions"]);
	  } else {
	    ga_suggestions = NumericMatrix(1, x1_n);
	    ga_suggestions(0,_) = x1;
	  }

	  // set lower and upper bounds for parameters space
	  if(opt_control.containsElementNamed("lower") & 
       opt_control.containsElementNamed("upper"))
	  {
	    ga_lower = opt_control["lower"];
	    ga_upper = opt_control["upper"];
	  } else {
  	  // bounds for the mean parameter
  	  
  	  double x1_mean;
  	  
  	  if(!is_z_mean_fixed)
  	  {
    	  x1_mean = x1[z_mean_ind];
    	  
    	  double ga_lower_mean = x1_mean - 2 * abs(x1_mean);
    	  double ga_upper_mean = x1_mean + 2 * abs(x1_mean);
  
    	  ga_lower[z_mean_ind] = ga_lower_mean;
    	  ga_upper[z_mean_ind] = ga_upper_mean;
  	  }
  	  
  	  // bounds for the sd parameter
  	  
  	  if(!is_z_sd_fixed)
  	  {
    	  double ga_lower_sd = x1[z_sd_ind] * 0.2;
    	  double ga_upper_sd = x1[z_sd_ind] * 5;
    	  
    	  ga_lower[z_sd_ind] = ga_lower_sd;
    	  ga_upper[z_sd_ind] = ga_upper_sd;
  	  }
  	  
  	  // bounds for the regression coefficients
  	  for(int i = z_coef_ind[0]; i < (z_coef_ind[0] + z_coef_ind.size()); i++)
  	  {
  	    double z_coef_ga = x1[i];
  	    ga_lower[i] = z_coef_ga - 2 * std::abs(z_coef_ga);
  	    ga_upper[i] = z_coef_ga + 2 * std::abs(z_coef_ga);
  	  }
  	  
  	  // bounds for the polynomial coefficients parameters
  	  ga_lower[pol_coefficients_ind] = -10;
  	  ga_upper[pol_coefficients_ind] = 10;
	  }
	  
	    // set maximum number of iterations
	  int ga_maxiter;
	  if(opt_control.containsElementNamed("maxiter"))
	  {
	    ga_maxiter = opt_control["maxiter"];
	  } else {
	    ga_maxiter = 50 * (1 + K) + 10 * z_d_col;
	  }
	  
	    // set population size
	  int ga_popSize;
	  if(opt_control.containsElementNamed("popSize"))
	  {
	    ga_popSize = opt_control["popSize"];
	  } else {
	    ga_popSize = 10 + (1 + K) * 5 + 2 * z_d_col;
	  }
	  
	    // set mutation probability
	  double ga_pmutation;
	  if(opt_control.containsElementNamed("pmutation"))
	  {
	    ga_pmutation = opt_control["pmutation"];
	  } else {
	    ga_pmutation = 0.2;
	  }
	  
	    // set crossover probability
	  double ga_pcrossover;
	  if(opt_control.containsElementNamed("pcrossover"))
	  {
	    ga_pcrossover = opt_control["pcrossover"];
	  } else {
	    ga_pcrossover = 0.8;
	  }
	  
	    // set elitism parameter
	  int ga_elitism;
	  if(opt_control.containsElementNamed("elitism"))
	  {
	    ga_elitism = opt_control["elitism"];
	  } else {
	    ga_elitism = 5 + (int)(std::round(ga_popSize * 0.1));
	  }
	  
	    // set optim parameter
	  bool ga_optim;
	  if(opt_control.containsElementNamed("optim"))
	  {
	    ga_optim = opt_control["optim"];
	  } else {
	    ga_optim = true;
	  }
	  
	    // set optim args
	  List ga_optimArgs;
	  if(opt_control.containsElementNamed("optimArgs"))
	  {
	    ga_optimArgs = opt_control["optimArgs"];
	  } else {
	    ga_optimArgs = List::create(
	      Named("control") = PGN_control,
	      Named("method") = "Nelder-Mead",
	      Named("poptim") = 0.2);
	  }
	  
	    // set seed parameter
	  int ga_seed;
	  if(opt_control.containsElementNamed("seed"))
	  {
	    ga_seed = opt_control["seed"];
	  } else {
	    ga_seed = 8;
	  }
	  
	     // set monitor parameter
	  bool ga_monitor;
	  if(opt_control.containsElementNamed("monitor"))
	  {
	    ga_monitor = opt_control["monitor"];
	  } else {
	    ga_monitor = true;
	  }

	    // apply genetic algorithm for global optimization
	  ga_List = List::create(Named("GA") = ga_R(
	    Rcpp::_["lower"] = ga_lower,
	    Rcpp::_["upper"] = ga_upper,
	    Rcpp::_["fitness"] = Rcpp::InternalFunction(&hpaBinaryLnLOptim),
	    Rcpp::_["type"] = "real-valued",
	    Rcpp::_["maxiter"] = ga_maxiter,
	    Rcpp::_["popSize"] = ga_popSize,
	    Rcpp::_["pmutation"] = ga_pmutation,
	    Rcpp::_["pcrossover"] = ga_pcrossover,
	    Rcpp::_["elitism"] = ga_elitism,
	    Rcpp::_["optim"] = ga_optim,
	    Rcpp::_["optimArgs"] = ga_optimArgs,
	    Rcpp::_["suggestions"] = ga_suggestions,
	    Rcpp::_["seed"] = ga_seed,
	    Rcpp::_["hpaBinary_args"] = hpaBinary_args,
	    Rcpp::_["monitor"] = ga_monitor));

	  // add \n in order to ensure that the prompt is at the next line
	  if(ga_monitor)
	  {
	    cat_R("\n");
	  }
	  
  	  // get ga function summary object
	  ga_summary = ga_summary_R(ga_List["GA"]);
	  x1 = ga_summary["solution"];
	  
	  // repeat local optimization
	  
	  optim_results = optim(
	    Rcpp::_["par"] = x1,
	    Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaBinaryLnLOptim),
	    Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaBinaryLnLOptim_grad),
	    Rcpp::_["control"] = PGN_control,
	    Rcpp::_["method"] = "BFGS",
	    Rcpp::_["hessian"] = true,
	    Rcpp::_["hpaBinary_args"] = hpaBinary_args);

	  // Reset optimal point
	  x1 = optim_results["par"];
	}
	
	double optim_value = optim_results["value"];

	NumericVector pol_coefficients = NumericVector(K);

	if (K != 0)
	{
		pol_coefficients = x1[pol_coefficients_ind];
	}

	pol_coefficients.push_front(1);

	NumericVector z_mean = NumericVector(1);
	if(!is_z_mean_fixed)
	{ 
		z_mean = NumericVector::create(x1[z_mean_ind]);
	}
	else
	{
		z_mean = NumericVector::create(z_mean_fixed);
	}

	NumericVector z_sd = NumericVector(1);
	if (!is_z_sd_fixed)
	{
		z_sd = NumericVector::create(x1[z_sd_ind]);
	}
	else
	{
		z_sd = NumericVector::create(z_sd_fixed);
	}

	NumericVector z_coef = x1[z_coef_ind];

	// Get covariance matrix estimate of "cov_type" type
	NumericMatrix cov_mat;
	arma::mat H_part;
	arma::mat J_part;
	
	// Estimate Jacobian for the inner part
	if ((cov_type == "gop") | (cov_type == "sandwich") | 
      (cov_type == "sandwichFD"))
	{
	  NumericMatrix my_jacobian = hpaBinaryLnLOptim_grad_ind(x1, hpaBinary_args);
	  J_part = as<arma::mat>(my_jacobian);
	}

	  // Estimate Hessian
	if ((cov_type == "hessian") | (cov_type == "sandwich"))
	{
	  NumericMatrix my_hessian = optim_results["hessian"];
	  H_part = as<arma::mat>(my_hessian).i();
	}

	if ((cov_type == "hessianFD") | (cov_type == "sandwichFD"))
	{
	  NumericMatrix my_hessian = optim_results["hessian"];
	  
	  try
	  {
	    my_hessian = hpaBinaryLnLOptim_hessian(x1, hpaBinary_args);
	  } catch (std::exception &ex) {
	    warning("Can't calculate Hessian via first difference method. Hessian from the optim function will be used instead.");
	    forward_exception_to_r(ex);
	  }
	  
	  H_part = as<arma::mat>(my_hessian).i();
	}

	  // Sandwich Estimate
	if ((cov_type == "sandwich") | (cov_type == "sandwichFD"))
	{
	  cov_mat = wrap(H_part * (J_part.t() * J_part) * H_part);
	}

	  // Inverse Hessian estimate
	if ((cov_type == "hessian") | (cov_type == "hessianFD"))
	{
	  cov_mat = wrap(-H_part);
	}

	  // Gradient (Jacobian) outer product estimate
	if (cov_type == "gop")
	{
	  cov_mat = wrap((J_part.t() * J_part).i());
	}

	// Apply bootstrap
	
	// store parameters for each iteration
	NumericMatrix boot_parameters = NumericMatrix(boot_iter, x1_n);
	
	// store standard deviation
	NumericVector sd_dev = NumericVector(x1_n);
	
	// temporal index matrix for each iteration
	NumericVector sample_ind = NumericVector(n_obs);
	
	// list to store bootstrap results
	List boot_List;
	
	if (cov_type == "bootstrap")
	{
	  for(int i = 0; i < boot_iter; i++)
	  {
	    // Generate sample with replacement
	    NumericVector sample_ind = floor(runif(n_obs, 0, n_obs));
	    
	    NumericVector z_boot = NumericVector(n_obs);
	    NumericMatrix z_d_boot = NumericMatrix(z_d.nrow(), z_d.ncol());
	    
	    for (int j = 0; j < n_obs; j++)
	    {
	      z_boot[j] = z[sample_ind[j]];
	      z_d_boot(j, _) = z_d(sample_ind[j], _);
	    }
	    
	    // Convert to arma
	    arma::vec z_arma_boot = as<arma::vec>(z_boot);
	    arma::mat z_d_arma_boot = as<arma::mat>(z_d_boot);
	    
	    // Divide into 0 and 1 samples
	    arma::vec z_1_boot = as<arma::vec>(z_boot[z_boot == 1]);
	    arma::mat z_d_1_boot = (
	      as<arma::mat>(z_d_boot)).rows(arma::find(z_arma_boot == 1));
	    
	    arma::vec z_0_boot = as<arma::vec>(z_boot[z_boot == 0]);
	    arma::mat z_d_0_boot = (
	      as<arma::mat>(z_d_boot)).rows(arma::find(z_arma_boot == 0));
	    
	    // Prepare arguments for bootstrap List
	    List hpaBinary_args_boot = hpaBinary_args;
	    hpaBinary_args_boot["z_1"] = z_1_boot;
	    hpaBinary_args_boot["z_0"] = z_0_boot;
	    hpaBinary_args_boot["z_d_1"] = z_d_1_boot;
	    hpaBinary_args_boot["z_d_0"] = z_d_0_boot;
	    
	    // Perform estimation
	    List boot_results = optim(
	      Rcpp::_["par"] = x1,
	      Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaBinaryLnLOptim),
	      Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaBinaryLnLOptim_grad),
	      Rcpp::_["control"] = PGN_control,
	      Rcpp::_["method"] = "BFGS",
	      Rcpp::_["hessian"] = true,
	      Rcpp::_["hpaBinary_args"] = hpaBinary_args_boot);
	    
	    // Store iteration results
	    NumericVector x1_new = boot_results["par"];
	    
	    boot_parameters(i, _) = x1_new;
	  }
	  
	  // Store bootstrap results
	  cov_mat = cov_R(boot_parameters);
	  sd_dev = sqrt(diag_R(cov_mat));
	  
	  boot_List = List::create(
	    Named("estimates") = boot_parameters,
	    Named("cov_mat") = cov_mat,
	    Named("sd") = sd_dev);
	}

	// Prepare beautifull results output
	NumericMatrix results(x1_n, 4);
	
	StringVector results_cols = StringVector::create("Estimate", "Std. Error", 
                                                   "z value", "P(>|z|)");
	StringVector results_rows(x1_n);

	NumericMatrix pol_ind = polynomialIndex(NumericVector::create(K), false);

		// for alpha
		for (int i = 0; i < K; i++)
		{
			// Convert double to string
			std::stringstream ss;
			ss << (i + 1); // to start from alpha_1
			std::string str2 = ss.str();
			//
			results_rows[i] = "a_" + str2;
			results(i, 0) = pol_coefficients[i + 1];
			results(i, 1) = sqrt(cov_mat(i, i));
			double z_stat = results(i, 0) / results(i, 1);
			NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
			results(i, 2) = z_stat;
			results(i, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
		}

		// for mean
	if (!is_z_mean_fixed)
	{
		results_rows[z_mean_ind] = "mean";
		results(z_mean_ind, 0) = x1[z_mean_ind];
		results(z_mean_ind, 1) = sqrt(cov_mat((z_mean_ind), (z_mean_ind)));
		double z_stat = results(z_mean_ind, 0) / results(z_mean_ind, 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(z_mean_ind, 2) = z_stat;
		results(z_mean_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

		// for sd
	if (!is_z_sd_fixed)
	{
		results_rows[z_sd_ind] = "sd";
		results(z_sd_ind, 0) = x1[z_sd_ind];
		results(z_sd_ind, 1) = sqrt(cov_mat((z_sd_ind), (z_sd_ind)));
		double z_stat = results(z_sd_ind, 0) / results(z_sd_ind, 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(z_sd_ind, 2) = z_stat;
		results(z_sd_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}
	
		// for z coefficients
	CharacterVector z_df_names = z_df.names();
	String first_coef_name = z_df_names[1];
	z_df_names.erase(0); // remove dependend variable name
	
	// for intercept if need
	if (!is_z_constant_fixed)
	{
		z_df_names.push_back("(Intercept)");
	}
	
	// remove first regressors if its coefficient is fixed
	if (is_z_coef_first_fixed)
	{
	  z_df_names.erase(0);
	}

	int coef_start = x1_n - n_coef;
	for (int i = coef_start; i < x1_n; i++)
	{
		results_rows(i) = z_df_names[i - coef_start]; // +1 because the first is 
		results(i, 0) = z_coef[i - coef_start];       // dependent variable
		results(i, 1) = sqrt(cov_mat(z_coef_ind[i - coef_start], 
                                 z_coef_ind[i - coef_start]));
		double z_stat = results(i, 0) / results(i, 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(i, 2) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
		results(i, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

		// for expectation and variance
	NumericVector z_e = ehpa(NumericMatrix(1, 1), pol_coefficients, 
                           NumericVector::create(K),
		LogicalVector::create(false), LogicalVector::create(false),
		z_mean, z_sd, 
		NumericVector::create(1), false, false);

	NumericVector z_e_2 = ehpa(NumericMatrix(1, 1), pol_coefficients, 
                             NumericVector::create(K),
		LogicalVector::create(0), LogicalVector::create(0),
		z_mean, z_sd,
		NumericVector::create(2), false, false);

		// assign names to the output
	rownames(results) = results_rows;
	colnames(results) = results_cols;

	x1.names() = results_rows;

	rownames(cov_mat) = results_rows;
	colnames(cov_mat) = results_rows;

	if (K != 0)
	{
		pol_coefficients.names() = c_R("a_0", results_rows[pol_coefficients_ind]);
	}

	StringVector z_coef_names = results_rows[z_coef_ind];

	if (is_z_coef_first_fixed)
	{
	  z_coef.push_front(1); // add identity coefficient for fixed value
	  z_coef_names.push_front(first_coef_name);
	}
	
	z_coef.names() = z_coef_names;

	// Estimate latent variable and probabilities

		// coefficients for independend variables
	arma::vec z_coef_arma = as<arma::vec>(z_coef);

		// get estimates for z*
	NumericMatrix z_latent = wrap(z_d_arma * z_coef_arma);

	if (is_z_constant_fixed)
	{
		z_latent = z_latent + z_constant_fixed;
	}

	NumericVector z_prob = 1 - phpa(-1.0 * z_latent, pol_coefficients,
                              		NumericVector::create(K),
                              		LogicalVector::create(0), 
                              		LogicalVector::create(0),
                              		z_mean, z_sd,
                              		is_parallel, false, false);
	
	// Estimate marginal effects
	NumericVector z_den = dhpa(-1.0 * z_latent, pol_coefficients,
                             NumericVector::create(K),
                             LogicalVector::create(0), 
                             LogicalVector::create(0),
                             z_mean, z_sd,
                             is_parallel, false, false);

	int n_coef_total = z_coef.size();
	
	NumericMatrix marginal_effects = NumericMatrix(z_latent.nrow(), n_coef_total);
	
	for (int i = 0; i < n_coef_total; i++)
	{
	  NumericVector me_vec = z_den * z_coef[i];
	  marginal_effects(_, i) = me_vec;
	}

	StringVector n_coefames = z_coef.names();
	colnames(marginal_effects) = n_coefames;

	// return results
	List model_Lists = List::create(Named("is_List") = is_List,
		Named("ind_List") = ind_List,
		Named("fixed_List") = fixed_List);

	if((opt_type == "GA") | (opt_type == "ga"))
	{
	  optim_results = List::create(
	    Named("optim") = optim_results,
	    Named("GA") = ga_List,
	    Named("GA_summary") = ga_summary);
	}
	
	List return_result = List::create(
	  Named("optim") = optim_results,
		Named("x1") = x1,
		Named("mean") = z_mean,
		Named("sd") = z_sd,
		Named("pol_coefficients") = pol_coefficients,
		Named("pol_degrees") = NumericVector::create(K),
		Named("coefficients") = z_coef,
		Named("results") = results,
		Named("errors_exp") = z_e[0],
		Named("errors_var") = (z_e_2[0] - z_e[0] * z_e[0]),
		Named("log-likelihood") = optim_value,
		Named("AIC") = 2 * (x1_n - optim_value),
		Named("n_obs") = z_latent.nrow(),
		Named("z_latent") = z_latent,
		Named("z_prob") = z_prob,
		Named("formula") = formula,
		Named("dataframe") = z_df,
		Named("model_Lists") = model_Lists,
		Named("cov_mat") = cov_mat,
		Named("marginal_effects") = marginal_effects);

	return_result.attr("class") = "hpaBinary";
	
	return(return_result);
}

// Perform semi-nonparametric log-likelihood 
// function estimation for binary choice model
List hpaBinaryLnLOptim_List(NumericVector x0, List hpaBinary_args) 
{
  // Get values from the hpaBinary_args list
  List is_List = Rcpp::as<Rcpp::List>(hpaBinary_args["is_List"]);
  List ind_List = Rcpp::as<Rcpp::List>(hpaBinary_args["ind_List"]);
  List fixed_List = Rcpp::as<Rcpp::List>(hpaBinary_args["fixed_List"]);
  arma::mat z_d_1 = hpaBinary_args["z_d_1"];
  arma::mat z_d_0 = hpaBinary_args["z_d_0"];
  double K = hpaBinary_args["K"];
  bool is_parallel = hpaBinary_args["is_parallel"];
  
	// Get values from the is_List
	bool is_z_coef_first_fixed = is_List["is_z_coef_first_fixed"];
	bool is_z_mean_fixed = is_List["is_z_mean_fixed"];
	bool is_z_sd_fixed = is_List["is_z_sd_fixed"];
	bool is_z_constant_fixed = is_List["is_z_constant_fixed"];

	// Get values from the ind_List
	double z_mean_ind = ind_List["z_mean_ind"];
	double z_sd_ind = ind_List["z_sd_ind"];
	NumericVector z_coef_ind = ind_List["z_coef_ind"];
	NumericVector pol_coefficients_ind = ind_List["pol_coefficients_ind"];
	
	// Get values from the fixed_List
	double z_mean_fixed = fixed_List["z_mean_fixed"];
	double z_sd_fixed = fixed_List["z_sd_fixed"];
	double z_constant_fixed = fixed_List["z_constant_fixed"];

	// Assign estimated parameters values to corresponding vectors

		// polynomial coefficients and degrees
	NumericVector pol_coefficients = NumericVector(K);

	if (K != 0)
	{
		pol_coefficients = x0[pol_coefficients_ind];
	} 

	pol_coefficients.push_front(1);

	NumericVector pol_degrees = NumericVector(1);

	pol_degrees[0] = K;

		// mean value
	NumericVector z_mean = NumericVector(1);

	if (is_z_mean_fixed)
	{
		z_mean[0] = z_mean_fixed;
	}
	else {
		z_mean[0] = x0[z_mean_ind];
	}

		// sd value
	NumericVector z_sd = NumericVector(1);

	if (is_z_sd_fixed)
	{
		z_sd[0] = z_sd_fixed;
	}
	else {
		z_sd[0] = x0[z_sd_ind];
	}

		// coefficients for independend variables
	NumericVector z_coef_R = x0[z_coef_ind];

	if (is_z_coef_first_fixed)
	{
		z_coef_R.push_front(1); // add identity coefficient for fixed value
	}

	arma::vec z_coef = as<arma::vec>(z_coef_R);

	// get estimates for z*
	NumericMatrix z_h_1 = wrap(z_d_1 * z_coef);
	NumericMatrix z_h_0 = wrap(z_d_0 * z_coef);

	  if (is_z_constant_fixed)
  	{
  		z_h_1 = z_h_1 + z_constant_fixed;
	    z_h_0 = z_h_0 + z_constant_fixed;
  	}
	  
	// get number of 1 observations
	int n_obs_1 = z_h_1.nrow();

	// Likelihood calculation
	NumericMatrix inf_vec = NumericMatrix(n_obs_1, 1);
	std::fill(inf_vec.begin(), inf_vec.end(), R_PosInf);
	
	NumericVector lnL_z_1 = ihpa(-1.0 * z_h_1, inf_vec,
                               pol_coefficients, pol_degrees,
                               LogicalVector{false}, LogicalVector{false},
                               z_mean, z_sd, 
                               is_parallel, true, false);
	
	NumericVector lnL_z_0 = phpa(-1.0 * z_h_0,
                               pol_coefficients, pol_degrees,
                               LogicalVector{false}, LogicalVector{false},
                               z_mean, z_sd, 
                               is_parallel, true, false);

	// Initialize list to store calculation results
	double aggregate_0 = 0.0;
	double aggregate_1 = 0.0;
	
	List return_List = List::create(
	  Named("individual_1") = NumericVector::create(0.0),
    Named("individual_0") = NumericVector::create(0.0),
    Named("aggregate_1") = aggregate_1,
    Named("aggregate_0") = aggregate_0);
	
	// Store calculation results
	return_List["individual_1"] = lnL_z_1;
	aggregate_1 = sum(lnL_z_1);
	return_List["aggregate_1"] = aggregate_1;

	return_List["individual_0"] = lnL_z_0;
	aggregate_0 = sum(lnL_z_0);
	return_List["aggregate_0"] = aggregate_0;
	
	return(return_List);
}

// Perform semi-nonparametric log-likelihood 
// function estimation for binary choice model
double hpaBinaryLnLOptim(NumericVector x0, List hpaBinary_args) 
{
  List return_List = hpaBinaryLnLOptim_List(x0, hpaBinary_args);

  double aggregate_0 = return_List["aggregate_0"];
  double aggregate_1 = return_List["aggregate_1"];
  
  double return_aggregate = aggregate_0 + aggregate_1;

  if(any(is_na(NumericVector::create(return_aggregate))))
  {
    return_aggregate = R_NegInf;
  }
  return(return_aggregate);
}

// Perform semi-nonparametric log-likelihood 
// function estimation for binary choice model
NumericVector hpaBinaryLnLOptim_ind(NumericVector x0, List hpaBinary_args) 
{ 
   List return_List = hpaBinaryLnLOptim_List(x0, hpaBinary_args);
   
   NumericVector individual_0 = return_List["individual_0"];
   NumericVector individual_1 = return_List["individual_1"];
   
   int n_obs_0 = individual_0.size();
   int n_obs_1 = individual_1.size();
   int n_obs = n_obs_0 + n_obs_1;
   
   NumericVector return_individual = NumericVector(n_obs);
  
  return_individual[Range(0, n_obs_1 - 1)] = individual_1;
  return_individual[Range(n_obs_1, n_obs - 1)] = individual_0;

  return(return_individual);
}                                          

List hpaBinaryLnLOptim_grad_List(NumericVector x0, List hpaBinary_args) 
{
  // Get values from the hpaBinary_args list
  List is_List = Rcpp::as<Rcpp::List>(hpaBinary_args["is_List"]);
  List ind_List = Rcpp::as<Rcpp::List>(hpaBinary_args["ind_List"]);
  List fixed_List = Rcpp::as<Rcpp::List>(hpaBinary_args["fixed_List"]);
  arma::vec z_1 = hpaBinary_args["z_1"];
  arma::vec z_0 = hpaBinary_args["z_0"];
  arma::mat z_d_1 = hpaBinary_args["z_d_1"];
  arma::mat z_d_0 = hpaBinary_args["z_d_0"];
  double K = hpaBinary_args["K"];
  bool is_parallel = hpaBinary_args["is_parallel"];
  
  // Get values from the is_List
  bool is_z_coef_first_fixed = is_List["is_z_coef_first_fixed"];
  bool is_z_mean_fixed = is_List["is_z_mean_fixed"];
  bool is_z_sd_fixed = is_List["is_z_sd_fixed"];
  bool is_z_constant_fixed = is_List["is_z_constant_fixed"];
  
  // Get values from the ind_List
  double z_mean_ind = ind_List["z_mean_ind"];
  double z_sd_ind = ind_List["z_sd_ind"];
  NumericVector z_coef_ind = ind_List["z_coef_ind"];
  NumericVector pol_coefficients_ind = ind_List["pol_coefficients_ind"];
  
  // Get values from the fixed_List
  double z_mean_fixed = fixed_List["z_mean_fixed"];
  double z_sd_fixed = fixed_List["z_sd_fixed"];
  double z_constant_fixed = fixed_List["z_constant_fixed"];
  
  // Get parameters number
  int n_param = x0.size();
  
  // Get estimated regressors number
  int n_reg = z_coef_ind.size();
  
  // Assign estimated parameters values to corresponding vectors
  
    // polynomial coefficients and degrees
  NumericVector pol_coefficients = NumericVector(K);
  
  if (K != 0)
  {
    pol_coefficients = x0[pol_coefficients_ind];
  } 
  
  pol_coefficients.push_front(1);
  
  NumericVector pol_degrees = NumericVector(1);
  
  pol_degrees[0] = K;
  
    // mean value
  NumericVector z_mean = NumericVector(1);
  
  if (is_z_mean_fixed)
  {
    z_mean[0] = z_mean_fixed;
  }
  else {
    z_mean[0] = x0[z_mean_ind];
  }
  
    // sd value
  NumericVector z_sd = NumericVector(1);
  
  if (is_z_sd_fixed)
  {
    z_sd[0] = z_sd_fixed;
  }
  else {
    z_sd[0] = x0[z_sd_ind];
  }
  
    // coefficients for independend variables
  NumericVector z_coef_R = x0[z_coef_ind];
  
  if (is_z_coef_first_fixed)
  {
    z_coef_R.push_front(1); // add identity coefficient for fixed value
  }
  
  arma::vec z_coef = as<arma::vec>(z_coef_R);
  
  // get estimates for z*
  NumericMatrix z_h_1 = wrap(z_d_1 * z_coef);
  NumericMatrix z_h_0 = wrap(z_d_0 * z_coef);
  
  if (is_z_constant_fixed)
  {
    z_h_1 = z_h_1 + z_constant_fixed;
    z_h_0 = z_h_0 + z_constant_fixed;
  }
  
  // calculate observations numbers
  int n_obs_1 = z_h_1.nrow();
  int n_obs_0 = z_h_0.nrow();
  int n_obs = n_obs_0 + n_obs_1;
  
  // Initialize vector to store gradient values
  NumericMatrix my_grad = NumericMatrix(n_obs_0 + n_obs_1, n_param);
  
  // Analytical part of gradient

    // for polynomial coefficients
  int pol_coefficients_n = pol_coefficients.size();
  
  // lower tail negative infinity matrix
  NumericMatrix neg_inf_vec_0 = NumericMatrix(n_obs_0, 1);
  std::fill(neg_inf_vec_0.begin(), neg_inf_vec_0.end(), R_NegInf);
  
  NumericMatrix inf_vec_1 = NumericMatrix(n_obs_1, 1);
  std::fill(inf_vec_1.begin(), inf_vec_1.end(), R_PosInf);
  
  // Calculate gradient

  NumericMatrix all_grad_1 = ihpaDiff(-1.0 * z_h_1, inf_vec_1,
                                      pol_coefficients, pol_degrees,
                                      LogicalVector{false}, 
                                      LogicalVector{false},
                                      z_mean, z_sd,
                                      "all",
                                      is_parallel, true, false);

  NumericMatrix all_grad_0 = ihpaDiff(neg_inf_vec_0, -1.0 * z_h_0,
                                      pol_coefficients, pol_degrees,
                                      LogicalVector{false}, 
                                      LogicalVector{false},
                                      z_mean, z_sd,
                                      "all",
                                      is_parallel, true, false);

  // Store gradients respect to
  
    // polynomial coefficients
  for (int i = 0; i < (pol_coefficients_n - 1); i++) // for each parameter
  {
    NumericVector my_grad_tmp = NumericVector(n_obs);

    my_grad_tmp[Range(0, n_obs_1 - 1)] = all_grad_1(_, i + 1);
    my_grad_tmp[Range(n_obs_1, n_obs - 1)] = all_grad_0(_, i + 1);

    my_grad(_, i) = my_grad_tmp;
  }
  
    // mean
    if(!is_z_mean_fixed)
    {
      NumericVector my_grad_tmp_mean = NumericVector(n_obs);
      
      my_grad_tmp_mean[Range(0, n_obs_1 - 1)] = all_grad_1(_, pol_coefficients_n);
      my_grad_tmp_mean[Range(n_obs_1, n_obs - 1)] = all_grad_0(_, pol_coefficients_n);
      
      my_grad(_, z_mean_ind) = my_grad_tmp_mean;
    }
    
    // sd
    if(!is_z_sd_fixed)
    {
      NumericVector my_grad_tmp_sd = NumericVector(n_obs);
      
      my_grad_tmp_sd[Range(0, n_obs_1 - 1)] = all_grad_1(_, pol_coefficients_n + 1);
      my_grad_tmp_sd[Range(n_obs_1, n_obs - 1)] = all_grad_0(_, pol_coefficients_n + 1);
      
      my_grad(_, z_sd_ind) = my_grad_tmp_sd;
    }

    // regression coefficients
  NumericMatrix z_d_0_adj = wrap(-1 * z_d_0);
  NumericMatrix z_d_1_adj = wrap(-1 * z_d_1);

  for (int i = 0; i < n_reg; i++)
  {
    NumericVector my_grad_tmp = NumericVector(n_obs);
    
    my_grad_tmp[Range(0, n_obs_1 - 1)] = z_d_1_adj(_, i + is_z_coef_first_fixed) *
                                         all_grad_1(_, pol_coefficients_n + 2);
    my_grad_tmp[Range(n_obs_1, n_obs - 1)] = z_d_0_adj(_, i + is_z_coef_first_fixed) *
                                             all_grad_0(_, pol_coefficients_n + 3);
    
    my_grad(_, z_coef_ind[i]) = my_grad_tmp;
  }

  // Return the results
  List return_List = List::create(Named("aggregate") = colSums(my_grad),
                                  Named("individual") = my_grad);
  
  return(return_List);
}

NumericVector hpaBinaryLnLOptim_grad(NumericVector x0, List hpaBinary_args) 
{
  List return_List = hpaBinaryLnLOptim_grad_List(x0, hpaBinary_args);

  NumericVector return_aggregate = return_List["aggregate"];

  if(any(is_na(return_aggregate)))
  {
    std::fill(return_aggregate.begin(), 
              return_aggregate.end(), 
              R_NegInf);
  }

  return(return_aggregate);
}

NumericMatrix hpaBinaryLnLOptim_grad_ind(NumericVector x0, List hpaBinary_args) 
{
  List return_List = hpaBinaryLnLOptim_grad_List(x0, hpaBinary_args);
  
  NumericMatrix return_individual = return_List["individual"];
  
  return(return_individual);
}

// Perform log-likelihood function hessian estimation 
// for Phillips-Gallant-Nychka distribution at point
NumericMatrix hpaBinaryLnLOptim_hessian(NumericVector x0, List hpaBinary_args)
{
  // Get values from the hpaBinary_args list
  List is_List = Rcpp::as<Rcpp::List>(hpaBinary_args["is_List"]);
  List ind_List = Rcpp::as<Rcpp::List>(hpaBinary_args["ind_List"]);
  List fixed_List = Rcpp::as<Rcpp::List>(hpaBinary_args["fixed_List"]);
  arma::vec z_1 = hpaBinary_args["z_1"];
  arma::vec z_0 = hpaBinary_args["z_0"];
  arma::mat z_d_1 = hpaBinary_args["z_d_1"];
  arma::mat z_d_0 = hpaBinary_args["z_d_0"];
  double K = hpaBinary_args["K"];
  
  // Get parameters number
  int n_param = x0.size();
  
  // Initialize vector to store hessian values
  NumericMatrix my_hessian = NumericMatrix(n_param, n_param);
  
  // Get values from the is_List
  bool is_z_coef_first_fixed = is_List["is_z_coef_first_fixed"];
  bool is_z_mean_fixed = is_List["is_z_mean_fixed"];
  bool is_z_sd_fixed = is_List["is_z_sd_fixed"];
  
  // Get values from the ind_List
  double z_mean_ind = ind_List["z_mean_ind"];
  double z_sd_ind = ind_List["z_sd_ind"];
  NumericVector z_coef_ind = ind_List["z_coef_ind"];
  NumericVector pol_coefficients_ind = ind_List["pol_coefficients_ind"];
  
  // Get values from the fixed_List
  double z_mean_fixed = fixed_List["z_mean_fixed"];
  double z_sd_fixed = fixed_List["z_sd_fixed"];
  
  // Assign estimated parameters values to corresponding vectors
  
  // polynomial coefficients and degrees
  NumericVector pol_coefficients = NumericVector(K);
  
  if (K != 0)
  {
    pol_coefficients = x0[pol_coefficients_ind];
  } 
  
  pol_coefficients.push_front(1);
  
  NumericVector pol_degrees = NumericVector(1);
  
  pol_degrees[0] = K;
  
  // mean value
  NumericVector z_mean = NumericVector(1);
  
  if (is_z_mean_fixed)
  {
    z_mean[0] = z_mean_fixed;
  }
  else {
    z_mean[0] = x0[z_mean_ind];
  }
  
  // sd value
  NumericVector z_sd = NumericVector(1);
  
  if (is_z_sd_fixed)
  {
    z_sd[0] = z_sd_fixed;
  }
  else {
    z_sd[0] = x0[z_sd_ind];
  }
  
  // coefficients for independend variables
  NumericVector z_coef_R = x0[z_coef_ind];
  
  if (is_z_coef_first_fixed)
  {
    z_coef_R.push_front(1); // add identity coefficient for fixed value
  }
  
  // Prepare precision related values for numeric differentiation
  double machinePrecision = std::numeric_limits<double>::epsilon();
  double my_precision = std::sqrt(machinePrecision);
  
  NumericVector eps = abs(x0 * my_precision);
  
  // Control for zero values
  eps[eps < (machinePrecision * 100)] = my_precision;
  
  // Estimate the gradient itself
  NumericVector x0_eps = clone(x0);
  
  // Values to store function values
  // given small increment
  NumericVector g_plus;
  NumericVector g_minus;
  
  // Perform hessian numeric estimation
  for (int i = 0; i < n_param; i++) // for each parameter
  {
    // Calculate g(x + eps)
    x0_eps[i] = x0[i] + eps[i];
    g_plus = hpaBinaryLnLOptim_grad(x0_eps, hpaBinary_args);
    
    // Calculate g(x - eps)
    x0_eps[i] = x0[i] - eps[i];
    g_minus = hpaBinaryLnLOptim_grad(x0_eps, hpaBinary_args);
    
    // Store the results to hessian
    my_hessian(_, i) = (g_plus - g_minus) / (2 * eps[i]);
    
    // Set x0_eps value to default
    x0_eps[i] = x0[i];
  }
  
  // Make hessian to be symmetric
  for (int i = 0; i < n_param; i++)   // for each parameter
  {
    for (int j = i; j < n_param; j++) // for each parameter
    {
      double hessian_val_1 = my_hessian(i, j);
      double hessian_val_2 = my_hessian(j, i);
      my_hessian(i, j) = (hessian_val_1 + hessian_val_2) * 0.5;
      my_hessian(j, i) = my_hessian(i, j);
    }
  }
  
  return(my_hessian);
}

//' Predict method for hpaBinary
//' @param object Object of class "hpaBinary"
//' @template newdata_Template
//' @param is_prob logical; if TRUE (default) 
//' then function returns predicted probabilities. 
//' Otherwise latent variable
//' (single index) estimates will be returned.
//' @return This function returns predicted probabilities 
//' based on \code{\link[hpa]{hpaBinary}} estimation results.
//' @export
// [[Rcpp::export]]
NumericVector predict_hpaBinary(List object, DataFrame newdata = R_NilValue, 
                                bool is_prob = true)
{
	List model = object;

	// Add additional environments
	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function model_frame = stats_env["model.frame"];
	Rcpp::Function na_omit_R = stats_env["na.omit"];

	Rcpp::Environment base_env("package:base");
	Rcpp::Function as_data_frame = base_env["as.data.frame"];

	// Extract variables from model

		// extract is values
	List model_Lists = model["model_Lists"];
	List is_List = model_Lists["is_List"];
	bool is_z_constant_fixed = is_List["is_z_constant_fixed"];

		// extract fixed values
	List fixed_List = model_Lists["fixed_List"];
	double z_constant_fixed = fixed_List["z_constant_fixed"];

		// extract coefficients
	NumericVector pol_coefficients = model["pol_coefficients"];
	NumericVector z_mean = model["mean"];
	NumericVector z_sd = model["sd"];
	NumericVector z_coef = model["coefficients"];

		// extract polynomial coefficients
	double K = pol_coefficients.size() - 1;

	// Check wheather new data frame has been supplied
	DataFrame data = newdata;

	if (newdata.size() == 0)
	{
		newdata = as_data_frame(model["dataframe"]);
	}

	// Remove NA values

	data = na_omit_R(newdata);

	// Working with Data

		// Extract data frame from formula
	Formula formula = model["formula"];
	DataFrame z_df = model_frame(Rcpp::_["formula"] = formula, 
                               Rcpp::_["data"] = data);
	int z_df_n = z_df.size();

	// Extract binary dependend variable
	NumericVector z = z_df[0]; // it is reference
	int n_obs = z.size();
	
	// Extract independend variables
	
	NumericMatrix z_d(n_obs, (z_df_n - 1) + 
	                         !is_z_constant_fixed);      // -1 because of 
	int z_d_col = z_d.ncol();                            // dependent variable
	
	  // the constant located in last column of regressors matrix
	if (!is_z_constant_fixed)
	{
	  z_d(_, z_d_col -  1) = (NumericVector(n_obs) + 1); // add constant
	  
	}
	
	for (int i = 0; i < (z_d_col - !is_z_constant_fixed); i++)
	{
	  z_d(_, i) = as<NumericVector>(z_df[i + 1]);        // +1 because of
	}                                                    // dependent variable

		// Convert to arma
	arma::vec z_arma = as<arma::vec>(z);
	arma::mat z_d_arma = as<arma::mat>(z_d);

	// Estimate latent variable and probabilities

		// coefficients for independend variables
	arma::vec z_coef_arma = as<arma::vec>(z_coef);

	// get estimates for z*
	NumericMatrix z_latent = wrap(z_d_arma * z_coef_arma);

	if (is_z_constant_fixed)
	{
		z_latent = z_latent + z_constant_fixed;
	}

	if (!is_prob)
	{
		NumericVector z_latent_vec = z_latent(_, 0);
		return(z_latent_vec);
	}

	NumericVector z_prob = 1 - phpa((-1) * z_latent, pol_coefficients,
		NumericVector::create(K),
		LogicalVector::create(0), LogicalVector::create(0),
		z_mean, z_sd, false, false, false);

	return(z_prob);
}

//' Summarizing hpaBinary Fits
//' @param object Object of class "hpaBinary"
//' @return This function returns the same list as 
//' \code{\link[hpa]{hpaBinary}} function changing 
//' its class to "summary.hpaBinary".
//' @export
// [[Rcpp::export]]
List summary_hpaBinary(List object)
{

	List return_result = clone(object); // in order to preserve model class

	return_result.attr("class") = "summary.hpaBinary";

	return(return_result);
}

//' Summary for hpaBinary output
//' @param x Object of class "hpaML"
//' @export	
// [[Rcpp::export]]
void print_summary_hpaBinary(List x)
{
	// Extract the model
	List model = x;

	// Load additional environments
	Rcpp::Environment base_env("package:base");
	Rcpp::Function as_table = base_env["as.table"];
	Rcpp::Function cbind = base_env["cbind"];
	Rcpp::Function round_R = base_env["round"];
	Rcpp::Function as_data_frame = base_env["as.data.frame"];
	Rcpp::Function print_R = base_env["print"];
	Rcpp::Function cat_R = base_env["cat"];

	NumericMatrix results = model["results"];
	results = round_R(Rcpp::_["x"] = results, Rcpp::_["digits"] = 5);

	List model_Lists = model["model_Lists"];

	List ind_List = model_Lists["ind_List"];
	NumericVector z_coef_ind = ind_List["z_coef_ind"];

	// Extract is values
	List is_List = model_Lists["is_List"];
	bool is_z_coef_first_fixed = is_List["is_z_coef_first_fixed"];
	bool is_z_mean_fixed = is_List["is_z_mean_fixed"];
	bool is_z_sd_fixed = is_List["is_z_sd_fixed"];
	bool is_z_constant_fixed = is_List["is_z_constant_fixed"];

	// Extract fixed values
	List fixed_List = model_Lists["fixed_List"];
	double z_constant_fixed = fixed_List["z_constant_fixed"];

	// Other stuff
	NumericVector x1 = model["x1"];

	DataFrame data = as_data_frame(model["dataframe"]);
	StringVector data_names = data.names();

	NumericVector p_values = results(_, 3);
	StringVector stars = starVector(p_values);

	NumericVector z_coef = model["coefficients"];

	StringVector results_rownames = rownames(results);
	StringVector results_colnames = colnames(results);

	double mean = model["mean"];
	double sd = model["sd"];

	double lnL = model["log-likelihood"];
	double AIC = model["AIC"];
	int n_obs = model["n_obs"];
	double n_obs_double = n_obs * 1.0;
	double BIC = log(n_obs_double) * x1.size() - 2 * lnL;
	int df = x1.size();
	std::string lnL_string = "Log-Likelihood: " + std::to_string(lnL) + "\n";
	std::string AIC_string = "AIC: " + std::to_string(AIC) + "\n";
	std::string BIC_string = "BIC: " + std::to_string(BIC) + "\n";
	std::string n_obs_string = "Observations: " + std::to_string(n_obs) + "\n";
	std::string df_string = std::to_string(df) + 
	                        " free parameters (df = " + 
	                        std::to_string(n_obs - df) + ")" + "\n";

	cat_R("--------------------------------------------------------------\n");
	cat_R("Semi-nonparametric binary choice model estimation\n");
	cat_R(lnL_string.c_str());
	cat_R(AIC_string.c_str());
	cat_R(n_obs_string.c_str());
	cat_R(df_string.c_str());
	cat_R("---\n");
	cat_R("Coefficients:\n");
	
	int z_coef_first = z_coef_ind[0];
	int z_coef_last = z_coef_ind[z_coef_ind.size() - 1];
	NumericMatrix z_coef_results = results(Range(z_coef_first, z_coef_last), _);
	rownames(z_coef_results) = results_rownames[z_coef_ind];
	colnames(z_coef_results) = results_colnames;
	print_R(as_table(cbind(z_coef_results, stars[z_coef_ind])));

	cat_R("---\n");
	cat_R("Distribution parameters:\n");
	
	int distr_first = 0;
	int distr_last = df - z_coef_ind.size() - 1;
	NumericMatrix distr_results = results(Range(distr_first, distr_last), _);
	StringVector distr_rownames = results_rownames[Range(distr_first, distr_last)];
	rownames(distr_results) = distr_rownames;
	colnames(distr_results) = results_colnames;
	print(as_table(cbind(distr_results, stars[Range(distr_first, distr_last)])));

	cat_R("---\n");
	cat_R("Fixed Coefficients:\n");
	
	if (is_z_constant_fixed)
	{
		std::string new_str = "(Intercept) = " + 
		                      std::to_string(z_constant_fixed) + "\n";
	  cat_R(new_str.c_str());
	}

	if (is_z_coef_first_fixed)
	{
		String new_str_names = data_names(1);
		std::string new_str_names_str = new_str_names;
		std::string new_str = new_str_names_str + " = 1" + "\n";
		cat_R(new_str.c_str());
	}

	cat_R("---\n");
	cat_R("Fixed Distribution Parameters:\n");
	cat_R("a_0 = 1\n");
	
	if (is_z_mean_fixed)
	{
		std::string new_str = "mean = " + std::to_string(mean) + "\n";
	  cat_R(new_str.c_str());
	}

	if (is_z_sd_fixed)
	{
		std::string new_str = "sd = " + std::to_string(sd) + "\n";
	  cat_R(new_str.c_str());
	}

	cat_R("---\n");
	cat_R("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n");
	cat_R("--------------------------------------------------------------\n");
}

//' Plot hpaBinary random errors approximated density
//' @param x Object of class "hpaBinary"
//' @export	
// [[Rcpp::export]]
void plot_hpaBinary(List x) 
{
  // Extract the model
	List model = x;

	// Load additional environments
	Rcpp::Environment base_env("package:base");
	Rcpp::Function round_R = base_env["round"];
	Rcpp::Function Version_R = base_env["R.Version"];
	Rcpp::Function plot_R = base_env["c"];

	// For compatibility with earlier versions of R
	List Version_R_List = Version_R();
	CharacterVector major_vec = Version_R_List["major"];
	if(major_vec[0] == "4")
	{
	  plot_R = base_env["plot"];
	} else {
	  Rcpp::Environment graphics_env("package:graphics");
	  plot_R = graphics_env["plot"];
	}
	
	// Load data from the model
	NumericVector pol_coefficients = model["pol_coefficients"];
	NumericVector pol_degrees = model["pol_degrees"];

	NumericVector mean = model["mean"];
	NumericVector sd = model["sd"];

	// Adjust precision
	double errors_exp = model["errors_exp"];
	double errors_var = model["errors_var"];

	double plot_min = errors_exp - 3 * sqrt(errors_var);
	double plot_max = errors_exp + 3 * sqrt(errors_var);

	int n = 10000;

	double precise = (plot_max - plot_min) / n;

	NumericMatrix x_matr = NumericMatrix(n, 1);
	x_matr(0, 0) = plot_min;
	
	for (int i = 1; i < n; i++)
	{
		x_matr(i, 0) = x_matr(i - 1, 0) + precise;
	}

	NumericVector x_vec = x_matr(_, 0);

	// Calculate densities
	NumericVector den = dhpa(x_matr,
		pol_coefficients, pol_degrees,
		LogicalVector::create(false),
		LogicalVector::create(false),
		mean, sd, false, false, false);

	double den_min = min(den);
	double den_max = max(den);

	// Build the plot
	plot_R(Rcpp::_["x"] = x_vec, Rcpp::_["y"] = den,
		Rcpp::_["xlim"] = NumericVector::create(plot_min, plot_max),
		Rcpp::_["xaxp"] = NumericVector::create(plot_min, plot_max, 5),
		Rcpp::_["yaxp"] = NumericVector::create(den_min, den_max, 5),
		Rcpp::_["type"] = "l",
		Rcpp::_["lwd"] = 3,
		Rcpp::_["main"] = "Random Errors Density Approximation Plot",
		Rcpp::_["xlab"] = "value", 
		Rcpp::_["ylab"] = "density");
}

//' Calculates log-likelihood for "hpaBinary" object
//' @description This function calculates log-likelihood for "hpaBinary" object
//' @param object Object of class "hpaBinary"
//' @export	
// [[Rcpp::export]]
double logLik_hpaBinary(List object) 
{
	double lnL = object["log-likelihood"];

	return(lnL);
}
