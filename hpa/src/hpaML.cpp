#include "hpaMain.h"
#include "hpaML.h"
#include "polynomialIndex.h"
#include "hpaValidation.h"
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace RcppParallel;

// [[Rcpp::depends(RcppArmadillo)]]

//' Semi-nonparametric maximum likelihood estimation
//' @description This function performs semi-nonparametric (SNP)
//' maximum likelihood estimation of unknown (possibly truncated) multivariate  
//' density using Hermite polynomial based approximating function proposed by 
//' Gallant and Nychka in 1987. Please, see \code{\link[hpa]{dhpa}} 'Details' 
//' section to get more information concerning this approximating function.
//' @template x_ML_Template
//' @template pol_degrees_Template
//' @template tr_left_vec_Template
//' @template tr_right_vec_Template
//' @template given_ind_Template
//' @template omit_ind_Template
//' @template x0_ML_Template
//' @template cov_type_Template
//' @template boot_iter_Template
//' @template is_parallel_Template
//' @template opt_type_Template
//' @template opt_control_Template
//' @template is_validation_Template
//' @template GN_details_Template
//' @template hpaML_formula_Template
//' @template parametric_paradigm_Template
//' @template optim_details_Template
//' @template opt_control_details_Template
//' @template opt_control_details_hpaML_Template
//' @return This function returns an object of class "hpaML".\cr \cr
//' An object of class "hpaML" is a list containing the following components:
//' \itemize{
//' \item \code{optim} - \code{\link[stats]{optim}} function output. 
//' If \code{opt_type = "GA"} then it is the list containing 
//' \code{\link[stats]{optim}} and \code{\link[GA]{ga}} functions outputs.
//' \item \code{x1} - numeric vector of distribution parameters estimates.
//' \item \code{mean} - density function mean vector estimate.
//' \item \code{sd} - density function sd vector estimate.
//' \item \code{pol_coefficients} - polynomial coefficients estimates.
//' \item \code{tr_left }- the same as \code{tr_left} input parameter.
//' \item \code{tr_right} - the same as \code{tr_right} input parameter.
//' \item \code{omit_ind }- the same as \code{omit_ind} input parameter.
//' \item \code{given_ind} - the same as \code{given_ind} input parameter.
//' \item \code{cov_mat} - covariance matrix estimate.
//' \item \code{results} - numeric matrix representing estimation results.
//' \item \code{log-likelihood} - value of Log-Likelihood function.
//' \item \code{AIC} - AIC value.
//' \item \code{data} - the same as \code{x} input parameter but without \code{NA} observations.
//' \item \code{n_obs} - number of observations.
//' \item \code{bootstrap} - list where bootstrap estimation results are stored.}
//' @seealso \link[hpa]{summary.hpaML}, \link[hpa]{predict.hpaML}, 
//' \link[hpa]{logLik.hpaML}, \link[hpa]{plot.hpaML}
//' @template hpaML_examples_Template
//' @export
// [[Rcpp::export]]
List hpaML(NumericMatrix x,
	NumericVector pol_degrees = NumericVector(0),
	NumericVector tr_left = NumericVector(0),
	NumericVector tr_right = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector x0 = NumericVector(0),
	String cov_type = "sandwich",
	int boot_iter = 100,
	bool is_parallel = false,
	String opt_type = "optim",
	List opt_control = R_NilValue,
	bool is_validation = true)
{
  // Validation Stuff
  
  if(is_validation)
  {
      // Validate covariance matrix type
    if((cov_type != "sandwich") & (cov_type != "sandwichFD") &
       (cov_type != "bootstrap") & (cov_type != "gop") & 
       (cov_type != "hessian") & (cov_type != "hessianFD"))
    {
      stop("Incorrect cov_type argument value.");
    }
    
      // Validate opt_type
    if((opt_type != "optim") & (opt_type != "GA"))
    {
      stop("Incorrect opt_type argument value.");
    }
    
      // Warning concerning large number of bootstrap iterations
    if(boot_iter > 1000)
    {
      warning("Since boot_iter is large estimation may take lots of time.");
    }
  
    int target_dim = x.ncol();
    
    // Validate polynomial degrees
    pol_Validate(pol_degrees, NumericVector(0));

    if (pol_degrees.size() != target_dim)
    {
      stop("pol_degrees length should be the same as the number of x columns.");
    }
    
    // Validate conditional and omitted values
    ind_Validate(given_ind, omit_ind);
    
    int n_given_ind = given_ind.size();
    if((n_given_ind != 0) & (n_given_ind != target_dim))
    {
      stop("given_ind length should be the same as the length of pol_degrees");
    }
    
    int n_omit_ind = omit_ind.size();
    if((n_omit_ind != 0) & (n_omit_ind != target_dim))
    {
      stop("omit_ind length should be the same as the length of pol_degrees");
    }
    
    // Validate truncation points
    if((tr_left.size() != tr_right.size()) & 
       (tr_left.size() != 0) & (tr_right.size() != 0))
    {
      stop("tr_left and tr_right should be vectors of the same dimensions");
    }
  }
	// Load additional environments

		// stats environment
	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function na_omit_R = stats_env["na.omit"];
	Rcpp::Function optim = stats_env["optim"];
	Rcpp::Function cov_R = stats_env["cov"];

		// base environment
	Rcpp::Environment base_env("package:base");
	Rcpp::Function c_R = base_env["c"];
	Rcpp::Function diag_R = base_env["diag"];
	Rcpp::Function requireNamespace_R = base_env["requireNamespace"];
	Rcpp::Function cat_R = base_env["cat"];
	
	  // GA environment
	Rcpp::Function ga_R = stats_env["optim"];
	Rcpp::Function ga_summary_R = stats_env["optim"];
	
	    // should be LogicalVector not bool otherwise not working
	LogicalVector is_ga_installed = requireNamespace_R(Rcpp::_["package"] = "GA", 
                                                     Rcpp::_["quietly"] = true);

	    // check weather GA package has been installed
	if((opt_type == "GA") | (opt_type == "ga"))
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
	
	// Remove NA values from data
	x = clone(x);
	x = na_omit_R(x);
	
	// Get the number of observations
	int n_obs = x.nrow();
	
	// Initialize polynomial structure 
	// related values
	int pol_degrees_n = pol_degrees.size();  // random vector dimensionality
	int pol_coefficients_n = 1;              // number of polynomial coefficients

	// Calculate the number of 
	// polynomial coefficients
	for (int i = 0; i < pol_degrees_n; i++)
	{
		pol_coefficients_n *= (pol_degrees[i] + 1); // +1 because starts from 0
	}
	pol_coefficients_n -= 1;                      // because of a(0...0) = 1

	// Initialize conditions and marginals
	if (given_ind.size() == 0)                    // if there is no 
	{                                             //conditioned components
		given_ind = LogicalVector(pol_degrees_n);   // false by default
	}

	if (omit_ind.size() == 0)                     // if there is no 
	{                                             // marginalized components
		omit_ind = LogicalVector(pol_degrees_n);    // false by default
	}

	// Determine whether initial values have been manually provided
	bool x0_given = true; // initially assume that they have been

	if (x0.size() == 0) // if x0 has not been provided manually
	{
		x0_given = false; // then initial values have not been provided
	  
		x0 = NumericVector(pol_coefficients_n + // 2 * pol_degrees_n since every 
		                                        // random vector components
		                   2 * pol_degrees_n);  // introduces additional mean and 
	}                                         // sd parameters pair
	
	// Initialize additional variable which helps
	// to assign indices of parameters in x0
	int k = 0;

	// Assign indices for

		// polynomial coefficients
	NumericVector pol_coefficients_ind(pol_coefficients_n);

	for (int i = 0; i < pol_coefficients_n; i++)
	{
	  if (!x0_given) // if user has not provided x0 manually
	  {              // then set mean parameter to sample mean
	    x0[i] = 0;
	  }
	  
		pol_coefficients_ind[i] = i;
	}

		// mean vector
	NumericVector mean_ind(pol_degrees_n);

	for (int i = pol_coefficients_n; i < (pol_coefficients_n + pol_degrees_n); i++)
	{
		mean_ind[k] = i;
	  
		if (!x0_given) // if user has not provided x0 manually then 
		               // set mean parameter to sample mean
		{
			x0[i] = mean(x(_, k));
		}
		
		k++;
	}
	
		// sd vector
	k = 0;

	NumericVector sd_ind(pol_degrees_n);

	for (int i = (pol_coefficients_n + pol_degrees_n); 
       i < (pol_coefficients_n + 2 * pol_degrees_n); i++)
	{
		sd_ind[k] = i;
	  
		if (!x0_given) // if user has not provided x0 manually then set
		               // sd parameter to sample standard deviation
		{
			x0[i] = sd(x(_, k));
		}
		
		k++;
	}

	// Deal with truncation
	NumericMatrix tr_left_mat(1, pol_degrees_n);
	NumericMatrix tr_right_mat(1, pol_degrees_n);

	if ((tr_left.size() > 0) | (tr_right.size() > 0))
	{
		if (tr_left.size() == 0) // if there is no left truncation
		{                        // set it to negative infinity
			tr_left = NumericVector(pol_degrees_n, R_NegInf);
		}

		if (tr_right.size() == 0) // if there is no right truncation
		{                         // set it to infinity
			tr_right = NumericVector(pol_degrees_n, R_PosInf);
		}

		tr_left_mat(0, _) = tr_left;
		tr_right_mat(0, _) = tr_right;
	} else {
		std::fill(tr_left_mat.begin(), tr_left_mat.end(), NA_REAL);
		std::fill(tr_right_mat.begin(), tr_right_mat.end(), NA_REAL);
	}
	
	// Apply optimization routine
	
	  // Set optim control parameters
	List PGN_control = List::create(
	     Named("maxit") = 100000000, 
       Named("fnscale") = -1.0,
       Named("abstol") = std::sqrt(
         std::numeric_limits<double>::epsilon()) * 0.01,
       Named("reltol") = std::sqrt(
         std::numeric_limits<double>::epsilon()) * 0.01);
	
	List hpaML_args = List::create(Named("x_data") = x,
                        Named("pol_coefficients_ind") = pol_coefficients_ind,
                        Named("pol_degrees") = pol_degrees,
                        Named("given_ind") = given_ind,
                        Named("omit_ind") = omit_ind,
                        Named("mean_ind") = mean_ind,
                        Named("sd_ind") = sd_ind,
                        Named("tr_left") = tr_left_mat,
                        Named("tr_right") = tr_right_mat,
                        Named("is_parallel") = is_parallel);
	
	  // Perform the optimization
	List optim_results = optim(
	    Rcpp::_["par"] = x0,
	    Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaLnLOptim),
	    Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaLnLOptim_grad),
	    Rcpp::_["control"] = PGN_control,
	    Rcpp::_["method"] = "BFGS",
	    Rcpp::_["hessian"] = true,
	    Rcpp::_["hpaML_args"] = hpaML_args);

	// Extract the point of maximum from optim function
	NumericVector x1 = optim_results["par"];
	
	int x1_n = x1.size();
	
	  // Genetic algorithm
	List ga_List;
	List ga_summary;
	  
	if(opt_type == "GA")
	{
  	  //suggest initial solution based on local optimization
  	NumericMatrix ga_suggestions;
	  if(opt_control.containsElementNamed("suggestions"))
	  {
	    ga_suggestions = Rcpp::as<Rcpp::NumericMatrix>(
	      opt_control["suggestions"]);
	  } else {
	    ga_suggestions= NumericMatrix(1, x1_n);
  	  ga_suggestions(0,_) = x1;
	  }
  	
  	  // set lower and upper bounds for parameters space
  	NumericVector ga_lower = NumericVector(x1_n);
  	NumericVector ga_upper = NumericVector(x1_n);
  	
  	if(opt_control.containsElementNamed("lower") & 
       opt_control.containsElementNamed("upper"))
  	{
  	  ga_lower = opt_control["lower"];
  	  ga_upper = opt_control["upper"];
  	} else {
  	    // bounds for the mean parameter
  	  NumericVector x1_mean = x1[mean_ind];
  	  
  	  NumericVector ga_lower_mean = x1_mean - 2 * abs(x1_mean);
  	  NumericVector ga_upper_mean = x1_mean + 2 * abs(x1_mean);
  	  
  	  ga_lower[mean_ind] = ga_lower_mean;
  	  ga_upper[mean_ind] = ga_upper_mean;
  	  
  	    // bounds for the sd parameter
  	  NumericVector ga_lower_sd = abs(x1[sd_ind]) * 0.2;
  	  NumericVector ga_upper_sd = abs(x1[sd_ind]) * 5;
  	  
  	  ga_lower[sd_ind] = ga_lower_sd;
  	  ga_upper[sd_ind] = ga_upper_sd;
  	  
  	    // bounds for thepolynomial coefficients parameters
  	  ga_lower[pol_coefficients_ind] = -10;
  	  ga_upper[pol_coefficients_ind] = 10;
  	}
  	
      // set maximum number of iterations
  	int ga_maxiter;
  	if(opt_control.containsElementNamed("maxiter"))
  	{
  	  ga_maxiter = opt_control["maxiter"];
  	} else {
  	  ga_maxiter = 50 * (1 + pol_coefficients_n);
  	}

  	  // set population size
  	int ga_popSize;
  	if(opt_control.containsElementNamed("popSize"))
  	{
  	  ga_popSize = opt_control["popSize"];
  	} else {
  	  ga_popSize = 10 + pol_coefficients_n * 2;
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
  	  ga_elitism = 2 + (int)(std::round(0.1 * ga_popSize));
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
  	  Rcpp::_["fitness"] = Rcpp::InternalFunction(&hpaLnLOptim),
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
      Rcpp::_["hpaML_args"] = hpaML_args,
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
  	  Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaLnLOptim),
  	  Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaLnLOptim_grad),
  	  Rcpp::_["control"] = PGN_control,
  	  Rcpp::_["method"] = "BFGS",
  	  Rcpp::_["hessian"] = true,
  	  Rcpp::_["hpaML_args"] = hpaML_args);

  	// Reset optimal point
  	x1 = optim_results["par"];
	}
	
	// Extract optimization results and assign them to the variables
	// representing estimated parameters
	double lnL = optim_results["value"];

	NumericVector mean = x1[mean_ind];
	NumericVector sd = x1[sd_ind];

	NumericVector pol_coefficients = x1[pol_coefficients_ind];
	pol_coefficients.push_front(1);

	// Get covariance matrix estimate of "cov_type" type
	NumericMatrix cov_mat;                                // covariance matrix
	
	arma::mat H_part;                                     // butter of sandwich 
	arma::mat J_part;                                     // bread of sandwich
	
	NumericMatrix my_hessian = optim_results["hessian"];  // Hessian matrix
	
	  // Estimate hessian
	if ((cov_type == "hessianFD") | (cov_type == "sandwichFD"))
	{
	  try
	  {
  	  my_hessian = hpaLnLOptim_hessian(x1, hpaML_args);
	  } catch (std::exception &ex) {
	    warning("Can't calculate Hessian via first difference method. Hessian from the optim function will be used instead.");
	    forward_exception_to_r(ex);
	  }
	}
	
	  // estimate inverse Hessian
	if ((cov_type == "hessian") | (cov_type == "sandwich") |
      (cov_type == "hessianFD") | (cov_type == "sandwichFD"))
	{
	  
	  H_part = as<arma::mat>(my_hessian);
	  
	  int H_rank = rank(H_part);
	  if(H_rank == x1_n)
	  {
	    H_part = (as<arma::mat>(my_hessian)).i();
	  } else {
	    cov_type = "gop";
	    warning("Can't get inverse Hessian. GOP covariance matrix estimate will be returned.");
	  }
	}
	
	  // Estimate jacobian
	if ((cov_type == "gop") | (cov_type == "sandwich") |
      (cov_type == "sandwichFD"))
	{
	  NumericMatrix my_jacobian = hpaLnLOptim_grad_ind(x1, hpaML_args);
	  
	  J_part = as<arma::mat>(my_jacobian);
	}

	  // Sandwich estimate
	if ((cov_type == "sandwich") | (cov_type == "sandwichFD"))
	{
	  cov_mat = wrap(H_part * (J_part.t() * J_part) * H_part);
	}
	
	  // Inverse hessian estimate
	if ((cov_type == "hessian") | (cov_type == "hessianFD"))
	{
	  cov_mat = wrap(-H_part);
	}
	
	  // Gradient outer product estimate (should be the last)
	if (cov_type == "gop")
	{
	  arma::mat GOP_mat = J_part.t() * J_part;
	  
	  int H_rank = rank(GOP_mat);
	  if(H_rank == x1_n)
	  {
	    cov_mat = wrap(GOP_mat.i());
	  } else {
	    cov_mat = NumericMatrix(x1_n, x1_n);
	    cov_mat.fill_diag(1.0);
	    warning("Can't calculate covariance matrix. Identity matrix will be used instead.");
	  }
	}
	
	  // Apply bootstrap
	  
	    // store parameters for each iteration
	  NumericMatrix boot_parameters = NumericMatrix(boot_iter, x1_n);
	
	    // store standard deviation
	  NumericVector sd_dev = NumericVector(x1_n);
	  
	    // list to store bootstrap results
	  List boot_List;
	  
	    // bootstrap procedure
	 if (cov_type == "bootstrap")
	 {
	   for(int i = 0; i < boot_iter; i++)
	   {
	     // Generate sample with replacement
	     NumericVector sample_ind = floor(runif(n_obs, 0, n_obs));
	     NumericMatrix boot_sample = NumericMatrix(n_obs, pol_degrees_n);
	     
	     for (int j = 0; j < n_obs; j++)
	     {
	       boot_sample(j, _) = x(sample_ind[j], _);
	     }
	     
	     // Prepare arguments for bootstrap List
	     List hpaML_args_boot = hpaML_args;
	     hpaML_args_boot["x_data"] = boot_sample;

	     // Perform estimaton
	     List boot_results = optim(
	       Rcpp::_["par"] = x1,
	       Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaLnLOptim),
	       Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaLnLOptim_grad),
	       Rcpp::_["control"] = PGN_control,
	       Rcpp::_["method"] = "BFGS",
	       Rcpp::_["hessian"] = false,
	       Rcpp::_["hpaML_args"] = hpaML_args_boot);
	     
	     // Store iteration results
	     NumericVector x1_new = boot_results["par"];
	     boot_parameters(i, _) = x1_new;
	   }
	   
	   // Store bootstrap results
	   cov_mat = cov_R(boot_parameters);
	   sd_dev = sqrt(diag_R(cov_mat));
	   
	   // Store the results into the list
	   boot_List = List::create(
	     Named("estimates") = boot_parameters,
	     Named("cov_mat") = cov_mat,
	     Named("sd") = sd_dev);
	 }
	
	// Prepare beautiful results output
	
	  // Matrix to store the results
	NumericMatrix results(x1_n, 3);

	  // Provide columns names for the matrix
	StringVector results_cols = StringVector::create(
	  "Estimate", "Std. Error", "P(>|z|)");
	StringVector results_rows(x1_n);

		// Get vector index matrix for polynomial coefficients
	NumericMatrix pol_ind = polynomialIndex(pol_degrees, false);

		// Assign results matrix columns with values for

			// polynomial coefficients
	for (int i = 1; i < (pol_coefficients_n + 1); i++)
	{
		results_rows[(i - 1)] = "a";
		for (int j = 0; j < pol_degrees_n; j++)
		{
			int my_int = pol_ind(j, i);
			results_rows[(i - 1)] = as<std::string>(results_rows[(i - 1)]) + 
			                                        "_" + std::to_string(my_int);
		}
		results((i - 1), 0) = pol_coefficients[i];
		results((i - 1), 1) = sqrt(cov_mat((i - 1), (i - 1)));
		double t_stat = results((i - 1), 0) / results((i - 1), 1);
		NumericVector F_t_stat = pnorm(NumericVector::create(t_stat));
		results((i - 1), 2) = 2 * std::min(F_t_stat[0], 1 - F_t_stat[0]);
	}

			// mean
	for (int i = 0; i < pol_degrees_n; i++)
	{
		if (pol_degrees_n == 1)
		{
			results_rows[mean_ind[i]] = "mean";
		} else {
			results_rows[mean_ind[i]] = "mean_" + std::to_string(i + 1);
		}
		results(mean_ind[i], 0) = x1[mean_ind[i]];
		results(mean_ind[i], 1) = sqrt(cov_mat((mean_ind[i]), (mean_ind[i])));
		double z_stat = results(mean_ind[i], 0) / results(mean_ind[i], 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(mean_ind[i], 2) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

			// sd
	for (int i = 0; i < pol_degrees_n; i++)
	{
		if (pol_degrees_n == 1)
		{
			results_rows[sd_ind[i]] = "sd";
		} else {
			results_rows[sd_ind[i]] = "sd_" + std::to_string(i + 1);
		}
		results(sd_ind[i], 0) = x1[sd_ind[i]];
		results(sd_ind[i], 1) = sqrt(cov_mat((sd_ind[i]), (sd_ind[i])));
		double z_stat = results(sd_ind[i], 0) / results(sd_ind[i], 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(sd_ind[i], 2) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

		// assign names to rows and some vectors
	rownames(results) = results_rows;
	colnames(results) = results_cols;

	x1.names() = results_rows;

	rownames(cov_mat) = results_rows;
	colnames(cov_mat) = results_rows;

	pol_coefficients.names() = c_R("a_0", results_rows[pol_coefficients_ind]);
	
	if (cov_type == "bootstrap")
	{
  	colnames(boot_parameters) = results_rows;
  	
  	sd_dev.names() = results_rows;
	}

	// Collect the results
	
	  // Store the results related to optimization routines
	if((opt_type == "GA"))
	{
	  optim_results = List::create(
	    Named("optim") = optim_results,
	    Named("GA") = ga_List,
	    Named("GA_summary") = ga_summary);
	}

	  // Make the list containing the objects to return
	List return_result = List::create(
		Named("optim") = optim_results, 
		Named("x1") = x1,
		Named("mean") = mean, 
		Named("sd") = sd,
		Named("pol_coefficients") = pol_coefficients, 
		Named("pol_degrees") = pol_degrees,
		Named("tr_left") = tr_left_mat, 
		Named("tr_right") = tr_right_mat,
		Named("omit_ind") = omit_ind, 
		Named("given_ind") = given_ind,
		Named("cov_matrix") = cov_mat,
		Named("results") = results, 
		Named("log-likelihood") = lnL, 
		Named("AIC") = 2 * (x1_n - lnL),
		Named("data") = x,
		Named("n_obs") = n_obs,
		Named("bootstrap") = boot_List);

	// Assign the class to the output list
	return_result.attr("class") = "hpaML";

	return(return_result);
}

// Perform log-likelihood function estimation for 
// Phillips-Gallant-Nychka distribution at point
List hpaLnLOptim_List(NumericVector x0, List hpaML_args)
{
  // Get arguments from the hpaML_args
  NumericMatrix x_data = hpaML_args["x_data"];
  NumericVector pol_coefficients_ind = hpaML_args["pol_coefficients_ind"];
  NumericVector pol_degrees = hpaML_args["pol_degrees"];
  LogicalVector given_ind = hpaML_args["given_ind"];
  LogicalVector omit_ind = hpaML_args["omit_ind"];
  NumericVector mean_ind = hpaML_args["mean_ind"];
  NumericVector sd_ind = hpaML_args["sd_ind"];
  NumericMatrix tr_left = hpaML_args["tr_left"];
  NumericMatrix tr_right = hpaML_args["tr_right"];
  bool is_parallel = hpaML_args["is_parallel"];
  
  // Get number of observations
  int n = x_data.nrow();
  
  // Initialize values to return 
  List return_List;
  
  NumericVector return_individual = NumericVector(n);
  double return_aggregate;
  
	// Assign values based on their indices
	NumericVector mean = x0[mean_ind];
	NumericVector sd = x0[sd_ind];
	
	// Control for positive standard deviations
	if(sum(sd <= 0) > 0)
	{
	  return_aggregate = R_NegInf;
	  std::fill(return_individual.begin(), return_individual.end(), R_NegInf);
	  
	  return_List = List::create(Named("aggregate") = R_NegInf,
                               Named("individual") = return_individual);
	  return(return_List);
	}

	NumericVector pol_coefficients = x0[pol_coefficients_ind];
	pol_coefficients.push_front(1);

	// Perform calculations
	if (!(R_IsNA(tr_left(0, 0))) & !(R_IsNA(tr_right(0, 0))))
	{
	  return_individual = dtrhpa(x_data,
                               tr_left, tr_right,
                               pol_coefficients, pol_degrees,
                               given_ind, omit_ind,
                               mean, sd,
                               is_parallel, true, false);
	  
	  return_aggregate = sum(return_individual);
	  
	  // if there is some problems related to precision when truncation
	  // has been incorporated then return great negative number
	  
	  if(ihpa(tr_left, tr_right,
	    pol_coefficients, pol_degrees,
	    given_ind, omit_ind, 
	    mean, sd,
	    is_parallel, false, false)[0] < std::sqrt(
	        std::numeric_limits<double>::epsilon()))
	  {
	    std::fill(return_individual.begin(), return_individual.end(), R_NegInf);
	    
	    return_List = List::create(Named("aggregate") = R_NegInf,
                                 Named("individual") = return_individual);
	    
	    return(return_List);
	  }
	  
	  return_List = List::create(Named("aggregate") = return_aggregate,
                               Named("individual") = return_individual);

		return(return_List);
	}
	
	// Calculate log-likelihood values for each observations
	return_individual = dhpa(x_data,
                  		     pol_coefficients, pol_degrees,
                  		     given_ind, omit_ind,
                  		     mean, sd, 
                  		     is_parallel, true, false);
	
	// Calculate log-likelihood function value
	return_aggregate = sum(return_individual);
	
	// Aggregate and return the results
	return_List = List::create(Named("aggregate") = return_aggregate,
                             Named("individual") = return_individual);

	return(return_List);
}

// Get aggregate component (log-likelihood function value) from hpaLnLOptim_List
double hpaLnLOptim(NumericVector x0, List hpaML_args)
{
  List return_List = hpaLnLOptim_List(x0, hpaML_args);
  
  double return_aggregate = return_List["aggregate"];
  
  if(any(is_na(NumericVector::create(return_aggregate))))
  {
    return_aggregate = R_NegInf;
  }
  
  return(return_aggregate);
}

// Get individual component (log-lkelihood function contributions) 
// from hpaLnLOptim_List
NumericVector hpaLnLOptim_ind(NumericVector x0, List hpaML_args)
{
  List return_List = hpaLnLOptim_List(x0, hpaML_args);
  
  NumericVector return_individual = return_List["individual"];
  
  return(return_individual);
}

// Perform log-likelihood function gradient estimation 
// for Phillips-Gallant-Nychka distribution at point
List hpaLnLOptim_grad_List(NumericVector x0, List hpaML_args)
{
  // Get arguments from the hpaML_args
  NumericMatrix x_data = hpaML_args["x_data"];
  NumericVector pol_coefficients_ind = hpaML_args["pol_coefficients_ind"];
  NumericVector pol_degrees = hpaML_args["pol_degrees"];
  LogicalVector given_ind = hpaML_args["given_ind"];
  LogicalVector omit_ind = hpaML_args["omit_ind"];
  NumericVector mean_ind = hpaML_args["mean_ind"];
  NumericVector sd_ind = hpaML_args["sd_ind"];
  NumericMatrix tr_left = hpaML_args["tr_left"];
  NumericMatrix tr_right = hpaML_args["tr_right"];
  bool is_parallel = hpaML_args["is_parallel"];
  
  // Get parameters number
  int n_param = x0.size();
  int n_obs = x_data.nrow();
  
  // Initialize values to return 
  List return_List;
  
  NumericVector return_individual = NumericVector(n_obs);
  double return_aggregate;
  
  // Initialize vector to store gradient (jacobian) values
  NumericMatrix my_grad = NumericMatrix(n_obs, n_param);
  
  // Assign values based on their indecies
  NumericVector mean = x0[mean_ind];
  NumericVector sd = x0[sd_ind];
  
  NumericVector pol_coefficients = x0[pol_coefficients_ind];
  pol_coefficients.push_front(1);
  
  // Control for non-negative standard deviations
  if(sum(sd < 0) > 0)
  {
    return_aggregate = R_NegInf;
    
    std::fill(return_individual.begin(), 
              return_individual.end(), 
              R_NegInf);
    
    return_List = List::create(Named("aggregate") = return_aggregate,
                               Named("individual") = return_individual);
    
    return(return_List);
  }
  
  // Gradient estimation
  
   // while there is excess calculation of gradient for
   // x it is still much faster then calculating gradients
   // for pol_coefficients, mean and sd separately since their
   // gradient share many common parts while gradient of x
   // is fast to calculate given this common parts
  NumericMatrix all_grad = dhpaDiff(x_data,
                                    pol_coefficients, pol_degrees,
                                    given_ind, omit_ind,
                                    mean, sd,
                                    "all",
                                    is_parallel, true, false);
  
  for(int i = 0; i < n_param; i++)
  {
    my_grad(_, i) = all_grad(_, i + 1);
  }
  
    // Deal with truncation
  if (!(R_IsNA(tr_left(0, 0))) & !(R_IsNA(tr_right(0, 0))))
  {
    NumericMatrix tr_grad = ihpaDiff(tr_left, tr_right,
                                     pol_coefficients, pol_degrees,
                                     given_ind, omit_ind,
                                     mean, sd,
                                     "all",
                                     is_parallel, true, false);

    for (int i = 0; i < n_param; i++)
    {
      my_grad(_, i) = my_grad(_, i) - sum(tr_grad(_, i + 1));
    }
  }
  // Return the results
  
  return_List = List::create(Named("aggregate") = colSums(my_grad),
                             Named("individual") = my_grad);

  return(return_List);
}

// Get aggregaate component (log-lkelihood function gradient) 
// from hpaLnLOptim_grad_List
NumericVector hpaLnLOptim_grad(NumericVector x0, List hpaML_args)
{
  List return_List = hpaLnLOptim_grad_List(x0, hpaML_args);
  
  NumericVector return_aggregate = return_List["aggregate"];
  
  if(any(is_na(return_aggregate)))
  {
    std::fill(return_aggregate.begin(), 
              return_aggregate.end(), 
              R_NegInf);
  }
  
  return(return_aggregate);
}

// Get individual component (log-lkelihood function gradient contributions) 
// from hpaLnLOptim_grad_List
NumericMatrix hpaLnLOptim_grad_ind(NumericVector x0, List hpaML_args)
{
  List return_List = hpaLnLOptim_grad_List(x0, hpaML_args);
  
  NumericMatrix return_individual = return_List["individual"];
  
  return(return_individual);
}

// Perform log-likelihood function hessian estimation 
// for Phillips-Gallant-Nychka distribution at point
NumericMatrix hpaLnLOptim_hessian(NumericVector x0, List hpaML_args)
{
  // Get arguments from the hpaML_args
  NumericMatrix x_data = hpaML_args["x_data"];
  NumericVector pol_coefficients_ind = hpaML_args["pol_coefficients_ind"];
  NumericVector pol_degrees = hpaML_args["pol_degrees"];
  LogicalVector given_ind = hpaML_args["given_ind"];
  LogicalVector omit_ind = hpaML_args["omit_ind"];
  NumericVector mean_ind = hpaML_args["mean_ind"];
  NumericVector sd_ind = hpaML_args["sd_ind"];
  NumericMatrix tr_left = hpaML_args["tr_left"];
  NumericMatrix tr_right = hpaML_args["tr_right"];
  
  // Get parameters number
  int n_param = x0.size();
  
  // Initialize vector to store hessian values
  NumericMatrix my_hessian = NumericMatrix(n_param, n_param);
  
  // Assign values based on their indices
  NumericVector mean = x0[mean_ind];
  NumericVector sd = x0[sd_ind];
  
  NumericVector pol_coefficients = x0[pol_coefficients_ind];
  pol_coefficients.push_front(1);
  
  // Prepare precision related values for numeric differentiation
  double machinePrecision = std::numeric_limits<double>::epsilon();
  double my_precision = std::sqrt(machinePrecision);
  
  NumericVector eps = abs(x0 * my_precision);
  
  // Control for zero values
  eps[eps < (machinePrecision * 100)] = machinePrecision * 100;
  
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
    g_plus = hpaLnLOptim_grad(x0_eps, hpaML_args);

    // Calculate g(x - eps)
    x0_eps[i] = x0[i] - eps[i];
    g_minus = hpaLnLOptim_grad(x0_eps, hpaML_args);
    
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

//' Predict method for hpaML
//' @param object Object of class "hpaML"
//' @template newdata_Template
//' @return This function returns predictions based 
//' on \code{\link[hpa]{hpaML}} estimation results.
//' @export
// [[Rcpp::export]]
NumericVector predict_hpaML(List object, 
                            NumericMatrix newdata = NumericMatrix(1, 1))
{
	List model = object;

	// Load additional environments

		// stats environment
	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function na_omit_R = stats_env["na.omit"];

	// Get distribution parameters
	NumericVector pol_coefficients = model["pol_coefficients"];
	NumericVector pol_degrees = model["pol_degrees"];

	NumericVector mean = model["mean"];
	NumericVector sd = model["sd"];

	NumericMatrix tr_left = model["tr_left"];
	NumericMatrix tr_right = model["tr_right"];

	LogicalVector omit_ind = model["omit_ind"];
	LogicalVector given_ind = model["given_ind"];

	NumericMatrix x = model["data"];

	// Get data
	if ((newdata.ncol() == 1) & (newdata.nrow() == 1))
	{
		newdata = x;
	} else {
		newdata = na_omit_R(newdata);
	}

	// Estimate
	if (!(R_IsNA(tr_left(0, 0))) & !(R_IsNA(tr_right(0, 0))))
	{
		return(dtrhpa(newdata,
			tr_left, tr_right,
			pol_coefficients, pol_degrees,
			given_ind, omit_ind,
			mean, sd, false, false, false));
	}

	return(dhpa(newdata,
		pol_coefficients, pol_degrees,
		given_ind, omit_ind,
		mean, sd, false, false, false));
}

//' Summarizing hpaML Fits
//' @param object Object of class "hpaML"
//' @return This function returns the same 
//' list as \code{\link[hpa]{hpaML}} function changing 
//' its class to "summary.hpaML".
//' @export
// [[Rcpp::export]]
List summary_hpaML(List object)
{
	List return_result = clone(object); // in order to preserve model class

	return_result.attr("class") = "summary.hpaML";

	return(return_result);
}

//' Summary for hpaML output
//' @param x Object of class "hpaML"
//' @export
// [[Rcpp::export]]
void print_summary_hpaML(List x)
{
	List model = x;

	// Load additional environments
	Rcpp::Environment base_env("package:base");
	Rcpp::Function as_table = base_env["as.table"];
	Rcpp::Function cbind = base_env["cbind"];
	Rcpp::Function round_R = base_env["round"];
	Rcpp::Function print_R = base_env["print"];
	Rcpp::Function cat_R = base_env["cat"];

	// Do other obvious stuff
	NumericMatrix results = model["results"];
	results = round_R(Rcpp::_["x"] = results, Rcpp::_["digits"] = 5);
	NumericVector x1 = model["x1"];
	
	NumericVector p_values = results(_, 2);
	StringVector stars = starVector(p_values);

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
	  " free parameters (df = " + std::to_string(n_obs - df) + ")" + "\n";

	cat_R("--------------------------------------------------------------\n");
	
	cat_R("Semi-nonparametric maximum likelihood estimation\n");
	
	cat_R("---\n");

	cat_R(lnL_string);
	cat_R(AIC_string);
	cat_R(n_obs_string);
	cat_R(df_string);
	
	cat_R("---\n");

	cat_R("Distribution parameters:\n");
	print_R(as_table(cbind(results, stars)));

	cat_R("---\n");
	cat_R("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n");

	cat_R("--------------------------------------------------------------\n");
}

//' Calculates log-likelihood for "hpaML" object
//' @description This function calculates log-likelihood for "hpaML" object
//' @param object Object of class "hpaML"
//' @export
// [[Rcpp::export]]
double logLik_hpaML(List object)
{
	double lnL = object["log-likelihood"];
  
	return(lnL);
}

// Create star vector
StringVector starVector(NumericVector p_values)
{
  int n = p_values.size();
  
  StringVector stars(n);
  
  for (int i = 0; i < n; i++)
  {
    if (!NumericVector::is_na(p_values[i]))
    {
      if (p_values[i] <= 0.001)
      {
        stars[i] = "***";
      } else {
        if ((0.001 < p_values[i]) & (p_values[i] <= 0.01))
        {
          stars[i] = "**";
        } else {
          if ((0.01 < p_values[i]) & (p_values[i] <= 0.05))
          {
            stars[i] = "*";
          } else {
            if ((0.05 < p_values[i]) & (p_values[i] <= 0.1))
            {
              stars[i] = ".";
            } else {
              stars[i] = " ";
            }
          }
        }
      }
    } else {
      stars[i] = " ";
    }
  }
  
  return(stars);
}

//' Calculates multivariate empirical cumulative distribution function
//' @description This function calculates multivariate 
//' empirical cumulative distribution function
//' at each point of the sample
//' @param x numeric matrix which rows are observations
//' @export
// [[Rcpp::export]]
NumericVector mecdf(NumericMatrix x)
{
  int n = x.nrow();
  int m = x.ncol();
  
  NumericVector my_cdf = NumericVector(n);
  
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      int n_greater = 0;
      
      for (int t = 0; t < m; t++)
      {
        if (x(j, t) <= x(i, t))
        {
          n_greater++;
        }
      }
      
      if (n_greater == m)
      {
        my_cdf[i]++;
      }
      
      if (n_greater == 0)
      {
        my_cdf[j]++;
      }
      
    }
  }
  
  my_cdf = my_cdf / n;
  
  return(my_cdf);
}
