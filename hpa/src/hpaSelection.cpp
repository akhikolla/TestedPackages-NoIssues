#include "hpaMain.h"
#include "hpaML.h"
#include "hpaSelection.h"
#include "hpaBinary.h"
#include "polynomialIndex.h"
#include "hpaValidation.h"

#include <RcppArmadillo.h>

using namespace RcppArmadillo;

// [[Rcpp::depends(RcppArmadillo)]]

//' Perform semi-nonparametric selection model estimation
//' @description This function performs semi-nonparametric (SNP) maximum 
//' likelihood estimation of sample selection model 
//' using Hermite polynomial based approximating function proposed by Gallant 
//' and Nychka in 1987. Please, see \code{\link[hpa]{dhpa}} 'Details' section to 
//' get more information concerning this approximating function.
//' @param selection an object of class "formula" 
//' (or one that can be coerced to that class): a symbolic description of the 
//' selection equation form. All variables in \code{selection} should be numeric 
//' vectors of the same length.
//' @param outcome an object of class "formula" (or one that can be coerced 
//' to that class): a symbolic description of the outcome equation form. 
//' All variables in \code{outcome} should be numeric vectors of the 
//' same length.
//' @template data_Template
//' @template z_K_Template
//' @template y_K_Template
//' @param pol_elements number of conditional expectation approximating terms 
//' for Newey's method. If \code{is_Newey_loocv} is \code{TRUE} then determines 
//' maximum number of these terms during leave-one-out cross-validation.
//' @param is_Newey logical; if TRUE then returns only Newey's method 
//' estimation results (default value is FALSE).
//' @param is_Newey_loocv logical; if TRUE then number of conditional 
//' expectation approximating terms for Newey's method will be selected
//' based on leave-one-out cross-validation criteria iterating through 0 
//' to pol_elements number of these terms.
//' @template x0_selection_Template
//' @template cov_type_Template
//' @template boot_iter_Template
//' @template is_parallel_Template
//' @template opt_type_Template
//' @template opt_control_Template
//' @template is_validation_Template
//' @template GN_details_Template
//' @template hpaSelection_formula_Template
//' @details Note that coefficient for the first
//' independent variable in \code{selection} will be fixed
//' to 1 i.e. \eqn{\gamma_{1}=1}.
//' @template is_numeric_selection_Template
//' @template parametric_paradigm_Template
//' @template Newey_details_Template
//' @details Note that selection equation dependent variables should have 
//' exactly two levels (0 and 1) where "0" states for the selection results 
//' which leads to unobservable values of dependent variable in 
//' outcome equation.
//' @template Mroz_reference_Template
//' @template optim_details_Template
//' @template opt_control_details_Template
//' @template opt_control_details_hpaSelection_Template
//' @return This function returns an object of class "hpaSelection".\cr \cr
//' An object of class "hpaSelection" is a list containing the 
//' following components:
//' \itemize{
//' \item \code{optim} - \code{\link[stats]{optim}} function output. 
//' If \code{opt_type = "GA"} then it is the list containing 
//' \code{\link[stats]{optim}} and \code{\link[GA]{ga}} functions outputs.
//' \item \code{x1} - numeric vector of distribution parameters estimates.
//' \item \code{Newey} - list containing information concerning Newey's 
//' method estimation results.
//' \item \code{z_mean} - estimate of the hermite polynomial mean parameter 
//' related to selection equation random error marginal distribution.
//' \item \code{y_mean} - estimate of the hermite polynomial mean parameter 
//' related to outcome equation random error marginal distribution.
//' \item \code{z_sd} - estimate of sd parameter related to selection equation 
//' random error marginal distribution.
//' \item \code{y_sd} - estimate of the hermite polynomial sd parameter related 
//' to outcome equation random error marginal distribution.
//' \item \code{pol_coefficients} - polynomial coefficients estimates.
//' \item \code{pol_degrees} - numeric vector which first element is \code{z_K} 
//' and the second is \code{y_K}.
//' \item \code{z_coef} - selection equation regression coefficients estimates.
//' \item \code{y_coef} - outcome equation regression coefficients estimates.
//' \item \code{cov_mat} - covariance matrix estimate.
//' \item \code{results} - numeric matrix representing estimation results.
//' \item \code{log-likelihood} - value of Log-Likelihood function.
//' \item \code{re_moments} - list which contains information about random 
//' errors expectations, variances and correlation.
//' \item \code{data_List} - list containing model variables and their 
//' partition according to outcome and selection equations.
//' \item \code{n_obs} - number of observations.
//' \item \code{ind_List} - list which contains information about parameters 
//' indexes in \code{x1}.
//' \item \code{selection_formula} - the same as \code{selection} 
//' input parameter.
//' \item \code{outcome_formula} - the same as \code{outcome} input parameter.}
//' Abovementioned list \code{Newey} has class "hpaNewey" and contains 
//' the following components:
//' \itemize{
//' \item \code{y_coef} - regression coefficients estimates (except 
//' constant term which is part of conditional expectation 
//' approximating polynomial).
//' \item \code{z_coef} - regression coefficients estimates related 
//' to selection equation.
//' \item \code{constant_biased} - biased estimate of constant term.
//' \item \code{inv_mills} - inverse mills ratios estimates and their 
//' powers (including constant).
//' \item \code{inv_mills_coef} - coefficients related to \code{inv_mills}.
//' \item \code{pol_elements} - the same as \code{pol_elements} 
//' input parameter. However if \code{is_Newey_loocv} is \code{TRUE}
//' then it will equal to the number of conditional expectation 
//' approximating terms for Newey's method which minimize leave-one-out 
//' cross-validation criteria.
//' \item \code{outcome_exp_cond} - dependent variable conditional 
//' expectation estimates.
//' \item \code{selection_exp} - selection equation random error 
//' expectation estimate.
//' \item \code{selection_var} - selection equation random error 
//' variance estimate.
//' \item \code{hpaBinaryModel} - object of class "hpaBinary" which 
//' contains selection equation estimation results.}
//' Abovementioned list \code{re_moments} contains the following components:
//' \itemize{
//' \item \code{selection_exp} - selection equation random errors 
//' expectation estimate.
//' \item \code{selection_var} - selection equation random errors 
//' variance estimate.
//' \item \code{outcome_exp} - outcome equation random errors 
//' expectation estimate.
//' \item \code{outcome_var} - outcome equation random errors 
//' variance estimate.
//' \item \code{errors_covariance} - outcome and selection equation 
//' random errors covariance estimate.
//' \item \code{rho} - outcome and selection equation random errors 
//' correlation estimate.
//' \item \code{rho_std} - outcome and selection equation random 
//' errors correlation estimator standard error estimate.}
//' @seealso \link[hpa]{summary.hpaSelection}, 
//' \link[hpa]{predict.hpaSelection}, \link[hpa]{plot.hpaSelection}, 
//' \link[hpa]{logLik.hpaSelection}
//' @template hpaSelection_examples_Template
//' @export	
// [[Rcpp::export]]
Rcpp::List hpaSelection(Rcpp::Formula selection,
	Rcpp::Formula outcome,
	DataFrame data,
	int z_K = 1,
	int y_K = 1,
	int pol_elements = 3,
	bool is_Newey = false,
	NumericVector x0 = NumericVector(0),
	bool is_Newey_loocv = false,
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
    
     // Validate polynomial degrees
     pol_Validate(NumericVector::create(z_K, y_K), NumericVector(0));
  }
  
	// Load additional environments

	  // stats environment
	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function optim = stats_env["optim"];
	Rcpp::Function model_frame = stats_env["model.frame"];
	Rcpp::Function na_omit_R = stats_env["na.omit"];
	Rcpp::Function na_pass = stats_env["na.pass"];
	Rcpp::Function complete_cases = stats_env["complete.cases"];
	Rcpp::Function cov_R = stats_env["cov"];

	  // base environment
	Rcpp::Environment base_env("package:base");
	Rcpp::Function solve = base_env["solve"];
	Rcpp::Function class_R = base_env["class"];
	Rcpp::Function c_R = base_env["c"];
	Rcpp::Function cbind_R = base_env["cbind"];
	Rcpp::Function subset_R = base_env["subset"];
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
	int pol_coefficients_n = (z_K + 1) * (y_K + 1) - 1;      // -1 because of 
	NumericVector pol_degrees = {(double)z_K, (double)y_K};  // restriction a(0...0) = 1
	NumericMatrix polIndex_mat = polynomialIndex(pol_degrees, 
                                               false);

	// Working with Data

		// Extract data frame from formula
	DataFrame z_df = model_frame(Rcpp::_["formula"] = selection, 
                               Rcpp::_["data"] = data, 
                               Rcpp::_["na.action"] = na_pass);
	
	DataFrame y_df = model_frame(Rcpp::_["formula"] = outcome, 
                               Rcpp::_["data"] = data, 
                               Rcpp::_["na.action"] = na_pass);

	DataFrame z_y_df = cbind_R(z_df, y_df);

	LogicalVector is_z_y_df_complete = complete_cases(z_y_df);
	LogicalVector is_z_df_complete = complete_cases(z_df);
	NumericVector z_temporal = z_y_df[0];
	LogicalVector is_y_unobs = (z_temporal == 0);

	LogicalVector df_cond = is_z_y_df_complete | 
	                        (is_z_df_complete & is_y_unobs);

	z_df = subset_R(Rcpp::_["x"] = z_df, 
                  Rcpp::_["subset"] = df_cond);
	
	y_df = subset_R(Rcpp::_["x"] = y_df, 
                  Rcpp::_["subset"] = df_cond);

	CharacterVector z_df_names = z_df.names();
	CharacterVector y_df_names = y_df.names();

	int z_df_n = z_df.size();
	int y_df_n = y_df.size();

	DataFrame y_df_no_y = y_df;

		// Extract dependend variables (regressors)
	NumericVector z = z_df[0]; // it is reference
	NumericVector y = y_df[0]; // it is reference

	int n = z.size();

		// Extract independend variable
	NumericMatrix z_d(n, z_df_n - 1); // -1 because there is no constant term
	NumericMatrix y_d(n, y_df_n - 1); // -1 because there is no constant term

	int z_d_col = z_d.ncol();
	int y_d_col = y_d.ncol();

	for (int i = 0; i < z_d_col; i++)
	{
		z_d(_, i) = NumericVector(z_df[i+1]);
	}
	for (int i = 0; i < y_d_col; i++)
	{
		y_d(_, i) = NumericVector(y_df[i + 1]);
	}

	// Create initial values vector
	bool x0_given = true;

	if (x0.size() == 0)
	{
		x0_given = false;
		x0 = NumericVector(pol_coefficients_n + 3 + z_d_col + y_d_col); // +2 for mean 
	}                                                                 // and sd

	// Assign indexes

		// Initialize additional index and upper value for some loops
	int lower_ind = 0;
	int upper_ind = pol_coefficients_n;

		// for polynomial coefficients
	NumericVector pol_coefficients_ind(pol_coefficients_n);

	for (int i = lower_ind; i < upper_ind; i++)
	{
		pol_coefficients_ind[i] = i;
	}

		// for mean vector
	int z_mean_ind = pol_coefficients_n;
	int y_mean_ind = z_mean_ind + 1;

		// for sd vector
	int z_sd_ind = y_mean_ind + 1;
	int y_sd_ind = z_sd_ind + 1;

		// for z coefficients
	NumericVector z_coef_ind(z_d_col - 1);

	lower_ind = y_sd_ind + 1;
	upper_ind = lower_ind + z_d_col - 1;

	for (int i = lower_ind; i < upper_ind; i++)
	{
		z_coef_ind[i - lower_ind] = i;
	}

		// for y coefficients
	NumericVector y_coef_ind(y_d_col);

	lower_ind = upper_ind;
	upper_ind = lower_ind + y_d_col;

	for (int i = lower_ind; i < upper_ind; i++)
	{
		y_coef_ind[i - lower_ind] = i;
	}

	// Convert to arma
	arma::vec y_arma = as<arma::vec>(y);
	arma::vec z_arma = as<arma::vec>(z);
	arma::mat y_d_arma = as<arma::mat>(y_d);
	arma::mat z_d_arma = as<arma::mat>(z_d);

	// Divide into observable and unobservable samples

		// observable
	arma::vec y_1 = as<arma::vec>(y[z == 1]);
	int n_1 = y_1.size();
	arma::mat y_d_1 = (as<arma::mat>(y_d)).rows(arma::find(z_arma == 1));
	arma::vec z_1 = as<arma::vec>(z[z == 1]);
	arma::mat z_d_1 = (as<arma::mat>(z_d)).rows(arma::find(z_arma == 1));

		// unobservable
	arma::vec y_0 = as<arma::vec>(y[z == 0]);
	arma::mat y_d_0 = (as<arma::mat>(y_d)).rows(arma::find(z_arma == 0));
	arma::vec z_0 = as<arma::vec>(z[z == 0]);
	arma::mat z_d_0 = (as<arma::mat>(z_d)).rows(arma::find(z_arma == 0));

	// Create List to store Newey's method post estimates
	List Newey;
	
	// Get initial values from hpaBinary and Newey
	if (!x0_given | is_Newey)
	{
		// Estimate selection equation parameters via hpaBinary
		List modelBinary;
		try
		{
			modelBinary = hpaBinary(selection,
				data, z_K,
				NA_REAL, NA_REAL, 0,
				true, true, false, NumericVector(0), 
				"sandwich", 100, is_parallel, "optim", R_NilValue, false);
		} catch (std::exception &ex) {
			warning("Can't get initial values from semi-nonparametric binary choice model");
			forward_exception_to_r(ex);
		}

		// Store hpaBinary estimates to x0

		  // store polynomial coefficients
		NumericVector z_pol_coef_temporal = modelBinary["pol_coefficients"];

		int z_pol_ind = 1;
		for (int i = 1; i < pol_coefficients_n; i++)
		{
			if (polIndex_mat(1, i) == 0)
			{
				x0[pol_coefficients_ind[i - 1]] = z_pol_coef_temporal[z_pol_ind];
				z_pol_ind++;
			}
		}

		  // store mean and sd estimates
		  x0[z_mean_ind] = modelBinary["mean"];
		  x0[z_sd_ind] = modelBinary["sd"];

		  // store coefficients
		NumericVector z_coef_temporal = modelBinary["coefficients"];
		NumericVector z_coef_ind_temporal = z_coef_temporal[Rcpp::Range(
		                                                    1, z_d_col - 1)];
		x0[z_coef_ind] = z_coef_ind_temporal;

		  // get latent variable value
		NumericVector z_latent = wrap(z_d_arma * as<arma::vec>(z_coef_temporal));
		double z_exp = modelBinary["errors_exp"];
		double z_var = modelBinary["errors_var"];
		NumericVector z_latent_std = (z_latent + z_exp) / sqrt(z_var); // standardize for
		                                                               // mills ratio

		  // separate latent variable values depending on the outcome
		NumericVector z_latent_1 = z_latent[z == 1];
		NumericVector z_latent_std_1 = z_latent_std[z == 1];

		  // estimate mills ratios
		NumericVector z_mills_std_dnorm = dnorm(z_latent_std_1);
		NumericVector z_mills_std_pnorm = pnorm(z_latent_std_1);
		NumericVector z_mills_base = z_mills_std_dnorm / z_mills_std_pnorm;
		
		NumericMatrix z_mills = NumericMatrix(n_1, pol_elements + 1);

		for (int i = 0; i < (pol_elements + 1); i++)
		{
			z_mills(_, i) = pow(z_mills_base, i);
		}

		// Newey method with pol_elements 
		// approximating polynomial elements
		int pol_elements_i = pol_elements;

		if (is_Newey_loocv)
		{
			pol_elements_i = 0;
		}
		
			// Values to return from cycle based on loocv
		float loocv_best = R_PosInf;
		float pol_elements_best = pol_elements_i;
		NumericVector y_coef_Newey;
		NumericVector residuals_ls;

			// Iterate through all polynomial degrees
		for(int t = pol_elements_i; t <= pol_elements; t++)
		{
			// Prepare the data
			arma::mat y_d_1_Newey = arma::mat(n_1, y_d_col + t + 1, arma::fill::ones);

			for (int i = 0; i < y_d_col; i++)
			{
				y_d_1_Newey.col(i) = y_d_1.col(i);
			}

			for (int i = y_d_col; i < (y_d_col + t + 1); i++)
			{
				NumericVector z_mills_temporal = z_mills(_, i - y_d_col);
				y_d_1_Newey.col(i) = as<arma::vec>(z_mills_temporal);
			}

			  // Estimate hat matrix
			NumericMatrix hat_Newey_cycle = wrap(y_d_1_Newey * 
			                                inv(y_d_1_Newey.t() * y_d_1_Newey) * 
			                                y_d_1_Newey.t());

			  // Estimate coefficients
			NumericVector y_coef_Newey_cycle = wrap(inv(y_d_1_Newey.t() * 
			                                            y_d_1_Newey) * 
			                                            y_d_1_Newey.t() * y_1);

			  // Estimate loocv
			NumericVector residuals_ls_cycle = wrap(
			  y_d_1_Newey * as<arma::vec>(y_coef_Newey_cycle) - y_1);

			float loocv = 0.0;

			for (int j = 0; j < n_1; j++)
			{
				float val_1 = residuals_ls_cycle[j];
				float val_2 = hat_Newey_cycle(j, j);
				float val_ratio = val_1 / (1 - val_2);
				loocv += val_ratio * val_ratio;
			}

			loocv /= (float)n_1;
			
			if (loocv < loocv_best)
			{
				loocv_best = loocv;
				pol_elements_best = t;
				y_coef_Newey = y_coef_Newey_cycle;
				residuals_ls = residuals_ls_cycle;
			}

		}

		// Assign values from the cycle

		  // Assign coefficients to x0
		NumericVector y_coef_Newey_temporal = y_coef_Newey[Rcpp::Range(
		                                                   0, y_d_col - 1)];
		x0[y_coef_ind] = y_coef_Newey_temporal;
		x0[y_mean_ind] = y_coef_Newey[y_d_col];

		arma::mat residuals_ls_arma = as<arma::mat>(residuals_ls);
		arma::mat residuals_squared = residuals_ls_arma.t() * residuals_ls_arma;

		x0[y_sd_ind] = accu(sqrt(residuals_squared / (n_1 - 1 - y_d_col)));

		  // Get additional values for inverse mills ratios
		arma::mat y_d_1_Newey = arma::mat(n_1, y_d_col + pol_elements_best + 
		                                       1, arma::fill::ones);
		NumericVector y_coef_Newey_mills = y_coef_Newey[Rcpp::Range(
		                                                y_d_col, y_d_col + 
		                                                pol_elements_best)];
		NumericVector z_expect = NumericVector(n_1);

		for (int i = 0; i < (pol_elements_best + 1); i++)
		{
		  NumericVector z_mills_temporal = z_mills(_, i);
		  y_d_1_Newey.col(i) = as<arma::vec>(z_mills_temporal);
			z_expect = z_expect + y_coef_Newey_mills[i] * z_mills(_, i);
		}
		
		  // Summarize the output for Newey's method
		Newey = List::create(Named("y_coef") = y_coef_Newey_temporal,
			Named("z_coef") = z_coef_ind_temporal,
			Named("constant_biased") = x0[y_mean_ind], 
			Named("sd_biased") = x0[y_sd_ind],  // substitute for cov_mat in future
			Named("inv_mills") = z_mills,
			Named("inv_mills_coef") = y_coef_Newey_mills,
			Named("pol_elements") = pol_elements_best,
			Named("outcome_exp_cond") = z_expect,
			Named("selection_exp") = z_exp,
			Named("selection_var") = z_var,
			Named("z_mean") = x0[z_mean_ind],
			Named("z_sd") = x0[z_sd_ind],
			Named("y_mean") = x0[y_mean_ind],
			Named("y_sd") = x0[y_sd_ind],
			Named("hpaBinaryModel") = modelBinary);

		Newey.attr("class") = "hpaNewey";

		if (is_Newey)
		{
			return(Newey);
		}
	}

	// Create list for some variables because unfortunately optim 
	// function has limitation for the parameters number 
	// (parameters itself not estimated)
	
	  // Collect some values to lists since there are 
	  // limited number of objects could be stored in 
	  // the list in Rcpp
	List ind_List = List::create(
	  Named("pol_coefficients_ind") = pol_coefficients_ind,
		Named("z_mean_ind") = z_mean_ind,
		Named("y_mean_ind") = y_mean_ind,
		Named("y_sd_ind") = y_sd_ind,
		Named("z_sd_ind") = z_sd_ind,
		Named("y_coef_ind") = y_coef_ind,
		Named("z_coef_ind") = z_coef_ind);

	  // Store all the values into the list
	List hpaSelection_args = List::create(
	  Named("ind_List") = ind_List,
    Named("z_0") = z_0, Named("z_1") = z_1,
    Named("y_0") = y_0, Named("y_1") = y_1,
    Named("y_d_1") = y_d_1, Named("y_d_0") = y_d_0,
    Named("z_d_1") = z_d_1, Named("z_d_0") = z_d_0,
    Named("pol_degrees") = pol_degrees,
    Named("is_parallel") = is_parallel);

	// Apply optimization routine
	List PGN_control = List::create(
	  Named("maxit") = 10000000, 
	  Named("fnscale") = -1.0,
	  Named("abstol") = std::sqrt(std::numeric_limits<double>::epsilon()) * 0.01,
	  Named("reltol") = std::sqrt(std::numeric_limits<double>::epsilon()) * 0.01);

	List optim_results = optim(
		Rcpp::_["par"] = x0,
		Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaSelectionLnLOptim),
		Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaSelectionLnLOptim_grad),
		Rcpp::_["control"] = PGN_control,
		Rcpp::_["method"] = "BFGS",
		Rcpp::_["hessian"] = true,
		Rcpp::_["hpaSelection_args"] = hpaSelection_args);

	// Extract coefficients and function value
	
	NumericVector x1 = optim_results["par"];
	
	int x1_n = x0.size();
	
	// Evolutionary algorithm
	
	List ga_List;
	List ga_summary;
	
	// set lower and upper bounds for parameters space
	NumericVector ga_lower = NumericVector(x1_n);
	NumericVector ga_upper = NumericVector(x1_n);
	
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
	  
	  if(opt_control.containsElementNamed("lower") & 
       opt_control.containsElementNamed("upper"))
	  {
	    ga_lower = opt_control["lower"];
	    ga_upper = opt_control["upper"];
	  } else {

	    // Set lower and upper bounds for parameters space
	    
	      // bounds for the z_mean parameter
  	  double ga_lower_z_mean = x1[z_mean_ind] - 2 * abs(x1[z_mean_ind]);
  	  double ga_upper_z_mean = x1[z_mean_ind] + 2 * abs(x1[z_mean_ind]);
  	    
  	  ga_lower[z_mean_ind] = ga_lower_z_mean;
  	  ga_upper[z_mean_ind] = ga_upper_z_mean;
  	  
  	    // bounds for the y_mean parameter
  	  double ga_lower_y_mean = x1[y_mean_ind] - 2 * abs(x1[y_mean_ind]);
  	  double ga_upper_y_mean = x1[y_mean_ind] + 2 * abs(x1[y_mean_ind]);
  	  
  	  ga_lower[y_mean_ind] = ga_lower_y_mean;
  	  ga_upper[y_mean_ind] = ga_upper_y_mean;
  	  
  	    // bounds for the y_sd parameter
  	  double ga_lower_y_sd = x1[y_sd_ind] * 0.2;
  	  double ga_upper_y_sd = x1[y_sd_ind] * 5;

  	  ga_lower[y_sd_ind] = ga_lower_y_sd;
  	  ga_upper[y_sd_ind] = ga_upper_y_sd;

  	    // bounds for the z regression coefficients
  	  for(int i = z_coef_ind[0]; i < (z_coef_ind[0] + z_coef_ind.size()); i++)
  	  {
  	    double z_coef_ga = x1[i];
  	    ga_lower[i] = z_coef_ga - 2 * z_coef_ga;
  	    ga_upper[i] = z_coef_ga + 2 * z_coef_ga;
  	  }

  	    // bounds for the y regression coefficients
  	  for(int i = y_coef_ind[0]; i < (y_coef_ind[0] + y_coef_ind.size()); i++)
  	  {
  	    double y_coef_ga = x1[i];
  	    ga_lower[i] = y_coef_ga - 2 * y_coef_ga;
  	    ga_upper[i] = y_coef_ga + 2 * y_coef_ga;
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
	    ga_maxiter = 50 * (1 + z_K) * (1 + y_K) + 10 * (z_d_col + y_d_col);
	  }
	  
  	  // set population size
	  int ga_popSize;
	  if(opt_control.containsElementNamed("popSize"))
	  {
	    ga_popSize = opt_control["popSize"];
	  } else {
	    ga_popSize = 10 + (1 + z_K) * (1 + y_K) * 5 + 2 * (z_d_col + y_d_col);
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
	    Rcpp::_["fitness"] = Rcpp::InternalFunction(&hpaSelectionLnLOptim),
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
	    Rcpp::_["hpaSelection_args"] = hpaSelection_args,
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
	    Rcpp::_["par"] = x0,
	    Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaSelectionLnLOptim),
	    Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaSelectionLnLOptim_grad),
	    Rcpp::_["control"] = PGN_control,
	    Rcpp::_["method"] = "BFGS",
	    Rcpp::_["hessian"] = true,
	    Rcpp::_["hpaSelection_args"] = hpaSelection_args);
	  
	  // Reset optimal point
	  x1 = optim_results["par"];
	}

	// Calculate additional values

		// calculate log-likelihood, AIC and BIC
	double lnL = optim_results["value"];

		// get polynomial coefficients
	NumericVector pol_coefficients = NumericVector(pol_coefficients_n);

	if ((z_K != 0) | (y_K != 0))
	{
		pol_coefficients = x1[pol_coefficients_ind];
	}

	pol_coefficients.push_front(1);

		// get mean and sd values
	double z_mean = x1[z_mean_ind];
	double y_mean = x1[y_mean_ind];
	double y_sd = x1[y_sd_ind];
	double z_sd = x1[z_sd_ind];

		// get coefficients
	NumericVector z_coef = x1[z_coef_ind];
	z_coef.push_front(1);
	NumericVector y_coef = x1[y_coef_ind];

	// Get covariance matrix estimate of "cov_type" type
	
	NumericMatrix cov_mat;
	
	arma::mat H_part;
	arma::mat J_part;
	
	// Estimate jacobian for inner part
	if ((cov_type == "gop") | (cov_type == "sandwich") | 
      (cov_type == "sandwichFD"))
	{
	  NumericMatrix my_jacobian = hpaSelectionLnLOptim_grad_ind(
	    x1, hpaSelection_args);
	  
	  J_part = as<arma::mat>(my_jacobian);
	}
	
	// Estimate hessian matrix
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
	    my_hessian = hpaSelectionLnLOptim_hessian(x1, hpaSelection_args);

	  } catch (std::exception &ex) {
	    warning("Can't calculate Hessian via first difference method. Hessian from the optim function will be used instead.");
	    forward_exception_to_r(ex);
	  }
	  
	  H_part = as<arma::mat>(my_hessian).i();
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
	
	// Gradient outer product estimate
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
	NumericVector sample_ind = NumericVector(n);
	
	  // list to store bootstrap results
	List boot_List;
	
	if (cov_type == "bootstrap")
	{
	  for(int i = 0; i < boot_iter; i++)
	  {
	    // Generate sample with replacement
	    NumericVector sample_ind = floor(runif(n, 0, n));
	    
	    NumericVector z_boot = NumericVector(n);
	    NumericVector y_boot = NumericVector(n);
	    
	    NumericMatrix z_d_boot = NumericMatrix(z_d.nrow(), z_d.ncol());
	    NumericMatrix y_d_boot = NumericMatrix(y_d.nrow(), y_d.ncol());
	    
	    for (int j = 0; j < n; j++)
	    {
	      z_boot[j] = z[sample_ind[j]];
	      y_boot[j] = y[sample_ind[j]];
	      
	      z_d_boot(j, _) = z_d(sample_ind[j], _);
	      y_d_boot(j, _) = y_d(sample_ind[j], _);
	    }
	    
	    // Convert to arma
	    arma::vec z_arma_boot = as<arma::vec>(z_boot);
	    arma::vec y_arma_boot = as<arma::vec>(y_boot);
	    
	    arma::mat z_d_arma_boot = as<arma::mat>(z_d_boot);
	    arma::mat y_d_arma_boot = as<arma::mat>(y_d_boot);
	    
	    // Divide into 0 and 1 samples
	    arma::vec z_1_boot = as<arma::vec>(z_boot[z_boot == 1]);
	    arma::vec y_1_boot = as<arma::vec>(y_boot[z_boot == 1]);
	    
	    arma::mat z_d_1_boot = (as<arma::mat>(z_d_boot)).rows(
	      arma::find(z_arma_boot == 1));
	    arma::mat y_d_1_boot = (as<arma::mat>(y_d_boot)).rows(
	      arma::find(z_arma_boot == 1));
	    
	    arma::vec z_0_boot = as<arma::vec>(z_boot[z_boot == 0]);
	    arma::vec y_0_boot = as<arma::vec>(y_boot[z_boot == 0]);
	    
	    arma::mat z_d_0_boot = (as<arma::mat>(z_d_boot)).rows(
	      arma::find(z_arma_boot == 0));
	    arma::mat y_d_0_boot = (as<arma::mat>(y_d_boot)).rows(
	      arma::find(z_arma_boot == 0));
	    
	    // Prepare arguments for bootstrap List
	    List hpaSelection_args_boot = hpaSelection_args;
	    hpaSelection_args_boot["z_1"] = z_1_boot;
	    hpaSelection_args_boot["z_0"] = z_0_boot;
	    hpaSelection_args_boot["z_d_1"] = z_d_1_boot;
	    hpaSelection_args_boot["z_d_0"] = z_d_0_boot;
	    hpaSelection_args_boot["y_1"] = y_1_boot;
	    hpaSelection_args_boot["y_0"] = y_0_boot;
	    hpaSelection_args_boot["y_d_1"] = y_d_1_boot;
	    hpaSelection_args_boot["y_d_0"] = y_d_0_boot;
	    
	    // Perform estimaton
	    List boot_results = optim(
	      Rcpp::_["par"] = x1,
	      Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaSelectionLnLOptim),
	      Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaSelectionLnLOptim_grad),
	      Rcpp::_["control"] = PGN_control,
	      Rcpp::_["method"] = "BFGS",
	      Rcpp::_["hessian"] = true,
	      Rcpp::_["hpaSelection_args"] = hpaSelection_args_boot);
	    
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

	// polIndex_mat

		// for alpha
	double z_stat = 0;
	NumericVector F_z_stat;
	for (int i = 1; i <= pol_coefficients_n; i++)
	{
		results_rows[(i - 1)] = "a_" + std::to_string((int)polIndex_mat(0, i)) + 
		                        "_" + std::to_string((int)polIndex_mat(1, i));
		results((i - 1), 0) = pol_coefficients[(i - 1) + 1];
		results((i - 1), 1) = sqrt(cov_mat((i - 1), (i - 1)));
		z_stat = results((i - 1), 0) / results((i - 1), 1);
		F_z_stat = pnorm(NumericVector::create(z_stat));
		results((i - 1), 2) = z_stat;
		results((i - 1), 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

		// for z_mean
	results_rows[z_mean_ind] = "z_mean";
	results(z_mean_ind, 0) = x1[z_mean_ind];
	results(z_mean_ind, 1) = sqrt(cov_mat((z_mean_ind - 1), (z_mean_ind - 1)));
	z_stat = results(z_mean_ind, 0) / results(z_mean_ind, 1);
	F_z_stat = pnorm(NumericVector::create(z_stat));
	results(z_mean_ind, 2) = z_stat;
	results(z_mean_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);

		// for y_mean
	results_rows[y_mean_ind] = "y_mean";
	results(y_mean_ind, 0) = x1[y_mean_ind];
	results(y_mean_ind, 1) = sqrt(cov_mat((y_mean_ind), (y_mean_ind)));
	z_stat = results(y_mean_ind, 0) / results(y_mean_ind, 1);
	F_z_stat = pnorm(NumericVector::create(z_stat));
	results(y_mean_ind, 2) = z_stat;
	results(y_mean_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	
	  // for z_sd
	results_rows[z_sd_ind] = "z_sd";
	results(z_sd_ind, 0) = x1[z_sd_ind];
	results(z_sd_ind, 1) = sqrt(cov_mat((z_sd_ind), (z_sd_ind)));
	z_stat = results(z_sd_ind, 0) / results(z_sd_ind, 1);
	F_z_stat = pnorm(NumericVector::create(z_stat));
	results(z_sd_ind, 2) = z_stat;
	results(z_sd_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);

		// for y_sd
	results_rows[y_sd_ind] = "y_sd";
	results(y_sd_ind, 0) = x1[y_sd_ind];
	results(y_sd_ind, 1) = sqrt(cov_mat((y_sd_ind), (y_sd_ind)));
	z_stat = results(y_sd_ind, 0) / results(y_sd_ind, 1);
	F_z_stat = pnorm(NumericVector::create(z_stat));
	results(y_sd_ind, 2) = z_stat;
	results(y_sd_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);

		// for z coefficients
	for (int i = 0; i < (z_d_col - 1); i++)
	{
		results_rows[z_coef_ind[i]] = z_df_names(i + 2);
		results(z_coef_ind[i], 0) = x1[z_coef_ind[i]];
		results(z_coef_ind[i], 1) = sqrt(cov_mat(z_coef_ind[i], z_coef_ind[i]));
		z_stat = results(z_coef_ind[i], 0) / results(z_coef_ind[i], 1);
		F_z_stat = pnorm(NumericVector::create(z_stat));
		results(z_coef_ind[i], 2) = z_stat;
		results(z_coef_ind[i], 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

		// for y coefficients
	for (int i = 0; i < y_d_col; i++)
	{
		results_rows[y_coef_ind[i]] = y_df_names(i + 1);
		results(y_coef_ind[i], 0) = x1[y_coef_ind[i]];
		results(y_coef_ind[i], 1) = sqrt(cov_mat(y_coef_ind[i], y_coef_ind[i]));
		double z_stat = results(y_coef_ind[i], 0) / results(y_coef_ind[i], 1);
		F_z_stat = pnorm(NumericVector::create(z_stat));
		results(y_coef_ind[i], 2) = z_stat;
		results(y_coef_ind[i], 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

		// assign names to the output
	rownames(results) = results_rows;
	colnames(results) = results_cols;

	x1.names() = results_rows;

	rownames(cov_mat) = results_rows;
	colnames(cov_mat) = results_rows;

	if ((z_K != 0) & (y_K != 0))
	{
		pol_coefficients.names() = c_R("a_0", results_rows[pol_coefficients_ind]);
	}

	z_coef.names() = results_rows[z_coef_ind];
	y_coef.names() = results_rows[y_coef_ind];

	// Calculate expectation and variance
	NumericVector z_e = ehpa(NumericMatrix(1, 1), pol_coefficients, pol_degrees,
		LogicalVector{false, false}, LogicalVector{false, false},
		NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, y_sd),
		NumericVector::create(1, 0), is_parallel, false);

	NumericVector z_e_2 = ehpa(NumericMatrix(1, 1), pol_coefficients, pol_degrees,
		LogicalVector{false, false}, LogicalVector{false, false},
		NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, y_sd),
		NumericVector::create(2, 0), is_parallel, false);

	NumericVector y_e = ehpa(NumericMatrix(1, 1), pol_coefficients, pol_degrees,
		LogicalVector{false, false}, LogicalVector{false, false},
		NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, y_sd),
		NumericVector::create(0, 1), is_parallel, false);

	NumericVector y_e_2 = ehpa(NumericMatrix(1, 1), pol_coefficients, pol_degrees,
		LogicalVector{false, false}, LogicalVector{false, false},
		NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, y_sd),
		NumericVector::create(0, 2), is_parallel, false);

	NumericVector z_y_e = ehpa(NumericMatrix(1, 1), pol_coefficients, pol_degrees,
		LogicalVector{false, false}, LogicalVector{false, false},
		NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, y_sd),
		NumericVector::create(1, 1), is_parallel, false);

	double z_v = z_e_2[0] - z_e[0] * z_e[0];
	double y_v = y_e_2[0] - y_e[0] * y_e[0];
	double z_y_c = z_y_e[0] - z_e[0] * y_e[0];
	double rho = z_y_c / sqrt(z_v * y_v);
	
	// apply delta method in order to estimate correlation estimator's
	// standard error
	
	  // apply numeric differentiation in order to calculate gradient
	NumericVector rho_grad = NumericVector(x1_n);
	
	NumericVector x1_eps = clone(x1);
	
	double machinePrecision = std::numeric_limits<double>::epsilon();
	double my_precision = std::sqrt(machinePrecision);
	
  for(int i = 0; i <= pol_coefficients_n; i++)
  {
    double eps = std::abs(x1[i] * my_precision);

    x1_eps[i] = x1[i] + eps;

    NumericVector pol_coefficients_eps = x1_eps[pol_coefficients_ind];
    pol_coefficients_eps.push_front(1);
    
    NumericVector z_e_eps = ehpa(NumericMatrix(1, 1), pol_coefficients_eps, 
                                 pol_degrees,
                                 LogicalVector{false, false}, 
                                 LogicalVector{false, false},
                                 NumericVector::create(x1_eps[z_mean_ind], 
                                                       x1_eps[y_mean_ind]), 
                                 NumericVector::create(z_sd, x1_eps[y_sd_ind]),
                                 NumericVector::create(1, 0), 
                                 is_parallel, false);
    
    NumericVector z_e_2_eps = ehpa(NumericMatrix(1, 1), 
                                   pol_coefficients_eps, pol_degrees,
                                   LogicalVector{false, false}, 
                                   LogicalVector{false, false},
                                   NumericVector::create(x1_eps[z_mean_ind], 
                                                         x1_eps[y_mean_ind]), 
                                   NumericVector::create(z_sd, 
                                                         x1_eps[y_sd_ind]),
                                   NumericVector::create(2, 0), 
                                   is_parallel, false);
    
    NumericVector y_e_eps = ehpa(NumericMatrix(1, 1), 
                                 pol_coefficients_eps, pol_degrees,
                                 LogicalVector{false, false}, 
                                 LogicalVector{false, false},
                                 NumericVector::create(x1_eps[z_mean_ind], 
                                                       x1_eps[y_mean_ind]), 
                                 NumericVector::create(z_sd, 
                                                       x1_eps[y_sd_ind]),
                                 NumericVector::create(0, 1), 
                                 is_parallel, false);
    
    NumericVector y_e_2_eps = ehpa(NumericMatrix(1, 1), 
                                   pol_coefficients_eps, pol_degrees,
                                   LogicalVector{false, false}, 
                                   LogicalVector{false, false},
                                   NumericVector::create(x1_eps[z_mean_ind], 
                                                         x1_eps[y_mean_ind]), 
                                   NumericVector::create(z_sd, 
                                                         x1_eps[y_sd_ind]),
                                   NumericVector::create(0, 2),
                                   is_parallel, false);
    
    NumericVector z_y_e_eps = ehpa(NumericMatrix(1, 1), 
                                   pol_coefficients_eps, pol_degrees,
                                   LogicalVector{false, false}, 
                                   LogicalVector{false, false},
                                   NumericVector::create(x1_eps[z_mean_ind], 
                                                         x1_eps[y_mean_ind]), 
                                   NumericVector::create(z_sd, 
                                                         x1_eps[y_sd_ind]),
                                   NumericVector::create(1, 1), 
                                   is_parallel, false);
    
    double z_v_eps = z_e_2_eps[0] - z_e_eps[0] * z_e_eps[0];
    double y_v_eps = y_e_2_eps[0] - y_e_eps[0] * y_e_eps[0];
    double z_y_c_eps = z_y_e_eps[0] - z_e_eps[0] * y_e_eps[0];
    double rho_eps = z_y_c_eps / sqrt(z_v_eps * y_v_eps);

    rho_grad[i] = (rho_eps - rho) / eps;

    x1_eps[i] = x1[i];
  }

  arma::mat rho_grad_arma = as<arma::vec>(rho_grad);
  arma::mat cov_mat_arma = as<arma::mat>(cov_mat);
  
  NumericVector rho_var = wrap(rho_grad_arma.t() * 
                               cov_mat_arma * rho_grad_arma);
  
  // store moments information

	List re_moments = List::create(Named("optim") = optim_results, 
		Named("selection_exp") = z_e,
		Named("outcome_exp") = y_e,
		Named("selection_var") = z_v,
		Named("outcome_var") = y_v,
		Named("errors_covariance") = z_y_c,
		Named("rho") = rho,
		Named("rho_std") = sqrt(rho_var[0]));

	List data_List = List::create(Named("data_z") = z_df,
		Named("data_y") = y_df,
		Named("dataframe") = data);
	
	if(opt_type == "GA")
	{
	  optim_results = List::create(
	    Named("optim") = optim_results,
	    Named("GA") = ga_List,
	    Named("GA_summary") = ga_summary);
	}

	List return_result = List::create(
	  Named("optim") = optim_results,
		Named("x1") = x1,
		Named("Newey") = Newey,
		Named("log-likelihood") = lnL,
		Named("cov_mat") = cov_mat,
		Named("n_obs") = n,
		Named("data_List") = data_List,
		Named("results") = results,
		Named("z_mean") = z_mean,
		Named("y_mean") = y_mean, 
		Named("z_sd") = z_sd, 
		Named("y_sd") = y_sd,
		Named("pol_coefficients") = pol_coefficients,
		Named("pol_degrees") = pol_degrees,
		Named("y_coef") = x1[y_coef_ind],
		Named("z_coef") = z_coef,
		Named("selection_formula") = selection,
		Named("outcome_formula") = outcome,
		Named("re_moments") = re_moments,
		Named("ind_List") = ind_List);

	return_result.attr("class") = "hpaSelection";

	return(return_result);
}

// Perform log-likelihood function estimation for selection model
List hpaSelectionLnLOptim_List(NumericVector x0, List hpaSelection_args)
{
  // Get values from the hpaSelection_args list
  List ind_List = Rcpp::as<Rcpp::List>(hpaSelection_args["ind_List"]);
  arma::vec z_0 = hpaSelection_args["z_0"];
  arma::vec z_1 = hpaSelection_args["z_1"];
  arma::vec y_0 = hpaSelection_args["y_0"];
  arma::vec y_1 = hpaSelection_args["y_1"];
  arma::mat z_d_0 = hpaSelection_args["z_d_0"];
  arma::mat z_d_1 = hpaSelection_args["z_d_1"];
  arma::mat y_d_0 = hpaSelection_args["y_d_0"];
  arma::mat y_d_1 = hpaSelection_args["y_d_1"];
  NumericVector pol_degrees = hpaSelection_args["pol_degrees"];
  bool is_parallel = hpaSelection_args["is_parallel"];
  
    // Get values from the ind_List
  NumericVector pol_coefficients_ind = ind_List["pol_coefficients_ind"];
  int z_mean_ind = ind_List["z_mean_ind"];
  int y_mean_ind = ind_List["y_mean_ind"];
  int y_sd_ind = ind_List["y_sd_ind"];
  int z_sd_ind = ind_List["z_sd_ind"];
  NumericVector z_coef_ind = ind_List["z_coef_ind"];
  NumericVector y_coef_ind = ind_List["y_coef_ind"];
  
  // Assign estimated parameters values to corresponding vectors
  
    // polynomial coefficients and degrees
  NumericVector pol_coefficients = x0[pol_coefficients_ind];
  pol_coefficients.push_front(1); //add alpha(0...0)
  
    // mean
  double z_mean = x0[z_mean_ind];
  double y_mean = x0[y_mean_ind];
  NumericVector mean = {z_mean, y_mean};
  
    // sd
  double y_sd = x0[y_sd_ind];
  double z_sd = x0[z_sd_ind];
  NumericVector sd = {z_sd, y_sd};
  
  if (y_sd <= 0)
  {
    y_sd = std::sqrt(std::numeric_limits<double>::epsilon());
  }
  
  if (z_sd <= 0)
  {
    z_sd = std::sqrt(std::numeric_limits<double>::epsilon());
  }

  // coefficients for independend variables
  arma::vec y_coef = as<arma::vec>(x0[y_coef_ind]);
  NumericVector z_coef_temporal = x0[z_coef_ind];
  z_coef_temporal.push_front(1); //for fixed coefficient
  arma::vec z_coef = as<arma::vec>(z_coef_temporal);

  // get estimates for z*
  NumericVector z_h_1 = wrap(z_d_1 * z_coef);
  NumericVector z_h_0 = wrap(z_d_0 * z_coef);
  
  // get estimates for y and random errors
  arma::mat y_h_1 = y_d_1 * y_coef;
  
  NumericVector e_h_1 = wrap(y_1 - y_h_1);
  
  // concatenate e_h and z_h and prepare to insert into function
  NumericMatrix z_y_1 = NumericMatrix(z_h_1.size(), 2);
  NumericMatrix z_y_0 = NumericMatrix(z_h_0.size(), 2);
  
  z_y_1(_, 0) = -1.0 * z_h_1;
  z_y_1(_, 1) = e_h_1;
  z_y_0(_, 0) = -1.0 * z_h_0;
  
  // get number of observations with 0 and 1 values
  int n_obs_0 = z_h_0.size();
  int n_obs_1 = z_h_1.size();
  
  // lower tail negative infinity matrix
  NumericMatrix inf_y_vec_1 = NumericMatrix(n_obs_1, 2);
  std::fill(inf_y_vec_1.begin(), inf_y_vec_1.end(), R_PosInf);
  inf_y_vec_1(_, 1) = e_h_1; // in order to condition on e_h_1
  NumericMatrix neg_inf_vec_0 = NumericMatrix(n_obs_0, 2);
  std::fill(neg_inf_vec_0.begin(), neg_inf_vec_0.end(), R_NegInf);
  
  // likelihood calculation
  NumericVector lnL_y_1 = dhpa(z_y_1,
                               pol_coefficients, pol_degrees,
                               LogicalVector{false, false}, 
                               LogicalVector{true, false},
                               mean, sd,
                               is_parallel, true, false);
  
  NumericVector lnL_z_y_1 = ihpa(z_y_1, inf_y_vec_1,
                                 pol_coefficients, pol_degrees,
                                 LogicalVector{false, true}, 
                                 LogicalVector{false, false},
                                 mean, sd,
                                 is_parallel, true, false);
  
  NumericVector lnL_z_y_0 = ihpa(neg_inf_vec_0, z_y_0,
                                 pol_coefficients, pol_degrees,
                                 LogicalVector{false, false}, 
                                 LogicalVector{false, true},
                                 mean, sd,
                                 is_parallel, true, false);

  // Initialize list to store calculation results
  double aggregate_y_1 = 0.0;
  double aggregate_z_y_1 = 0.0;
  double aggregate_z_y_0 = 0.0;
  
  List return_List = List::create(
    Named("individual_y_1") = NumericVector::create(0.0),
    Named("individual_z_y_1") = NumericVector::create(0.0),
    Named("individual_z_y_0") = NumericVector::create(0.0),
    Named("aggregate_y_1") = aggregate_y_1,
    Named("aggregate_z_y_1") = aggregate_z_y_1,
    Named("aggregate_z_y_0") = aggregate_z_y_0);
  
  // Store calculation results
  return_List["individual_y_1"] = lnL_y_1;
  return_List["individual_z_y_1"] = lnL_z_y_1;
    
  aggregate_y_1 = sum(lnL_y_1);
  aggregate_z_y_1 = sum(lnL_z_y_1);
    
  return_List["aggregate_y_1"] = aggregate_y_1;
  return_List["aggregate_z_y_1"] = aggregate_z_y_1;

  return_List["individual_z_y_0"] = lnL_z_y_0;
  aggregate_z_y_0 = sum(lnL_z_y_0);
  return_List["aggregate_z_y_0"] = aggregate_z_y_0;
  
  return(return_List);
}

// Perform log-likelihood function estimation for selection model
double hpaSelectionLnLOptim(NumericVector x0, List hpaSelection_args)
{
  List return_List = hpaSelectionLnLOptim_List(x0, hpaSelection_args);
  
  double aggregate_y_1 = return_List["aggregate_y_1"];
  double aggregate_z_y_1 = return_List["aggregate_z_y_1"];
  double aggregate_z_y_0 = return_List["aggregate_z_y_0"];
  
  double return_aggregate = 0.0;

  return_aggregate += aggregate_z_y_0;
  return_aggregate += aggregate_y_1 + aggregate_z_y_1;

  return(return_aggregate);
}

// Perform log-likelihood function estimation for selection model
NumericVector hpaSelectionLnLOptim_ind(NumericVector x0, 
                                       List hpaSelection_args)
{
  List return_List = hpaSelectionLnLOptim_List(x0, hpaSelection_args);
  
  NumericVector individual_y_1 = return_List["individual_y_1"];
  NumericVector individual_z_y_1 = return_List["individual_z_y_1"];
  NumericVector individual_z_y_0 = return_List["individual_z_y_0"];
  
  int n_obs_0 = individual_z_y_0.size();
  int n_obs_1 = individual_y_1.size();
  int n_obs = n_obs_0 + n_obs_1;
  
  NumericVector return_individual = NumericVector(n_obs);
  
  return_individual[Range(0, n_obs_1 - 1)] = individual_y_1 + individual_z_y_1;
  return_individual[Range(n_obs_1, n_obs - 1)] = individual_z_y_0;
  
  return(return_individual);
}

// Gradient of likelihood function
List hpaSelectionLnLOptim_grad_List(NumericVector x0, List hpaSelection_args)
{
  // Get values from the hpaSelection_args list
  List ind_List = Rcpp::as<Rcpp::List>(hpaSelection_args["ind_List"]);
  arma::vec z_0 = hpaSelection_args["z_0"];
  arma::vec z_1 = hpaSelection_args["z_1"];
  arma::vec y_0 = hpaSelection_args["y_0"];
  arma::vec y_1 = hpaSelection_args["y_1"];
  arma::mat z_d_0 = hpaSelection_args["z_d_0"];
  arma::mat z_d_1 = hpaSelection_args["z_d_1"];
  arma::mat y_d_0 = hpaSelection_args["y_d_0"];
  arma::mat y_d_1 = hpaSelection_args["y_d_1"];
  NumericVector pol_degrees = hpaSelection_args["pol_degrees"];
  bool is_parallel = hpaSelection_args["is_parallel"];
  
  // Get values from the ind_List
  NumericVector pol_coefficients_ind = ind_List["pol_coefficients_ind"];
  int z_mean_ind = ind_List["z_mean_ind"];
  int y_mean_ind = ind_List["y_mean_ind"];
  int y_sd_ind = ind_List["y_sd_ind"];
  int z_sd_ind = ind_List["z_sd_ind"];
  NumericVector z_coef_ind = ind_List["z_coef_ind"];
  NumericVector y_coef_ind = ind_List["y_coef_ind"];
  
  // Assign estimated parameters values to corresponding vectors
  
    // polynomial coefficients and degrees
  NumericVector pol_coefficients = x0[pol_coefficients_ind];
  pol_coefficients.push_front(1); //add alpha(0...0)
  
    // mean
  double z_mean = x0[z_mean_ind];
  double y_mean = x0[y_mean_ind];
  NumericVector mean = {z_mean, y_mean};
  
    // sd
  double y_sd = x0[y_sd_ind];
  double z_sd = x0[z_sd_ind];
  NumericVector sd = {z_sd, y_sd};
  
  if (y_sd <= 0)
  {
    y_sd = std::sqrt(std::numeric_limits<double>::epsilon());
  }
  
  if (z_sd <= 0)
  {
    z_sd = std::sqrt(std::numeric_limits<double>::epsilon());
  }
  
  // coefficients for independent variables
  arma::vec y_coef = as<arma::vec>(x0[y_coef_ind]);
  NumericVector z_coef_temporal = x0[z_coef_ind];
  z_coef_temporal.push_front(1); //for fixed coefficient
  arma::vec z_coef = as<arma::vec>(z_coef_temporal);
  
  // get estimates for z*
  NumericVector z_h_1 = wrap(z_d_1 * z_coef);
  NumericVector z_h_0 = wrap(z_d_0 * z_coef);
  
  // get estimates for y and random errors
  arma::mat y_h_1 = y_d_1 * y_coef;
  NumericVector e_h_1 = wrap(y_1 - y_h_1);
  
  // concatenate e_h and z_h and prepare to insert into function
  NumericMatrix z_y_1 = NumericMatrix(z_h_1.size(), 2);
  NumericMatrix z_y_0 = NumericMatrix(z_h_0.size(), 2);
  
  z_y_1(_, 0) = -1.0 * z_h_1;
  z_y_1(_, 1) = e_h_1;
  z_y_0(_, 0) = -1.0 * z_h_0;

  // get number of observations with 0 and 1 values
  int n_obs_0 = z_h_0.size();
  int n_obs_1 = z_h_1.size();
  int n_obs = n_obs_0 + n_obs_1;
  
  // get the number of estimated parameters
  int n_param = x0.size();
    
  // Initialize vector to store gradient values
  NumericMatrix my_grad = NumericMatrix(n_obs_0 + n_obs_1, n_param);

  // infinity matricies
  NumericMatrix inf_y_vec_1 = NumericMatrix(n_obs_1, 2);
  std::fill(inf_y_vec_1.begin(), inf_y_vec_1.end(), R_PosInf);
  inf_y_vec_1(_, 1) = e_h_1; // in order to condition on e_h_1
  NumericMatrix neg_inf_vec_0 = NumericMatrix(n_obs_0, 2);
  std::fill(neg_inf_vec_0.begin(), neg_inf_vec_0.end(), R_NegInf);
  
  // Store the number of polynomial coefficients
  int pol_coefficients_n = pol_coefficients.size();
  
  // Gradient calculation
  
  NumericMatrix grad_y_1 = dhpaDiff(z_y_1,
                                    pol_coefficients, pol_degrees,
                                    LogicalVector{false, false}, 
                                    LogicalVector{true, false},
                                    mean, sd,
                                    "all",
                                    is_parallel, true, false);
  
  NumericMatrix grad_z_y_1 = ihpaDiff(z_y_1, inf_y_vec_1,
                                      pol_coefficients, pol_degrees,
                                      LogicalVector{false, true}, 
                                      LogicalVector{false, false},
                                      mean, sd,
                                      "all",
                                      is_parallel, true, false);
  
  NumericMatrix grad_z_y_0 = ihpaDiff(neg_inf_vec_0, z_y_0,
                                      pol_coefficients, pol_degrees,
                                      LogicalVector{false, false}, 
                                      LogicalVector{false, true},
                                      mean, sd,
                                      "all",
                                      is_parallel, true, false);
  
  // Store gradients respect to
  
    // polynomial coefficients, mean and sd
  for (int i = 0; i < (pol_coefficients_n + 3); i++) // for each parameter
  {
    NumericVector my_grad_tmp = NumericVector(n_obs);
    
    my_grad_tmp[Range(0, n_obs_1 - 1)] = grad_z_y_1(_, i + 1) + 
                                         grad_y_1(_, i + 1);
    my_grad_tmp[Range(n_obs_1, n_obs - 1)] = grad_z_y_0(_, i + 1);
    
    my_grad(_, i) = my_grad_tmp;
  }

    // selection equation regression coefficients
  NumericMatrix z_d_0_adj = wrap(-1.0 * z_d_0);
  NumericMatrix z_d_1_adj = wrap(-1.0 * z_d_1);
  
  int z_coef_n = z_coef_ind.size();
  
  for (int i = 0; i < z_coef_n; i++) // for each regressor
  {
    NumericVector my_grad_tmp = NumericVector(n_obs);

    my_grad_tmp[Range(0, n_obs_1 - 1)] = z_d_1_adj(_, i + 1) * 
                                         grad_z_y_1(_, pol_coefficients_n + 4);
    my_grad_tmp[Range(n_obs_1, n_obs - 1)] = z_d_0_adj(_, i + 1) * 
                                             grad_z_y_0(_, pol_coefficients_n + 6);

    my_grad(_, z_coef_ind[i]) = my_grad_tmp;
  }

    // outcome equation regression coefficients
  NumericMatrix y_d_1_adj = wrap(-1.0 * y_d_1);
    
  int y_coef_n = y_coef_ind.size();

  for (int i = 0; i < y_coef_n; i++) // for each regressor
  {
    NumericVector my_grad_tmp = NumericVector(n_obs);

    my_grad_tmp[Range(0, n_obs_1 - 1)] = y_d_1_adj(_, i) *
                                         (grad_y_1(_, pol_coefficients_n + 5) +
                                          grad_z_y_1(_, pol_coefficients_n + 7));

    my_grad(_, y_coef_ind[i]) = my_grad_tmp;
  }

  // Return the results
  List return_List = List::create(Named("aggregate") = colSums(my_grad),
                                  Named("individual") = my_grad);

  return(return_List);
}

// Gradient of likelihood function
NumericVector hpaSelectionLnLOptim_grad(NumericVector x0,
                                        List hpaSelection_args)
{
  List return_List = hpaSelectionLnLOptim_grad_List(x0, hpaSelection_args);
  
  NumericVector return_aggregate = return_List["aggregate"];
  
  return(return_aggregate);
}

// Perform log-likelihood function hessian estimation 
// for Phillips-Gallant-Nychka distribution at point
NumericMatrix hpaSelectionLnLOptim_hessian(NumericVector x0, 
                                           List hpaSelection_args)
{
  // Get values from the hpaSelection_args list
  List ind_List = Rcpp::as<Rcpp::List>(hpaSelection_args["ind_List"]);
  NumericVector pol_degrees = hpaSelection_args["pol_degrees"];
  
  // Get parameters number
  int n_param = x0.size();
  
  // Initialize vector to store Hessian values
  NumericMatrix my_hessian = NumericMatrix(n_param, n_param);
  
  // Get values from the list
  NumericVector pol_coefficients_ind = ind_List["pol_coefficients_ind"];
  int z_mean_ind = ind_List["z_mean_ind"];
  int y_mean_ind = ind_List["y_mean_ind"];
  int y_sd_ind = ind_List["y_sd_ind"];
  int z_sd_ind = ind_List["z_sd_ind"];
  NumericVector z_coef_ind = ind_List["z_coef_ind"];
  NumericVector y_coef_ind = ind_List["y_coef_ind"];
  
  // Assign estimated parameters values to corresponding vectors
  
  // polynomial coefficients and degrees
  NumericVector pol_coefficients = x0[pol_coefficients_ind];
  pol_coefficients.push_front(1); //add alpha(0...0)
  
  // mean
  double z_mean = x0[z_mean_ind];
  double y_mean = x0[y_mean_ind];
  NumericVector mean = {z_mean, y_mean};
  
  // sd
  double y_sd = x0[y_sd_ind];
  double z_sd = x0[z_sd_ind];
  NumericVector sd = {z_sd, y_sd};
  
  if (y_sd <= 0)
  {
    y_sd = std::sqrt(std::numeric_limits<double>::epsilon());
  }
  
  if (z_sd <= 0)
  {
    z_sd = std::sqrt(std::numeric_limits<double>::epsilon());
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
    g_plus = hpaSelectionLnLOptim_grad(x0_eps, hpaSelection_args);
    
    // Calculate g(x - eps)
    x0_eps[i] = x0[i] - eps[i];
    g_minus = hpaSelectionLnLOptim_grad(x0_eps, hpaSelection_args);
    
    // Store the results to hessian
    my_hessian(_, i) = (g_plus - g_minus) / (2 * eps[i]);
    
    // Set x0_eps value to default
    x0_eps[i] = x0[i];
  }
  
  // Make hessian to be symmetric
  for (int i = 0; i < n_param; i++) // for each parameter
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

// Gradient of likelihood function
NumericMatrix hpaSelectionLnLOptim_grad_ind(NumericVector x0, 
                                            List hpaSelection_args)
{
  List return_List = hpaSelectionLnLOptim_grad_List(x0, 
                                                    hpaSelection_args);
  
  NumericMatrix return_individual = return_List["individual"];
  
  return(return_individual);
}

//' Predict outcome and selection equation values from hpaSelection model
//' @description This function predicts outcome and selection equation 
//' values from hpaSelection model.
//' @param object Object of class "hpaSelection"
//' @param method string value indicating prediction method based on hermite 
//' polynomial approximation "HPA" or Newey method "Newey".
//' @template newdata_Template
//' @param is_cond logical; if \code{TRUE} (default) then conditional 
//' predictions will be estimated. Otherwise unconditional predictions 
//' will be returned.
//' @param is_outcome logical; if \code{TRUE} (default) then predictions 
//' for selection equation will be estimated using "HPA" method.
//' Otherwise selection equation predictions (probabilities) will be returned.
//' @details Note that Newey method can't predict conditional outcomes for 
//' zero selection equation value. Conditional probabilities for 
//' selection equation could be estimated only when dependent variable from 
//' outcome equation is observable.
//' @return This function returns the list which structure depends 
//' on \code{method}, \code{is_probit} and \code{is_outcome} values.
//' @export
// [[Rcpp::export]]
List predict_hpaSelection(List object, DataFrame newdata = R_NilValue, 
                          std::string method = "HPA", 
	bool is_cond = true, bool is_outcome = true)
{
	List model = object;

	// Add additional environments
	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function model_frame = stats_env["model.frame"];
	Rcpp::Function na_omit_R = stats_env["na.omit"];
	Rcpp::Function na_pass = stats_env["na.pass"];

	Rcpp::Environment base_env("package:base");
	Rcpp::Function as_data_frame = base_env["as.data.frame"];

	// Working with Data
	List Newey = model["Newey"];

	double z_mean = model["z_mean"];
	double y_mean = model["y_mean"];
	double z_sd = model["z_sd"];
	double y_sd = model["y_sd"];

	NumericVector pol_degrees = model["pol_degrees"];
	NumericVector pol_coefficients = model["pol_coefficients"];

		// Check whether new data frame has been supplied
	Rcpp::Formula selection = model["selection_formula"];
	Rcpp::Formula outcome = model["outcome_formula"];

	DataFrame data;

	if (newdata.size() == 0)
	{
		List data_List = model["data_List"];
		newdata = as_data_frame(data_List["dataframe"]);
	}

	data = newdata;

	// Extract data frame from formula
	DataFrame z_df = model_frame(Rcpp::_["formula"] = selection, 
                               Rcpp::_["data"] = newdata, 
                               Rcpp::_["na.action"] = na_pass);
	DataFrame y_df = model_frame(Rcpp::_["formula"] = outcome, 
                               Rcpp::_["data"] = newdata, 
                               Rcpp::_["na.action"] = na_pass);

	CharacterVector z_df_names = z_df.names();
	CharacterVector y_df_names = y_df.names();

	int z_df_n = z_df.size();
	int y_df_n = y_df.size();

	// Extract dependent variables
	NumericVector z = z_df[0];       // it is reference
	NumericVector y = y_df[0];       // it is reference

	int n = z.size();

	// Extract independent variable
	NumericMatrix z_d(n, z_df_n - 1); // -1 because there is no constant term
	NumericMatrix y_d(n, y_df_n - 1); // -1 because there is no constant term

	int z_d_col = z_d.ncol();
	int y_d_col = y_d.ncol();

	for (int i = 0; i < z_d_col; i++)
	{
		z_d(_, i) = NumericVector(z_df[i + 1]);
	}
	for (int i = 0; i < y_d_col; i++)
	{
		y_d(_, i) = NumericVector(y_df[i + 1]);
	}

	// calculate latent variables values
	NumericVector z_latent = NumericVector(n);
	NumericVector z_coef;
	
	if (method == "HPA")
	{
		z_coef = model["z_coef"];
	}

	if (method == "Newey")
	{
		z_coef = Newey["z_coef"];
	}

	for (int i = 0; i < z_d_col; i++)
	{
		z_latent = z_latent + z_d(_, i) * z_coef[i];
	}

	// calculate unconditional y values without constant
	NumericVector y_uncond = NumericVector(n);
	
	NumericVector y_coef;

	if (method == "HPA")
	{
		y_coef = model["y_coef"];
	}

	if (method == "Newey")
	{
		y_coef = Newey["y_coef"];
	}

	for (int i = 0; i < y_d_col; i++)
	{
		y_uncond = y_uncond + y_d(_, i) * y_coef[i];
	}

	// calculate conditional expectations
	List results;

	NumericVector y_cond = NumericVector(n);

	if (!is_outcome)
	{
		NumericMatrix e_mat = NumericMatrix(n, 1);
		e_mat(_, 0) = y - y_uncond;
		NumericMatrix z_y = NumericMatrix(n, 2);
		z_y(_, 0) = z_latent * (-1);
		z_y(_, 1) = e_mat;
		if (is_cond)
		{
			NumericVector z_prob = 1 - phpa(z_y, pol_coefficients,
				pol_degrees,
				LogicalVector::create(false, true), LogicalVector::create(false, false),
				NumericVector::create(z_mean, y_mean), 
				                      NumericVector::create(z_sd, y_sd), 
				false, false, false);

			return(List::create(Named("prob") = z_prob));
		} else {
			NumericVector z_prob = 1 - phpa(z_y, pol_coefficients,
				pol_degrees,
				LogicalVector::create(false, false), LogicalVector::create(false, true),
				NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, 
                                                                     y_sd), 
				false, false, false);
		  
			return(List::create(Named("prob") = z_prob));
		}
	}

	if (method == "HPA")
	{
		if (!is_cond)
		{
			NumericVector errors_exp_vec = ehpa(NumericMatrix(1, 1),
				pol_coefficients, pol_degrees,
				LogicalVector::create(false, false),
				LogicalVector::create(false, false),
				NumericVector::create(z_mean, y_mean), 
				NumericVector::create(y_mean, y_sd),
				NumericVector::create(0, 1), false, false);
			double errors_exp = errors_exp_vec[0];

			return(List::create(Named("y") = y_uncond + errors_exp));
		}

		NumericVector NegInfVec = NumericVector(n);
		NumericVector PosInfVec = NumericVector(n);

		std::fill(NegInfVec.begin(), NegInfVec.end(), R_NegInf);
		std::fill(PosInfVec.begin(), PosInfVec.end(), R_PosInf);

		// for 1 outcome
		NumericMatrix lower_1 = NumericMatrix(n, 2);
		lower_1(_, 0) = (-1.0) * z_latent;
		lower_1(_, 1) = NegInfVec;

		NumericMatrix upper_1 = NumericMatrix(n, 2);
		upper_1(_, 0) = PosInfVec;
		upper_1(_, 1) = PosInfVec;

		NumericVector e_tr_1 = etrhpa(lower_1, upper_1,
			pol_coefficients, pol_degrees,
			NumericVector::create(z_mean, y_mean),
			NumericVector::create(z_sd, y_sd),
			NumericVector::create(0, 1), false, false);

		NumericVector y_cond_1 = y_uncond + e_tr_1;

		// for 0 outcome
		NumericMatrix lower_0 = NumericMatrix(n, 2);
		lower_0(_, 0) = NegInfVec;
		lower_0(_, 1) = NegInfVec;

		NumericMatrix upper_0 = NumericMatrix(n, 2);
		upper_0(_, 0) = (-1) * z_latent;
		upper_0(_, 1) = PosInfVec;

		NumericVector e_tr_0 = etrhpa(lower_0, upper_0,
			pol_coefficients, pol_degrees,
			NumericVector::create(z_mean, y_mean),
			NumericVector::create(z_sd, y_sd),
			NumericVector::create(0, 1), false, false);

		NumericVector y_cond_0 = y_uncond + e_tr_0;

		// aggregate result
		y_cond[z == 1] = y_cond_1[z == 1];
		y_cond[z == 0] = y_cond_0[z == 0];

		results = List::create(Named("y") = y_cond,
			Named("y_1") = y_cond_1,
			Named("y_0") = y_cond_0);
	}

	if (method == "Newey")
	{

		if (!is_cond)
		{
			return(List::create(Named("y") = y_uncond));
		}

		y_cond = y_uncond;

		double z_exp = Newey["selection_exp"];
		double z_var = Newey["selection_var"];

		int pol_elements = Newey["pol_elements"];

		NumericVector y_coef_mills = Newey["inv_mills_coef"];
		z_latent = (z_latent + z_exp) / sqrt(z_var);
		NumericVector z_mills = dnorm(z_latent) / pnorm(z_latent);

		for (int i = 0; i <= pol_elements; i++)
		{
			y_cond = y_cond + pow(z_mills, i) * y_coef_mills[i];
		}

		results = List::create(Named("y_1") = y_cond);
	}

	return(results);
}

//' Summarizing hpaSelection Fits
//' @description This function summarizing hpaSelection Fits.
//' @param object Object of class "hpaSelection".
//' @return This function returns the same list as 
//' \code{\link[hpa]{hpaSelection}} function changing its class 
//' to "summary.hpaSelection".
//' @export
// [[Rcpp::export]]
List summary_hpaSelection(List object) 
{

	List return_result = clone(object); //in order to preserve model class

	return_result.attr("class") = "summary.hpaSelection";

	return(return_result);
}

//' Summary for hpaSelection output
//' @param x Object of class "hpaSelection"
//' @export
// [[Rcpp::export]]
void print_summary_hpaSelection(List x) {

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

	// extract list of indexes
	List ind_List = model["ind_List"];
	NumericVector pol_coefficients_ind = ind_List["pol_coefficients_ind"];
	NumericVector z_mean_ind = ind_List["z_mean_ind"];
	NumericVector z_coef_ind = ind_List["z_coef_ind"];
	NumericVector y_mean_ind = ind_List["y_mean_ind"];
	NumericVector y_sd_ind = ind_List["y_sd_ind"];
	NumericVector y_coef_ind = ind_List["y_coef_ind"];

	// other stuff

	NumericVector x1 = model["x1"];

	List data_List = model["data_List"];

	DataFrame data = as_data_frame(data_List["dataframe"]);
	DataFrame data_z = as_data_frame(data_List["data_z"]);
	StringVector data_names_z = data_z.names();
	NumericVector z = data_z[0];
	int n_censored = sum(z);

	NumericVector p_values = results(_, 3);

	StringVector stars = starVector(p_values);

	NumericVector z_coef = model["z_coef"];
	NumericVector y_coef = model["y_coef"];

	StringVector results_rownames = rownames(results);
	StringVector results_colnames = colnames(results);
	
	int n_obs = model["n_obs"];
	int df = x1.size();

	double lnL = model["log-likelihood"];
	double AIC = 2 * (df - lnL);
	double n_obs_double = n_obs * 1.0;
	double BIC = log(n_obs_double) * x1.size() - 2 * lnL;

	std::string lnL_string = "Log-Likelihood: " + std::to_string(lnL) + "\n";
	std::string AIC_string = "AIC: " + std::to_string(AIC) + "\n";
	std::string BIC_string = "BIC: " + std::to_string(BIC) + "\n";
	std::string n_obs_string = "Observations: " + std::to_string(n_obs) + 
	                           " (" + std::to_string(n_censored) + 
	                           " observed)"+ "\n";
	std::string df_string = std::to_string(df) + 
	                        " free parameters (df = " + 
	                        std::to_string(n_obs - df) + ")" + "\n";

	cat_R("--------------------------------------------------------------\n");

	cat_R("Semi-nonparametric selection model estimation\n");

	cat_R(lnL_string.c_str());
	cat_R(AIC_string.c_str());
	cat_R(BIC_string.c_str());
	cat_R(n_obs_string.c_str());
	cat_R(df_string.c_str());
	
	cat_R("---\n");

	cat_R("Selection equation coefficients:\n");
	int z_coef_first = z_coef_ind[0];
	int z_coef_last = z_coef_ind[z_coef_ind.size() - 1];
	NumericMatrix z_coef_results = results(Range(z_coef_first, z_coef_last), _);
	rownames(z_coef_results) = results_rownames[z_coef_ind];
	colnames(z_coef_results) = results_colnames;
	print(as_table(cbind(z_coef_results, stars[z_coef_ind])));

	cat_R("---\n");

	cat_R("Outcome equation coefficients:\n");
	int y_coef_first = y_coef_ind[0];
	int y_coef_last = y_coef_ind[y_coef_ind.size() - 1];
	NumericMatrix y_coef_results = results(Range(y_coef_first, y_coef_last), _);
	rownames(y_coef_results) = results_rownames[y_coef_ind];
	colnames(y_coef_results) = results_colnames;
	print_R(as_table(cbind(y_coef_results, stars[y_coef_ind])));

	cat_R("---\n");

	cat_R("Distribution parameters:\n");
	int distr_first = 0;
	int distr_last = df - z_coef_ind.size() - y_coef_ind.size() - 1;
	NumericMatrix distr_results = results(Range(distr_first, distr_last), _);
	StringVector distr_rownames = results_rownames[Range(distr_first, 
                                                       distr_last)];
	rownames(distr_results) = distr_rownames;
	colnames(distr_results) = results_colnames;
	print_R(as_table(cbind(distr_results, stars[Range(distr_first, 
                                                    distr_last)])));

	cat_R("---\n");

	cat_R("Selection equation fixed coefficients:\n");
	String new_str_names = data_names_z(1);
	std::string new_str_names_str = new_str_names;
	std::string new_str = new_str_names_str + " = 1" + "\n";
	cat_R(new_str.c_str());

	cat_R("---\n");

	cat_R("Fixed Distribution Parameters:\n");
	cat_R("a_0 = 1\n");

	cat_R("---\n");

	List re_moments = model["re_moments"];
	double rho = re_moments["rho"];
	double rho_std = re_moments["rho_std"];
	NumericVector rho_F_z_stat = pnorm(NumericVector::create(rho / rho_std));
	double rho_p_value = std::min(rho_F_z_stat[0], 1 - rho_F_z_stat[0]);
	std::string new_str_rho = "Correlation between random errors (rho) is: \n Estimate: " +
	                          std::to_string(rho) +
	                          "\n Standard Error (delta method): " + 
	                          std::to_string(rho_std) +
	                          "\n p-value (H0: rho = 0): " + 
	                          std::to_string(rho_p_value) + "\n";
	cat_R(new_str_rho.c_str());

	cat_R("---\n");

	cat_R("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n");

	cat_R("--------------------------------------------------------------\n");
}

//' Plot hpaSelection random errors approximated density
//' @param x Object of class "hpaSelection"
//' @param is_outcome logical; if TRUE then function plots the graph for 
//' outcome equation random errors. 
//' Otherwise plot for selection equation random errors will be plotted.
//' @return This function returns the list containing random error 
//' expected value \code{errors_exp}
//' and variance \code{errors_var} estimates for selection 
//' (if \code{is_outcome = TRUE}) or outcome (if \code{is_outcome = FALSE}) 
//' equation.
//' @export
// [[Rcpp::export]]
List plot_hpaSelection(List x, bool is_outcome = true) {

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

	// load data from the model

	double z_mean = model["z_mean"];
	double y_mean = model["y_mean"];
	double z_sd = model["z_sd"];
	double y_sd = model["y_sd"];

	NumericVector mean = NumericVector::create(z_mean, y_mean);
	NumericVector sd = NumericVector::create(z_sd, y_sd);

	NumericVector pol_degrees = model["pol_degrees"];
	NumericVector pol_coefficients = model["pol_coefficients"];

	// get random errors expectation and variance
	int eq_ind = is_outcome;

		NumericVector errors_exp_vec = ehpa(NumericMatrix(1, 1), 
			pol_coefficients, pol_degrees,
			LogicalVector::create(false, false), 
			LogicalVector::create(false, false),
			mean, sd,
			NumericVector::create(1 - eq_ind, eq_ind), false, false);
		double errors_exp = errors_exp_vec[0];

		NumericVector errors_exp_2_vec = ehpa(NumericMatrix(1, 1), 
			pol_coefficients, pol_degrees,
			LogicalVector::create(false, false), 
			LogicalVector::create(false, false),
			mean, sd,
			NumericVector::create(2 * (1 - eq_ind), 2 * eq_ind), false, false);
		double errors_exp_2 = errors_exp_2_vec[0];

		double errors_var = errors_exp_2 - (errors_exp * errors_exp);

	// adjust precision

	double plot_min = errors_exp - 3 * sqrt(errors_var);
	double plot_max = errors_exp + 3 * sqrt(errors_var);

	int n = 10000;

	double precise = (plot_max - plot_min) / n;

	NumericMatrix x_matr = NumericMatrix(n, 2);
	x_matr(0, 0) = plot_min;
	x_matr(0, 1) = plot_min;

	for (int i = 1; i < n; i++)
	{
		x_matr(i, eq_ind) = x_matr(i - 1, eq_ind) + precise;
	}

	// calculate densities

	NumericVector den = dhpa(x_matr,
		pol_coefficients, pol_degrees,
		LogicalVector::create(false, false),
		LogicalVector::create(is_outcome, !is_outcome),
		mean, sd, false, false, false);

	double den_min = min(den);
	double den_max = max(den);

	NumericVector x_vec = x_matr(_, eq_ind);

	plot_R(Rcpp::_["x"] = x_vec, Rcpp::_["y"] = den,
		Rcpp::_["xlim"] = NumericVector::create(plot_min, plot_max),
		Rcpp::_["xaxp"] = NumericVector::create(plot_min, plot_max, 5),
		Rcpp::_["yaxp"] = NumericVector::create(den_min, den_max, 5),
		Rcpp::_["type"] = "l",
		Rcpp::_["lwd"] = 3,
		Rcpp::_["main"] = "Random Errors Density Approximation Plot",
		Rcpp::_["xlab"] = "value",
		Rcpp::_["ylab"] = "density");

	List moments = List::create(Named("errors_exp") = errors_exp,
		Named("errors_var") = errors_var);

	return(moments);
}

//' Calculates log-likelihood for "hpaSelection" object
//' @description This function calculates log-likelihood for 
//' "hpaSelection" object
//' @param object Object of class "hpaSelection"
//' @export
// [[Rcpp::export]]
double logLik_hpaSelection(List object) 
{

	double lnL = object["log-likelihood"];

	return(lnL);
}
