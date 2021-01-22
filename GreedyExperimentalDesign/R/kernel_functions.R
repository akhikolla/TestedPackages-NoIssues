#a set of the valid kernel types for the method build_kpca_object
VALID_KERNEL_TYPES = c("vanilla", "rbf", "poly", "tanh", "bessel", "laplace", "anova", "spline")



#' Gram Matrix Computation
#' 
#' Computes the Gram Matrix for a user-specified kernel using the library \code{kernlab}. Note that 
#' this function automatically standardizes the columns of the data entered.
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design.
#' @param kernel_type 		One of the following: "vanilla", "rbf", "poly", "tanh", "bessel", "laplace", 
#' 							"anova" or "spline".
#' @param params			A vector of numeric parameters. Each \code{kernel_type} has different numbers of
#' 							parameters required. For more information see documentation for the \code{kernlab} library.
#' @return					The \code{n x n} gram matrix for the given kernel on the given data. 
#' 
#' @author Adam Kapelner
#' @export
compute_gram_matrix = function(X, kernel_type, params = c()){
	#ensure that the kernel_type is valid
	if (!(kernel_type %in% VALID_KERNEL_TYPES)){
		stop("The kernel type must be one of the following:", paste(VALID_KERNEL_TYPES, collapse = ", "))
	}
	
	#now build the the correct kernel based on the kernel type. This is a giant switch statement
	if (kernel_type == "vanilla"){
		if (length(params) != 0){
			stop("vanilla kernels do not take parameters")
		}
		kernel = kernlab::vanilladot()
	} else if (kernel_type == "rbf"){
		if (length(params) != 1){
			stop("gaussian radial basis function kernels require one parameter: gamma")
		}
		kernel = kernlab::rbfdot(params[1])
	} else if (kernel_type == "poly"){
		if (length(params) != 3){
			stop("polynomial kernels require three parameters: degree, scale and offset in that order")
		}
		kernel = kernlab::polydot(params[1], params[2])
	} else if (kernel_type == "tanh"){
		if (length(params) != 2){
			stop("hyperbolic tangent kernels require two parameters: scale and offset in that order")
		}
		kernel = kernlab::tanhdot(params[1], params[2])
	} else if (kernel_type == "bessel"){
		if (length(params) != 3){
			stop("bessel kernels require three parameters: gamma, order and degree in that order")
		}
		kernel = kernlab::besseldot(params[1], params[2], params[3])
	} else if (kernel_type == "laplace"){
		if (length(params) != 1){
			stop("laplace kernels require one parameter: gamma")
		}
		kernel = kernlab::laplacedot(params[1])
	} else if (kernel_type == "anova"){
		if (length(params) != 2){
			stop("anova kernels require two parameters: gamma and degree in that order")
		}
		kernel = kernlab::anovadot(params[1], params[2])
	} else if (kernel_type == "spline"){
		if (length(params) != 0){
			stop("spline kernels do not take parameters")
		}
		kernel = kernlab::splinedot()
	}
	
	#compute the Gram matrix as type matrix
	kernlab::kernelMatrix(kernel, standardize_data_matrix(X))
}


#private
setGramMatrix = function(java_obj, Kgram){
	for (i in 1 : nrow(Kgram)){	
		.jcall(java_obj, "V", "setKgramRow", as.integer(i - 1), Kgram[i, , drop = FALSE]) #java indexes from 0...n-1
	}
}