#' Starts the parallelized greedy design search. 
#' 
#' Once begun, this function cannot be run again.
#' 
#' @param obj 		The \code{experimental_design} object that will be running the search
#' 
#' @author Adam Kapelner
#' @export
startSearch = function(obj){
	if (.jcall(obj$java_obj, "Z", "began")){
		stop("Search Already begun.")
	}
	.jcall(obj$java_obj, "V", "beginSearch")
}

#' Stops the parallelized greedy design search. 
#' 
#' Once stopped, it cannot be restarted.
#' 
#' @param obj 		The \code{experimental_design} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @export
stopSearch = function(obj){
	.jcall(obj$java_obj, "V", "stopSearch")
}

#' Generates a design matrix with standardized predictors. 
#' 
#' This function is useful for debugging.
#' 
#' @param n					Number of rows in the design matrix 
#' @param p 				Number of columns in the design matrix
#' @param covariate_gen		The function to use to draw the covariate realizations (assumed to be iid).
#' 							This defaults to \code{rnorm} for $N(0,1)$ draws.
#' @param ...				Optional arguments to be passed to the \code{covariate_dist} function.
#' @return 					THe design matrix
#' 
#' @author Adam Kapelner
#' @export
generate_stdzied_design_matrix = function(n = 50, p = 1, covariate_gen = rnorm, ...){
	X = matrix(covariate_gen(n * p, ...), nrow = n, ncol = p)
	#now standardize the matrix to make things easier later
	apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})	
}


#' Returns the amount of time elapsed
#' 
#' @param obj 		The \code{experimental_design} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @export
searchTimeElapsed = function(obj){
	.jcall(obj$java_obj, "I", "timeElapsedInSeconds")
}

#' Computes Objective Value From Allocation Vector
#' 
#' Returns the objective value given a design vector as well an an objective function.
#' This is sometimes duplicated in Java. However, within Java, tricks are played to make
#' optimization go faster so Java's objective values may not always be the same as the true
#' objective function (e.g. logs or constants dropped).
#' 
#' @param X 		 	The n x p design matrix
#' @param indic_T		The n-length binary allocation vector
#' @param objective		The objective function to use. Default is \code{abs_sum_diff} and the other option is 
#' 						\code{mahal_dist}.
#' @param inv_cov_X		Optional: the inverse sample variance covariance matrix. Use this
#' 						argument if you will be doing many calculations since passing this
#' 						in will cache this data.
#' 
#' @author Adam Kapelner
#' @export
compute_objective_val = function(X, indic_T, objective = "abs_sum_diff", inv_cov_X = NULL){
	if (!isTRUE(all.equal(sort(unique(indic_T)), c(0, 1)))){
		stop("indic_T must be binary")
	}
	X_T = X[indic_T == 1, , drop = FALSE] #coerce as matrix in order to matrix multiply later
	X_C = X[indic_T == 0, , drop = FALSE] #coerce as matrix in order to matrix multiply later
	X_T_bar = colMeans(X_T)
	X_C_bar = colMeans(X_C)	
	
	if (objective == "abs_sum_diff"){
		s_j_s = apply(X, 2, sd)
		sum(abs((X_T_bar - X_C_bar) / s_j_s))
	} else if (objective == "mahal_dist"){
		#saves computation to pass it in if you're doing a lot of them in a row
		if (is.null(inv_cov_X)){
			inv_cov_X = solve(var(X))
		}	
		X_T_bar_minus_X_C_bar = as.matrix(X_T_bar - X_C_bar) #need to matricize for next computation
		as.numeric(t(X_T_bar_minus_X_C_bar) %*% inv_cov_X %*% X_T_bar_minus_X_C_bar)
	} else {
		stop("objective invalid.")
	}
}

#' Standardizes the columns of a data matrix.
#' 
#' @param X 		 	The n x p design matrix
#' @return				The n x p design matrix with columns standardized
#' 
#' @author Adam Kapelner
#' @export
standardize_data_matrix = function(X){
	apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
}

#private
verify_objective_function = function(objective, Kgram = NULL, n = NULL){
	if (objective != "mahal_dist" && objective != "abs_sum_diff" && objective != "kernel"){
		stop("Objective function must be one of the following:\n  mahal_dist\n  abs_sum_diff\n  kernel\n\n")
	}
	if (objective == "kernel"){
		if (is.null(Kgram) || is.null(n)){
			stop("You must specify a gram matrix \"Kgram\" and \"n\".\n")
		}
		if (class(Kgram) != "kernelMatrix" && class(Kgram) != "matrix"){
			stop("The gram matrix must be type kernelMatrix or type matrix.\n")
		}
		if (!all.equal(dim(Kgram), c(n, n))){
			stop("The gram matrix must have dimension n x n.\n")
		}
	}
	if (!is.null(Kgram) && objective != "kernel"){
		stop("If you specify a gram matrix, you must specify the \"kernel\" objective.\n")
	}
}