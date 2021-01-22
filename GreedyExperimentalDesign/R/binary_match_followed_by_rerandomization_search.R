#' Begin a Search for Binary Matching Followed by Rerandomization
#' 
#' This method creates an object of type binary_then_rerandomization_experimental_design and will find optimal matched pairs which
#' are then rerandomized in order to further minimize a balance metric. You can then
#' use the function \code{resultsBinaryMatchThenRerandomizationSearch} to obtain the randomized allocation vectors. For one column
#' in X, the matching just sorts the values to find the pairs trivially.
#' 
#' @param X						The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 								(one for each measurement on the subject). This is the design matrix you wish 
#' 								to search for a more optimal design.
#' @param compute_dist_matrix	The function that computes the distance matrix between every two observations in \code{X}, 
#' 								its only argument. The default is \code{NULL} signifying euclidean squared distance optimized in C++.
#' @param ...					Arguments passed to \code{initGreedyExperimentalDesignObject}. It is recommended to set
#' 								\code{max_designs} otherwise it will default to 10,000.
#' @return						An object of type \code{binary_experimental_design} which can be further operated upon.
#' 
#' @author Adam Kapelner
#' @export
binaryMatchFollowedByRerandomizationDesignSearch = function(X, compute_dist_matrix = NULL, ...){
	n = nrow(X)
	p = ncol(X)
	
	if (n %% 4 != 0){
		stop("Design matrix must have number of rows divisible by four for this type of design.")
	}
	binary_match_design = binaryMatchExperimentalDesignSearch(X, compute_dist_matrix)
	#now we create a reduced matrix X by diffing the pairs
	Xdiffs = matrix(NA, nrow = nrow(X) / 2, ncol = ncol(X))
	for (i in 1 : (nrow(X) / 2)){		
		Xdiffs[i, ] = X[binary_match_design$indices_pairs[i, 1], ] - X[binary_match_design$indices_pairs[i, 2], ]
	}
	
	binary_then_rerandomization_experimental_design = list()
	binary_then_rerandomization_experimental_design$X = X
	binary_then_rerandomization_experimental_design$n = n
	binary_then_rerandomization_experimental_design$p = p
	binary_then_rerandomization_experimental_design$binary_match_design = binary_match_design
	binary_then_rerandomization_experimental_design$rerandomization_design = initRerandomizationExperimentalDesignObject(Xdiffs, ...)
	class(binary_then_rerandomization_experimental_design) = "binary_then_rerandomization_experimental_design"
	binary_then_rerandomization_experimental_design
}

#' Returns unique allocation vectors that are binary matched
#' 
#' @param obj 				The \code{binary_then_greedy_experimental_design} object where the pairs are computed.
#' @param num_vectors		How many random allocation vectors you wish to return. The default is \code{NULL} indicating you want all of them.
#' @param compute_obj_vals	Should we compute all the objective values for each allocation? Default is \code{FALSE}.
#' @param form				Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' 
#' @author Adam Kapelner
#' @export
resultsBinaryMatchThenRerandomizationSearch = function(obj, num_vectors = NULL, compute_obj_vals = FALSE, form = "zero_one"){
	assertClass(obj, "binary_then_rerandomization_experimental_design")
	assertCount(num_vectors, positive = TRUE, null.ok = TRUE)
	if (is.null(num_vectors)){
		num_vectors = obj$rerandomization_design$max_designs
	}
	num_vectors_completed = rerandomizationSearchCurrentProgress(obj$rerandomization_design)
	if (num_vectors > num_vectors_completed){
		warning("You requested ", num_vectors, " but only ", num_vectors_completed, " are available.")
		num_vectors = num_vectors_completed
	}
	assertLogical(compute_obj_vals)
	assertChoice(form, c("zero_one", "pos_one_min_one"))
	
	rerand_res = resultsRerandomizationSearch(obj$rerandomization_design, include_assignments = TRUE, "zero_one")
	#the allocation vectors returned here have entries = 1 if the pair is left unswitched and = 0 if we should switch the pair
	
	indicTs = matrix(NA, nrow = num_vectors, ncol = obj$n)
	for (r in 1 : num_vectors){
		#first we copy the binary indices starting point
		pair_matrix_copy = obj$binary_match_design$indices_pairs
		#now we pull out a w_diff
		w_diff = rerand_res$ending_indicTs[r, ]
		
		#split the pair matrix based on the greedy vector
		pair_matrix_T_is_first = pair_matrix_copy[w_diff == 1, ] #"zero_one" form was forced above
		pair_matrix_C_is_first = pair_matrix_copy[w_diff == 0, ] #"zero_one" form was forced above
		#now set all the entries as T if T is first and it's in the first column and 0 if second column
		indicTs[r, pair_matrix_T_is_first[, 1]] = 1
		indicTs[r, pair_matrix_T_is_first[, 2]] = 0
		#vice versa if C is first
		indicTs[r, pair_matrix_C_is_first[, 1]] = 0
		indicTs[r, pair_matrix_C_is_first[, 2]] = 1		
	}	
	
	obj_vals = NULL
	if (compute_obj_vals){
		if (obj$rerandomization_design$objective == "mahal_dist"){
			SinvX = solve(var(obj$X))
			obj_vals = apply(indicTs, 1, FUN = function(w){compute_objective_val(obj$X, w, objective = "mahal_dist", SinvX)})
		} else {
			obj_vals = apply(indicTs, 1, FUN = function(w){compute_objective_val(obj$X, w, objective = "abs_sum_diff")})	
		}	
	}

	if (form == "pos_one_min_one"){
		indicTs = (indicTs - 0.5) * 2
	}
	list(
		obj_vals = obj_vals,
		indicTs = indicTs,
		form = form
	)
}

#' Prints a summary of a \code{binary_then_rerandomization_experimental_design} object
#' 
#' @param x			The \code{binary_then_rerandomization_experimental_design} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print binary_then_rerandomization_experimental_design
#' @export
print.binary_then_rerandomization_experimental_design = function(x, ...){
	print.binary_then_rerandomization_experimental_design(x$rerandomization_design)
}

#' Prints a summary of a \code{binary_then_rerandomization_experimental_design} object
#' 
#' @param object		The \code{binary_then_rerandomization_experimental_design} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary binary_then_rerandomization_experimental_design
#' @export
summary.binary_then_rerandomization_experimental_design = function(object, ...){
	print(object, ...)
}
