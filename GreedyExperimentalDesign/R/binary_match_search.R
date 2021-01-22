#' Begin a Search for Binary Matching Designs
#' 
#' This method creates an object of type binary_experimental_design and will find pairs. You can then
#' use the function \code{resultsBinaryMatchSearch} to create randomized allocation vectors. For one column
#' in X, we just sort to find the pairs trivially.
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design.
#' @param compute_dist_matrix	The function that computes the distance matrix between every two observations in \code{X}, 
#' 								its only argument. The default is \code{NULL} signifying euclidean squared distance optimized in C++.
#' @return					An object of type \code{binary_experimental_design} which can be further operated upon.
#' 
#' @author Adam Kapelner
#' @export
binaryMatchExperimentalDesignSearch = function(X, compute_dist_matrix = NULL){
	assertClass(X, "matrix")
	assertClass(compute_dist_matrix, "function", null.ok = TRUE)
	n = nrow(X)
	p = ncol(X)
	if (n %% 2 != 0){
		stop("Design matrix must have even rows to have equal treatments and controls")
	}
	
	if (is.null(compute_dist_matrix) & p == 1){
		#we don't need to do anything except order them up
		indices_pairs = matrix(order(X[, 1]), ncol = 2, byrow = TRUE)
	} else {
		if (is.null(compute_dist_matrix)) {	#default is C++-optimized sqd euclidean distance function		
			D = compute_distance_matrix_cpp(X)
		} else {
			D = compute_dist_matrix(X)
		}
		#ensure diagonal is infinity
		diag(D) = .Machine$double.xmax
		
		#get the matching solution using the heuristic
		indices_pairs = as.matrix(nbpMatching::nonbimatch(nbpMatching::distancematrix(D))$matches[, c("Group1.Row", "Group2.Row")])
		
		for (i in 1 : n){
			indices_pairs[i, ] = sort(indices_pairs[i, ])
		}	
		indices_pairs = unique(indices_pairs)
	}

	#now return information as an object (just a list)
	binary_experimental_design = list()
	binary_experimental_design$X = X
	binary_experimental_design$n = n
	binary_experimental_design$p = p
	binary_experimental_design$compute_dist_matrix = compute_dist_matrix
	binary_experimental_design$D = D
	binary_experimental_design$indices_pairs = indices_pairs
	class(binary_experimental_design) = "binary_experimental_design"
	binary_experimental_design
}

#' Returns unique allocation vectors that are binary matched
#' 
#' @param obj 				The \code{binary_experimental_design} object where the pairs are computed.
#' @param num_vectors		How many random allocation vectors you wish to return. The default is 1000.
#' @param objective			Should we compute all the objective values for each allocation? Default is \code{NULL} for "no".
#' 							If non-null, it needs to either be "mahal_dist" or "abs_sum_diff".
#' @param form				Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' 
#' @author Adam Kapelner
#' @export
resultsBinaryMatchSearch = function(obj, num_vectors = 1000, objective = NULL, form = "zero_one"){
	assertClass(obj, "binary_experimental_design")
	assertCount(num_vectors, positive = TRUE)
	assertChoice(form, c("zero_one", "pos_one_min_one"))
	if (!is.null(objective)){
		verify_objective_function(objective)
	}	
	
	#now that we have the pairs, we can randomize for as many vectors as we wish
	n = obj$n
	if (2^(n / 2) < num_vectors){
		stop(paste("The total number of unique vectors is", 2^(n / 2), "which is less than the", num_vectors, "you requested."))
	}

	minus_half_plus_half = c(-.5, .5)
	
	indicTs = matrix(NA, nrow = 0, ncol = n)
	batch_size = ceiling(num_vectors / 16) * 4 #needs to be divisible by 4 because when divided by two, it must be even
	repeat {
		indicTs_batch = matrix(NA, nrow = batch_size, ncol = n)
		one_minus_one = matrix(
				sample(c(
					rep(1, batch_size / 2 * n / 2), 
					rep(-1, batch_size / 2 * n / 2)
				)), 
				nrow = batch_size)
		for (r in 1 : batch_size){
			for (i in 1 : (n / 2)){
				indicTs_batch[r, obj$indices_pairs[i, ]] = minus_half_plus_half * one_minus_one[r, i] + 0.5
			}
		}
		indicTs = rbind(indicTs, indicTs_batch)
		indicTs = unique(indicTs)
		if (nrow(indicTs) >= num_vectors){
			indicTs = indicTs[1 : num_vectors, , drop = FALSE]
			break
		}
	}
	
	obj_vals = NULL
	if (!is.null(objective)){
		verify_objective_function(objective)
		if (objective == "mahal_dist"){
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
		obj_vals_unordered = obj_vals,
		indicTs = indicTs,
		form = form
	)
}

#' Prints a summary of a \code{binary_experimental_design} object
#' 
#' @param x			The \code{binary_experimental_design} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print binary_experimental_design
#' @export
print.binary_experimental_design = function(x, ...){
	cat("The pairs have been computed. Now use the resultsBinaryMatchSearch to make allocations.\n")
}

#' Prints a summary of a \code{binary_experimental_design} object
#' 
#' @param object		The \code{binary_experimental_design} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary binary_experimental_design
#' @export
summary.binary_experimental_design = function(object, ...){
	print(object, ...)
}
