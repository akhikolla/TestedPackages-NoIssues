#' Begin a Search for the Optimal Solution
#' 
#' This method creates an object of type optimal_experimental_design and will immediately initiate
#' a search through $1_{T}$ space. Since this search takes exponential time, for most machines, 
#' this method is futile beyond 28 samples. You've been warned!
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design.
#' @param objective			The objective function to use when searching design space. This is a string
#' 							with valid values "\code{mahal_dist}" (the default), "\code{abs_sum_diff}" or "\code{kernel}".
#' @param Kgram				If the \code{objective = kernel}, this argument is required to be an \code{n x n} matrix whose
#' 							entries are the evaluation of the kernel function between subject i and subject j. Default is \code{NULL}.
#' @param wait				Should the \code{R} terminal hang until all \code{max_designs} vectors are found? The 
#' 							deafult is \code{FALSE}.
#' @param start				Should we start searching immediately (default is \code{TRUE}).
#' @param num_cores 		The number of CPU cores you wish to use during the search. The default is \code{1}.
#' @return					An object of type \code{optimal_experimental_design_search} which can be further operated upon
#' 
#' @author Adam Kapelner
#' @export
initOptimalExperimentalDesignObject = function(
		X = NULL,
		objective = "mahal_dist", 
		Kgram = NULL,
		wait = FALSE, 
		start = TRUE,
		num_cores = 1){
	
	verify_objective_function(objective, Kgram, n)
	
	if (!is.null(Kgram)){
		n = nrow(Kgram)
		p = NA
	} else {
		n = nrow(X)
		p = ncol(X)
	}
	if (n %% 2 != 0){
		stop("Design matrix must have even rows to have equal treatments and controls")
	}
	
	
	if (objective == "abs_sum_diff"){
		#standardize it -- much faster here
		Xstd = standardize_data_matrix(X)
	}
	if (objective == "mahal_dist"){
		if (p < n){
			SinvX = solve(var(X))
		}
	}
	
	#we are about to construct a OptimalExperimentalDesign java object. First, let R garbage collect
	#to clean up previous objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!	
	
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("OptimalExperimentalDesign.OptimalExperimentalDesign")
	.jcall(java_obj, "V", "setNumCores", as.integer(num_cores))
	.jcall(java_obj, "V", "setN", as.integer(n))
	if (objective != "kernel"){
		p = ncol(X)
		.jcall(java_obj, "V", "setP", as.integer(p))
	}
	.jcall(java_obj, "V", "setObjective", objective)
	if (wait){
		.jcall(java_obj, "V", "setWait")
	}	
	
	#feed in the gram matrix if applicable
	if (!is.null(Kgram)){
		setGramMatrix(java_obj, Kgram)
	} else {
		#feed in the data
		for (i in 1 : n){	
			if (objective == "abs_sum_diff"){
				.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), Xstd[i, , drop = FALSE]) #java indexes from 0...n-1
			} else {
				.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), X[i, , drop = FALSE]) #java indexes from 0...n-1
			}
		}
		
		#feed in the inverse var-cov matrix
		if (objective == "mahal_dist"){
			if (p < n){
				for (j in 1 : p){
					.jcall(java_obj, "V", "setInvVarCovRow", as.integer(j - 1), SinvX[j, , drop = FALSE]) #java indexes from 0...n-1
				}
			}
		}
	}
		
	#now return information as an object (just a list)
	optimal_experimental_design_search = list()
	optimal_experimental_design_search$start = start
	optimal_experimental_design_search$wait = wait
	optimal_experimental_design_search$X = X
	optimal_experimental_design_search$n = n
	optimal_experimental_design_search$p = p
	optimal_experimental_design_search$objective = objective
	optimal_experimental_design_search$java_obj = java_obj
	class(optimal_experimental_design_search) = "optimal_experimental_design_search"
	#if the user wants to run it immediately...
	if (start){
		startSearch(optimal_experimental_design_search)
	}
	#return the final object
	optimal_experimental_design_search
}

#' Prints a summary of a \code{optimal_experimental_design_search} object
#' 
#' @param x			The \code{optimal_experimental_design_search} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print optimal_experimental_design_search
#' @export
print.optimal_experimental_design_search = function(x, ...){
	progress = .jcall(x$java_obj, "D", "progress")
	time_elapsed = searchTimeElapsed(x)
	if (progress == 0){
		cat("No progress on the OptimalExperimentalDesign. Did you run \"startOptimalSearch?\"\n")
	} else if (progress < 1){
		cat("The search is", round(progress * 100, 2), "% complete.\n")
	} else {
		cat("The search is complete.\n")
	}
}

#' Prints a summary of a \code{optimal_experimental_design_search} object
#' 
#' @param object		The \code{optimal_experimental_design_search} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary optimal_experimental_design_search
#' @export
summary.optimal_experimental_design_search = function(object, ...){
	print(object, ...)
}

#' Returns the results (thus far) of the optimal design search
#' 
#' @param obj 				The \code{optimal_experimental_design} object that is currently running the search
#' @param num_vectors		How many allocation vectors you wish to return. The default is 1 meaning the best vector. If \code{Inf},
#' 							it means all vectors.
#' @param form				Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's.
#' 
#' @author Adam Kapelner
#' @export
resultsOptimalSearch = function(obj, num_vectors = 1, form = "one_zero"){
	obj_vals = .jcall(obj$java_obj, "[D", "getAllObjectiveVals", simplify = TRUE)
	ordered_indices = order(obj_vals)
	
	if (num_vectors == Inf){
		num_vectors = length(ordered_indices)
	}
	
	indicTs = .jcall(obj$java_obj, "[[I", "getIndicTs", as.integer(ordered_indices[1 : num_vectors] - 1), simplify = TRUE)
	if (form == "pos_one_min_one"){
		indicTs = (indicTs - 0.5) * 2
	}
	list(
		opt_obj_val = obj_vals[ordered_indices[1]],
		ordered_obj_vals = obj_vals[ordered_indices],
		obj_vals_unordered = obj_vals,
		orig_order = ordered_indices,
		indicTs = indicTs
	)
}
