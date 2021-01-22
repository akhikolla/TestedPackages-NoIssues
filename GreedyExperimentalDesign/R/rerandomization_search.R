#' Begin a Rerandomization Search
#' 
#' This method creates an object of type rerandomization_experimental_design and will immediately initiate
#' a search through $1_{T}$ space for forced-balance designs.
#' 
#' @param X							The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 									(one for each measurement on the subject). This is the design matrix you wish 
#' 									to search for a more optimal design.
#' @param obj_val_cutoff_to_include	Only allocation vectors with objective values lower than this threshold will be returned.
#' 									If the cutoff is infinity, you are doing BCRD and you should use the \code{complete_randomization_with_forced_balanced}
#' 									function instead.
#' @param max_designs 				The maximum number of designs to be returned. Default is 10,000. Make this large 
#' 									so you can search however long you wish as the search can be stopped at any time by
#' 									using the \code{\link{stopSearch}} method 
#' @param objective					The objective function to use when searching design space. This is a string
#' 									with valid values "\code{mahal_dist}" (the default), "\code{abs_sum_diff}" or "\code{kernel}".
#' @param Kgram						If the \code{objective = kernel}, this argument is required to be an \code{n x n} matrix whose
#' 									entries are the evaluation of the kernel function between subject i and subject j. Default is \code{NULL}.
#' @param wait						Should the \code{R} terminal hang until all \code{max_designs} vectors are found? The 
#' 									default is \code{FALSE}.
#' @param start						Should we start searching immediately (default is \code{TRUE}).
#' @param num_cores 				The number of CPU cores you wish to use during the search. The default is \code{1}.
#' @return							An object of type \code{rerandomization_experimental_design_search} which can be further operated upon.
#' 
#' @author Adam Kapelner
#' @export
initRerandomizationExperimentalDesignObject = function(
		X = NULL, 
		obj_val_cutoff_to_include,
		max_designs = 1000,
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
	
	#we are about to construct a RerandomizationExperimentalDesign java object. First, let R garbage collect
	#to clean up previous RerandomizationExperimentalDesign objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!	
	
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("RerandomizationExperimentalDesign.RerandomizationExperimentalDesign")
	.jcall(java_obj, "V", "setMaxDesigns", as.integer(max_designs))
	if (!is.null(obj_val_cutoff_to_include)){
		.jcall(java_obj, "V", "setObjValCutoffToInclude", as.numeric(obj_val_cutoff_to_include))
	}
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
		#feed in the raw data
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
	rerandomization_experimental_design_search = list()
	rerandomization_experimental_design_search$max_designs = max_designs
	rerandomization_experimental_design_search$obj_val_cutoff_to_include = obj_val_cutoff_to_include
	rerandomization_experimental_design_search$start = start
	rerandomization_experimental_design_search$wait = wait
	rerandomization_experimental_design_search$X = X
	rerandomization_experimental_design_search$n = n
	rerandomization_experimental_design_search$p = p
	rerandomization_experimental_design_search$objective = objective
	rerandomization_experimental_design_search$java_obj = java_obj
	class(rerandomization_experimental_design_search) = "rerandomization_experimental_design_search"
	#if the user wants to run it immediately...
	if (start){
		startSearch(rerandomization_experimental_design_search)
	}
	#return the final object
	rerandomization_experimental_design_search
}

#' Prints a summary of a \code{rerandomization_experimental_design_search} object
#' 
#' @param x			The \code{rerandomization_experimental_design_search} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print rerandomization_experimental_design_search
#' @export
print.rerandomization_experimental_design_search = function(x, ...){
	progress = rerandomizationSearchCurrentProgress(x)
	time_elapsed = searchTimeElapsed(x)
	if (progress == 0){
		cat("No progress on the RerandomizationExperimentalDesign. Did you run \"startSearch?\"\n")
	} else if (progress == x$max_designs){
		cat("The search completed in", time_elapsed, "seconds.", progress, "vectors have been found.\n")
	} else {
		cat("The search has found ", progress, " vectors thus far (", round(progress / x$max_designs * 100), "%) in ", time_elapsed," seconds.\n", sep = "")
	}
}

#' Prints a summary of a \code{rerandomization_experimental_design_search} object
#' 
#' @param object		The \code{rerandomization_experimental_design_search} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary rerandomization_experimental_design_search
#' @export
summary.rerandomization_experimental_design_search = function(object, ...){
	print(object, ...)
}

# Returns the number of vectors found by the rerandomization design search
# 
# @param obj 		The \code{rerandomization_experimental_design} object that is currently running the search
# 
# @author Adam Kapelner
rerandomizationSearchCurrentProgress = function(obj){
	.jcall(obj$java_obj, "I", "progress")
}


#' Returns the results (thus far) of the rerandomization design search
#' 
#' @param obj 					The \code{rerandomization_experimental_design} object that is currently running the search
#' @param include_assignments	Do we include the assignments (takes time) and default is \code{FALSE}.
#' @param form					Which form should the assignments be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' 
#' @author Adam Kapelner
#' @export
resultsRerandomizationSearch = function(obj, include_assignments = FALSE, form = "one_zero"){
	obj_vals = .jcall(obj$java_obj, "[D", "getObjectiveVals")
	
	ending_indicTs = NULL
	if (include_assignments){
		ending_indicTs = .jcall(obj$java_obj, "[[I", "getEndingIndicTs", simplify = TRUE)
		if (form == "pos_one_min_one"){
			ending_indicTs = (ending_indicTs - 0.5) * 2
		}
	}	
	
	rerandomization_experimental_design_search_results = list(
		obj_vals = obj_vals, 
		ending_indicTs = ending_indicTs
	)
	class(rerandomization_experimental_design_search_results) = "rerandomization_experimental_design_search_results"
	#return the final object
	rerandomization_experimental_design_search_results
}
#}