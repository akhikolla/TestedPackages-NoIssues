#' Begin Karp Search
#' 
#' This method creates an object of type karp_experimental_design and will immediately initiate
#' a search through $1_{T}$ space. Note that the Karp search only works 
#' for one covariate (i.e. $p=1$) and the objective "abs_sum_diff".
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more karp design.
#' @param balanced			Should the final vector be balanced? Default and recommended is \code{TRUE}.
#' @param wait				Should the \code{R} terminal hang until all \code{max_designs} vectors are found? The 
#' 							deafult is \code{FALSE}.
#' @param start				Should we start searching immediately (default is \code{TRUE}).
#' @return					An object of type \code{karp_experimental_design_search} which can be further operated upon
#' 
#' @author Adam Kapelner
#' @export
initKarpExperimentalDesignObject = function(X,
		wait = FALSE, 
		balanced = TRUE,
		start = TRUE){
	n = nrow(X)
	if (n %% 2 != 0 && balanced){
		stop("Design matrix must have even rows to have equal treatments and controls if you are requiring balance.")
	}
	p = ncol(X)
	
	if (p > 1){
		stop("Karp search only works for p = 1.")
	}
	
	#we are about to construct a KarpExperimentalDesign java object. First, let R garbage collect
	#to clean up previous objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!	
	
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("KarpExperimentalDesign.KarpExperimentalDesign")
	.jcall(java_obj, "V", "setN", as.integer(n))
	.jcall(java_obj, "V", "setP", as.integer(p))
	if (wait){
		.jcall(java_obj, "V", "setWait")
	}
	if (balanced){
		.jcall(java_obj, "V", "setBalanced")
	}
	
	#feed in the data
	for (i in 1 : n){		
		.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), X[i, , drop = FALSE]) #java indexes from 0...n-1
	}
	
	#now return information as an object (just a list)
	karp_experimental_design_search = list()
	karp_experimental_design_search$start = start
	karp_experimental_design_search$wait = wait
	karp_experimental_design_search$balanced = balanced
	karp_experimental_design_search$X = X
	karp_experimental_design_search$n = n
	karp_experimental_design_search$java_obj = java_obj
	class(karp_experimental_design_search) = "karp_experimental_design_search"
	#if the user wants to run it immediately...
	if (start){
		startSearch(karp_experimental_design_search)
	}
	#return the final object
	karp_experimental_design_search
}

#' Prints a summary of a \code{karp_experimental_design_search} object
#' 
#' @param x			The \code{karp_experimental_design_search} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print karp_experimental_design_search
#' @export
print.karp_experimental_design_search = function(x, ...){
	progress = .jcall(x$java_obj, "D", "progress")
	time_elapsed = searchTimeElapsed(x)
	if (progress == 0){
		cat("No progress on the KarpExperimentalDesign. Did you run \"startKarpSearch?\"\n")
	} else if (progress < 1){
		cat("The search is", round(progress * 100, 2), "% complete.\n")
	} else {
		cat("The search is complete.\n")
	}
}

#' Prints a summary of a \code{karp_experimental_design_search} object
#' 
#' @param object		The \code{karp_experimental_design_search} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary karp_experimental_design_search
#' @export
summary.karp_experimental_design_search = function(object, ...){
	print(object, ...)
}

#' Returns the results (thus far) of the karp design search
#' 
#' @param obj 			The \code{karp_experimental_design} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @export
resultsKarpSearch = function(obj){
	obj_val = .jcall(obj$java_obj, "D", "getKarpObjectiveVal")
	indicT = .jcall(obj$java_obj, "[I", "getKarpIndicT", simplify = TRUE)
	list(
		obj_val = obj_val,
		indicT = indicT
	)
}
