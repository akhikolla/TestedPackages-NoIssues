##' Curate More Orthogonal Vectors Via Numerical Optimization
##' 
##' This function takes a set of allocation vectors and pares them down
##' by eliminating the vectors that can result in the largest reduction in Avg[ |r_ij| ] via
##' Gurobi's numerical optimization for binary integer quadratic programs.
##' 
##' Make sure you setup Gurobi properly first. This means applying for a license, downloading, installing, then make sure 
##' you've registered the license on your computer using the \code{grbgetkey} command. Then you need to install the R package
##' from the local file (it is not on CRAN) via e.g. \code{install.packages("c:/gurobi902/win64/R/gurobi_9.0-2.zip", repos = NULL)}.
##' 
##' @param W 		A matrix in the set ${-1, 1}^{R x n}$ which have R allocation vectors for an experiment of sample size n.
##' @param R0 		The number of vectors in the pared-down design.
##' @param time_limit_seconds	The limit on how long in seconds the binary integer quadratic program should run (default is 2 * 60s).
##' @param num_cores				Number of cores to use during search. Default is \code{2}.
##' @param gurobi_params			A list of optional parameters to be passed to Gurobi (see their documentation online).
##' @param verbose				Should Gurobi log to console? Default is \code{TRUE}.
##' @return 			A list object which houses the results from Gurobi. Depending on the \code{gurobi_parms},
##' 					the data within will be different. The most relevant tags are \code{x} for the best found solution and \code{objval}
##' 					for the object
##' 
##' @author Adam Kapelner
#
#####EXPORT
#gurobi_numerical_optimization_orthogonalization_curation = function(W, R0, verbose = TRUE, num_cores = 2, time_limit_seconds = 2 * 60, gurobi_params = list()){
#	assertClass(W, "matrix")
#	R = nrow(W)
#	assertCount(R0, positive = TRUE)
#	assertNumeric(R0, upper = R)
#	assertCount(num_cores, positive = TRUE)
#	assertNumeric(time_limit_seconds, lower = 1)
#	assertLogical(verbose)
#	assertClass(gurobi_params, "list")
#	
#	n = ncol(W)
#	Rij = abs(W %*% t(W)) / n
#	diag(Rij) = 0
#
#	#  minimize
#	#        [s_1 s_2 ... s_R] Rij [s_1 s_2 ... s_R]^T
#	#  subject to
#	#        s_1 + s_2 + ... + s_R = R0,
#	#        s_1, s_2, ..., s_R binary
#	
#	model 		= list()
#	#quadratic component
#	model$Q     = Rij / (R0^2 - R0)
#	model$obj   = t(as.matrix(rep(0, R)))
#	#constraint component
#	model$A     = t(as.matrix(rep(1, R)))
#	model$rhs   = R0
#	model$sense = "="
#	model$vtype = "B"
#	gurobi_params$LogToConsole = as.numeric(verbose)
#	gurobi_params$Threads = num_cores
#	gurobi_params$TimeLimit = time_limit_seconds	
#	gurobi::gurobi(model, gurobi_params)
#}