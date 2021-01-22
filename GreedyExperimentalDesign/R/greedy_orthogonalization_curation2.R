#' Curate More Orthogonal Vectors Greedily
#' 
#' This function takes a set of allocation vectors and pares them down one-by-one
#' by eliminating the vector that can result in the largest reduction in Avg[ |r_ij| ].
#' It is recommended to begin with a set of unmirrored vectors for speed. Then add the mirrors later
#' for whichever subset you wish.
#' 
#' @param W 		A matrix in ${-1, 1}^{R x n}$ which have R allocation vectors for an experiment of sample size n.
#' @param R0 		The minimum number of vectors to consider in a design. The default is the true bottom, two.
#' @param verbose 	Default is \code{FALSE} but if not, it will print out a message for each iteration.
#' @return 			A list with two elements: (1) \code{avg_abs_rij_by_R} which is a data frame with R - Rmin + 1 rows and 
#' 					columns R and average absolute r_ij and (2) \code{Wsorted} which provides the collection of vectors in
#' 					sorted by best average absolute r_ij in row order from best to worst.
#' 
#' @author Adam Kapelner
#' @export
greedy_orthogonalization_curation2 = function(W, R0 = 100, verbose = FALSE){
	assertClass(W, "matrix")
	assertCount(R0)
	R = nrow(W)
	assertNumeric(R0, upper = R)
	assertLogical(verbose)
	
	Rij = abs(W %*% t(W)) / ncol(W)
	rownames(Rij) = 1 : R
	diag(Rij) = 0 #for speed
	Rij_colsums = colSums(Rij)
	
	left_outs = array(NA, R - R0)
	for (iter_num in 1 : (R - R0)){
		if (iter_num %% 100 == 0){
			if (verbose){cat("iter", iter_num, "\n")}
		}
		#record the min
		left_star = which.max(Rij_colsums)
		Rij_colsums[left_star] = NA
		Rij_colsums = Rij_colsums - Rij[left_star, ]
		left_outs[iter_num] = left_star
	}
	#now return the subset of original allocations in order of elimination
	W[-left_outs, ]
}
