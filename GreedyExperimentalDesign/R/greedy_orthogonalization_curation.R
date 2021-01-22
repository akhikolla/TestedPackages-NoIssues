#' Curate More Orthogonal Vectors Greedily
#' 
#' This function takes a set of allocation vectors and pares them down one-by-one
#' by eliminating the vector that can result in the largest reduction in Avg[ |r_ij| ].
#' It is recommended to begin with a set of unmirrored vectors for speed. Then add the mirrors later
#' for whichever subset you wish.
#' 
#' @param W 		A matrix in ${-1, 1}^{R x n}$ which have R allocation vectors for an experiment of sample size n.
#' @param Rmin 		The minimum number of vectors to consider in a design. The default is the true bottom, two.
#' @param verbose 	Default is \code{FALSE} but if not, it will print out a message for each iteration.
#' @return 			A list with two elements: (1) \code{avg_abs_rij_by_R} which is a data frame with R - Rmin + 1 rows and 
#' 					columns R and average absolute r_ij and (2) \code{Wsorted} which provides the collection of vectors in
#' 					sorted by best average absolute r_ij in row order from best to worst.
#' 
#' @author Adam Kapelner
#' @export
greedy_orthogonalization_curation = function(W, Rmin = 2, verbose = FALSE){
	assertClass(W, "matrix")
	assertCount(Rmin)
	R = nrow(W)
	assertNumeric(Rmin, upper = R)
	
	if (R > 1000){
		warning("If number of designs > 1000, it goes pretty slow.\n")
	}
	
	n = ncol(W)
	Rij = abs(W %*% t(W)) / n
	rownames(Rij) = 1 : R
	diag(Rij) = NA

	avg_abs_rijss = mean(c(Rij), na.rm = TRUE)
	if (verbose){cat("iter", 0, "avg abs rijs", round(avg_abs_rijss[1], 4), "\n")}
	eliminations = c()
	for (iter_num in 1 : (R - Rmin)){
		#leave each out
		avg_abs_rijss_left = compute_avg_abs_rijss_left(Rij)
		#record the min
		left_star = which.min(avg_abs_rijss_left)
		avg_abs_rijs_greedy = min(avg_abs_rijss_left)
		#if the min is greater than what happened last time, we jump ship
		if (avg_abs_rijs_greedy > avg_abs_rijss[length(avg_abs_rijss)]){
			break
		}
		#otherwise record this in our list and keep going
		#kill the one that would yield the lowest in objective function
		previous_observations_included = rownames(Rij)
		Rij = Rij[-left_star, -left_star, drop = FALSE]
		i_elim = setdiff(previous_observations_included, rownames(Rij))
		eliminations = c(eliminations, i_elim)
		avg_abs_rijss = c(avg_abs_rijss, avg_abs_rijs_greedy)
		if (verbose){cat("iter", iter_num, "obs eliminated:", i_elim, "avg abs rijs", round(tail(avg_abs_rijss, 1), 4), "\n")}
	}
	#first elimination is last
	eliminations = rev(as.numeric(eliminations))
	#now return the original allocations in order of elimination as well as the rho reductions
	list(
		avg_abs_rij_by_R = 
			data.frame(R = (R : Rmin)[1 : length(avg_abs_rijss)], avg_abs_rijs = avg_abs_rijss),
		Wsorted = 
			W[c(as.numeric(rownames(Rij)), eliminations), ]
	)
}
