#' Implements forced balanced randomization
#' 
#' @param n 		number of observations
#' @param r 		number of randomized designs you would like
#' @param form		Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' @return 			a matrix where each column is one of the \code{r} designs
#' 
#' @author Adam Kapelner
#' @export
complete_randomization_with_forced_balanced = function(n, r, form = "one_zero"){
	indicTs = matrix(NA, nrow = r, ncol = n)
	zero_one_vec = c(rep(0, n / 2), rep(1, n / 2))
	for (nsim in 1 : r){
		indicTs[nsim, ] = sample(zero_one_vec)
	}
	if (form == "pos_one_min_one"){
		indicTs = (indicTs - 0.5) * 2
	}
	indicTs
}


#' Implements complete randomization (without forced balance)
#' 
#' @param n 		number of observations
#' @param r 		number of randomized designs you would like
#' @param form		Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' @return 			a matrix where each column is one of the \code{r} designs
#' 
#' @author Adam Kapelner
#' @export
complete_randomization = function(n, r, form = "one_zero"){
	indicTs = matrix(NA, nrow = r, ncol = n)
	
	for (nsim in 1 : r){
		indicTs[nsim, ] = rbinom(n, 1, 0.5)
	}
	if (form == "pos_one_min_one"){
		indicTs = (indicTs - 0.5) * 2
	}
	indicTs
}
