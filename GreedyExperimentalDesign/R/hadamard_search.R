#' Create a Hadamard Design
#' 
#' This method returns unique designs according to a Hadamard matrix
#' 
#' @param X				The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 						(one for each measurement on the subject). The measurements aren't used to
#' 						compute the Hadamard designs, only the number of rows.
#' @param strict		Hadamard matrices are not available for all $n$.
#' @param form			Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' 
#' @return				An matrix of dimension $R$ x $n$ where $R$ is the number of Hadamard allocations.
#' 
#' @author Adam Kapelner
#' @export
hadamardExperimentalDesign = function(X, strict = TRUE, form = "zero_one"){
	assertMatrixWithEvenNumRows(X)
	n = nrow(X)
	W = t(survey::hadamard(n))[-1, ] #the first row is always just 1's, so remove it
	#also retain only those that are forced balance (I think this is all of them, but I want to just make sure
	W = W[rowMeans(W) == 0.5, ]
	if (strict){
		if (ncol(W) != n){
			stop("A Hadamard matrix for this sample size doesn't exist. Use strict = FALSE to approximate one randomly.")
		}		
	} else {
		Wcut = matrix(NA, nrow = nrow(W), ncol = n)
		for (r in 1 : nrow(W)){
			w = W[r, ]
			idx_Ts = which(w == 1)
			idx_Cs = which(w == 0)
			idx_trimmed = sort(c(
				sample(idx_Ts, n / 2), 
				sample(idx_Cs, n / 2)
			))
			Wcut[r, ] = w[idx_trimmed]
		}
		W = Wcut #overwrite
	}
	
	if (form == "pos_one_min_one"){
		W = (W - 0.5) * 2
	}
	#now return only those that are unique
	unique(W)
}

assertMatrixWithEvenNumRows = function(x, .var.name = vname(x), add = NULL) {
	assertClass(x, "matrix")
	makeAssertion(x, (function(x){if (nrow(x) %% 2 != 0) "Matrix must have an even number of rows" else TRUE})(x), .var.name, add)
}




