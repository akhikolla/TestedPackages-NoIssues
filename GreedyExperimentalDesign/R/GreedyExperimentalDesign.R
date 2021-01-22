#' A tool to find a priori experimental designs with good balance greedily
#'
#' @name 		GreedyExperimentalDesign
#' @docType 	package
#' @title 		Greedy Experimental Design Search
#' @author 		Adam Kapelner \email{kapelner@@qc.cuny.edu}
#' @references 	Kapelner, A
#' @keywords 	design optimize
#' @import      rJava stats graphics grDevices GreedyExperimentalDesignJARs
#' @importFrom  checkmate assertChoice assertTRUE assertSetEqual assertCount assertLogical assertNumeric assertClass vname makeAssertion
#' @importFrom 	nbpMatching distancematrix nonbimatch
#' @importFrom  Rcpp sourceCpp evalCpp
#' @useDynLib 	GreedyExperimentalDesign
##### Run "library(roxygen2); roxygenise("GreedyExperimentalDesign", clean = TRUE)" to regenerate all Rd files and NAMESPACE and DESCRIPTION file
##### but make sure you are in the root directory of the project. Make sure to add these two to namespace afterwards:
NULL