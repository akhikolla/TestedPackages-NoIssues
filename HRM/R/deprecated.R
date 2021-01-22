####################################################################################################################################
### Filename:    deprecated.R
### Description: Deprecated functions for compatibility with older versions; will be removed in of the next updates
### 
### 
###
###
####################################################################################################################################


#' @title Deprecated functions in package \pkg{HRM}.
#' @description The functions listed below are deprecated and will be defunct in
#'   the near future. When possible, alternative functions with similar
#'   functionality are also mentioned. Help pages for deprecated functions are
#'   available at \code{help("hrm.test.matrix-deprecated")} and \code{help("hrm.test.dataframe-deprecated")}.
#' @name HRM-deprecated
#' @keywords internal
NULL


#' @title hrm.test.matrix
#' @description  The function listed below is deprecated and will be defunct in the near future. Please use 'hrm.test' instead.
#' @usage hrm.test.matrix(data, alpha = 0.05)
#' @param data A list of matrices
#' @param alpha alpha level used for the test
#' @name hrm.test.matrix-deprecated
#' @seealso \code{\link{HRM-deprecated}}
#' @keywords internal
NULL


#' @title hrm.test.dataframe
#' @description  The function listed below is deprecated and will be defunct in the near future. Please use 'hrm.test' instead.
#' @usage hrm.test.dataframe(data, alpha = 0.05, group , subgroup,
#' factor1, factor2, factor3, subject, response)
#' @param data A data.frame containing the data
#' @param alpha alpha level used for the test
#' @param subject column name within the data frame identifying the subjects
#' @param group column name within the data frame specifying the first whole-plot factor
#' @param subgroup column name within the data frame specifying the second whole-plot factor
#' @param factor1 column name within the data frame specifying the first sub-plot factor
#' @param factor2 column name within the data frame specifying the second sub-plot factor
#' @param factor3 column name within the data frame specifying the third sub-plot factor
#' @param response column name within the data frame specifying the measurements
#' @name hrm.test.dataframe-deprecated
#' @seealso \code{\link{HRM-deprecated}}
#' @keywords internal
NULL


#' @rdname HRM-deprecated
#' @keywords export
hrm.test.matrix <- function(data, alpha = 0.05){
  .Deprecated("hrm_test", package=NULL,
              old = as.character(sys.call(sys.parent()))[1L])
  hrm_test(data, alpha)
}


#' @rdname HRM-deprecated
#' @details For \code{hrm.test.dataframe} and \code{hrm.test.matrix} use \code{\link{hrm_test}} instead.
#' @keywords export
hrm.test.dataframe <- function(data, alpha = 0.05, group , subgroup, factor1, factor2, factor3, subject, response ){
    
  .Deprecated("hrm_test", package=NULL,
               old = as.character(sys.call(sys.parent()))[1L])

  # generate formula object
  factors <- c()
  
  if(!missing(group)){
    factors <- paste( c(factors, group), collapse = "*" )
  }
  if(!missing(subgroup)){
    factors <- paste( c(factors, subgroup), collapse = "*" )
  }
  if(!missing(factor1)){
    factors <- paste( c(factors, factor1), collapse = "*" )
  }
  if(!missing(factor2)){
    factors <- paste( c(factors, factor2), collapse = "*" )
  }
  if(!missing(factor3)){
    factors <- paste( c(factors, factor3), collapse = "*" )
  }
  form <- as.formula(paste(response, factors, sep=" ~ "), env = parent.frame()) 

  # calculate test statistics with method hrm_test
  hrm_test(formula = form, alpha = alpha, subject = subject, data = data)
}


#' Graphical User Interface for Testing Multi-Factor High-Dimensional Repeated Measures
#' 
#' @description Graphical User Interface (R Package RGtk2 needed) for the Function 'hrm_test': Test for main effects and interaction effects of one or two between-subject factors and one, two or three within-subject factors (at most four factors can be used).
#' @return The results can be saved as LaTeX Code or as plain text. Additionally a plot of the group profiles an be saved when using one whole- and one subplot factor. 
#' @keywords internal
hrm.GUI <- function(){
  .Deprecated("hrm_GUI", package=NULL,
              old = as.character(sys.call(sys.parent()))[1L])
  hrm_GUI()
}
