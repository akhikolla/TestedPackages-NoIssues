####################################################################################################################################
### Filename:    S3methods.R
### Description: S3 Generic Methods for calculating the test statistic; Providing different ways to to 'input' the data
###              and methods for print/summary are defined here
###
###
###
####################################################################################################################################

#' Test for Multi-Factor High-Dimensional Repeated Measures
#'
#' @description Performing main and interaction effects of up to three whole- or subplot-factors. In total, a maximum of four factors can be used. There are two different S3 methods available. The first method requires a list of matrices in the wide table format. The second methodl requres a data.frame in the long table format.
#' @rdname hrm_test
#' @param data Either a data.frame (one observation per row) or a list with matrices (one subject per row) for all groups containing the data
#' @param formula A model formula object. The left hand side contains the response variable and the right hand side contains the whole- and subplot factors.
#' @param subject column name within the data frame X identifying the subjects
#' @param variable if not 'NULL' then multivariate tests are applied. We assume that for each factor level of 'variable', we observe several repated measurements. Currently only supports designs with 1 whole- and one sub-plot factor.
#' @param alpha alpha level used for calculating the critical value for the test
#' @param nonparametric Logical variable indicating wether the noparametric version of the test statistic should be used
#' @param np.correction Logical variable indicating wether a small sample size correction for the nonparametric test should be used (TRUE) or not (FALSE). By using NA, np.correction is used automatically in an high-dimensional setting.
#' @param character.only a logical indicating whether subject can be assumed to be a character string
#' @param ... Further arguments passed to 'hrm_test' will be ignored
#' @return Returns an object from class HRM containing
#' @return \item{result}{A dataframe with the results from the hypotheses tests.}
#' @return \item{formula}{The formula object which was used.}
#' @return \item{alpha}{The type-I error rate which was used.}
#' @return \item{subject}{The column name identifying the subjects.}
#' @return \item{factors}{A list containing the whole- and subplot factors.}
#' @return \item{data}{The data.frame or list containing the data.}
#' @example R/example_1.txt
#' @keywords export
hrm_test <- function(data, ...) {
  UseMethod("hrm_test")
}

#' @method hrm_test default
#' @keywords export
hrm_test.default <- function(data) {
  stop("Your data needs either to be a data.frame or a list containing matrices (a matrix for each group).")
}

#' @rdname hrm_test
#' @method hrm_test list
#' @keywords export
hrm_test.list <- function(data, alpha = 0.05, ...) {
  return(hrm.test.matrices(data, alpha))
}

#' @rdname hrm_test
#' @method hrm_test data.frame
#' @keywords export
hrm_test.data.frame <- function(data, formula, alpha = 0.05,  subject, variable = NULL, nonparametric = FALSE,
                                np.correction = NA, character.only = FALSE, ... ) {
  if(!character.only) {
    subject <- as.character(substitute(subject))
    variable <- as.character(substitute(variable))
    if(length(variable) == 0) {
      variable <- NULL
    }
  }

  return(hrm_test_internal(formula=formula, data=data, alpha=alpha,subject=subject, variable=variable,
                           nonparametric, np.correction ))
}


#' @keywords export
print.HRM <- function(x, ...) {
  if(x$nonparametric & is.null(x$variable)){
    cat("Nonparametric Repeated Measures Analysis\n")
    if(x$np.correction) {
      cat("With Bias Correction for Degrees of Freedom\n")
    } else if(!x$np.correction) {
      cat("Without Bias Correction for Degrees of Freedom\n")
    }
    cat("\n")
  }
  if(!is.null(x$variable)) {
    if(x$nonparametric) {
      cat("Nonparametric Multivariate Repeated Measures Analysis\n")
      cat("\n")
    } else {
      cat("Multivariate Repeated Measures Analysis\n")
      cat("\n")
    }
  }
  if(!is.null(x$formula)) {
    cat("Call:", "\n")
    print(x$formula)
    cat("\n")
  }

  print(x$result, row.names = FALSE)
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  cat("\n")
}


#' @keywords export
summary.HRM <- function(object, ...) {
  if(object$nonparametric & is.null(object$variable)) {
    cat("Nonparametric Repeated Measures Analysis\n")
    if(object$np.correction) {
      cat("With Bias Correction for Degrees of Freedom\n")
    } else if(!object$np.correction) {
      cat("Without Bias Correction for Degrees of Freedom\n")
    }
  }
  cat("\n")
  cat("Summary:\n")
  cat("\n")
  if(!is.null(object$formula)) {
    cat("Call:", "\n")
    print(object$formula)
    cat("\n")
    cat("between-subject factors: ")
    cat(object$factors[[1]], sep = ", ")
    cat("\n")
    cat("within-subject factors: ")
    cat(object$factors[[2]], sep = ", ")
    cat("\n")
    cat("\n")
  }

  print(object$result, row.names = FALSE)
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  cat("\n")
}
