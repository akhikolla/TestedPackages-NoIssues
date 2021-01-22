#' ludic
#' 
#' Linkage Using Diagnosis Codes
#' 
#' This package implements probabilistic record linkage methods that relies on the use of diagnosis 
#' codes only, in the absence of direct identifiers .
#' 
#' \tabular{ll}{
#' Package: \tab ludic\cr
#' Type: \tab Package\cr
#' Version: \tab ludic 0.1.8\cr
#' Date: \tab 2019-12-04\cr
#' License: \tab \href{https://cran.r-project.org/web/licenses/MIT}{The "MIT License" (MIT)}\cr
#' }
#' The main function of \code{ludic} is \code{\link{recordLink}}.
#' 
#' @author Boris P. Hejblum, Tianxi Cai
#' --- Maintainer: Boris P. Hejblum
#' 
#' @references Hejblum BP, Weber G, Liao KP, Palmer N, Churchill S, Szolovits P, Murphy S, Kohane I, Cai T 
#' Probabilistic Record Linkage of De-Identified Research Datasets Using Diagnosis Codes, \emph{submitted}, 2017.
#' 
#' @docType package
#' @name ludic-package
#' @aliases ludic
#' 
#' @useDynLib ludic, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' 
NULL