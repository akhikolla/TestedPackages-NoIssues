#' StMoSim: Plots a QQ-Norm Plot with Several Gaussian Simulations
#'
#' With this package you can simulate several lines into the QQ-Norm Plot under the assumption of Gaussian distribution. 
#' If the realised observations lie inside of the simulations tracks there is the possibility that the observations stem from a Gaussian distribution.
#' This can be very useful in residual analysis where you have to evaluate whether the model residuals fit the assumption of gaussian distributed terms or not.
#'
#' @section Changelog:
#' 
#' -----------<CHANGELOG>-----------
#'   
#' -------< v3.1.1 - 2018-11-19 >-------
#'   
#'  provide more (plot) arguments to the user.   
#'  
#'  updated documentation - added more expamples.
#'  
#'  added BugReports argument in DESCRIPTION.
#'  
#'  implemented all recommendations from RcppParallel package.
#'   
#' -------< v3.1 - 2018-11-13 >-------
#'   
#'  Minor bug fixes, due to CHECK changes on CRAN.   
#'  
#'  Moved documentation to roxygen2.
#'   
#' -------< v3.0 - 2014-10-16 >-------
#'   
#'  Computation intense code moved to C++. 
#' 
#'  Moved to parallel computation, thanks to Rcpp/RcppParallel !
#'   
#'  Minor bug fixes.
#' 
#' -------< v2.2 - 2012-02-24 >-------
#'   
#'  Minor bug fixes, due to CHECK changes on CRAN.
#' 
#' -------< v2.1 - 2012-02-24 >-------
#'   
#'  Minor bug fixes.
#' 
#' -------<v2.0 - 2011-03-31 >-------
#'   
#'  Moved to S4 Classes.
#' 
#' -------<v1.1 - 2010-05-03 >-------
#'   
#'  First Version on CRAN.
#' 
#' -----------</CHANGELOG>-----------
#'
#' @keywords package
#' 
#' @author Matthias Salvisberg <matthias.salvisberg@@gmail.com>
#'
#' @docType package
#' @name StMoSim
NULL