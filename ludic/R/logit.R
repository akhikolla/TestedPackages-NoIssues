#' Logit link function
#' 
#' inverse of the logistic transformation function
#' 
#' @param p a real number in [0,1] to be transformed to the real scale 
#' ( ]-infinity,infinity[ )
#'
#'@keywords internal
#'@export
logit <- function(p){
  log(p/(1-p))
}


#' Expit link function
#' 
#' logistic transformation function
#' 
#' @param x a real number in ]-infinity,infinity[ to be transformed to the 
#' probability scale ( [0,1] )
#'
#'@rdname logit
#'@keywords internal
#'
#'@export
expit <- function(x){
  exp(x)/(1+exp(x))
}


#' Expit link function first derivative
#' 
#' logistic transformation function
#' 
#' @param x a real number in ]-infinity,infinity[ to be transformed to the 
#' probability scale ( [0,1] )
#'
#'@rdname logit
#'@keywords internal
#'
#'@export
expit_dev1 <- function(x){
  exp(x)/(1+exp(x))^2
}


