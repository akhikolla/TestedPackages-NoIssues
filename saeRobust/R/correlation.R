#' Correlation Structure
#'
#' Various correlation structures. They can be used inside the \link{rfh}
#' function to supply an alterantive variance structure to be fitted. For
#' examples see the documentation of \link{rfh}.
#'
#' @details \code{corSAR1} can be used to model a simultanous autoregressive
#'   process of order one: spatial correlation.
#'
#' @slot W the row-standardised proximity matrix
#'
#' @rdname correlation
#' @export
corSAR1(W ~ matrix | Matrix) %type% {
  .Object
}

#' @name corSAR1
#' @usage corSAR1(W, ...)
#'
#' @param W the row-standardised proximity matrix
#' @param ... arguments passed to \code{new}. In the case of corSAR1AR1
#'   arguments \code{W} and \code{nTime} are expected.
#'
#' @rdname correlation
#' @export corSAR1
corSAR1

#' @details \code{corAR1} can be used to model a autoregressive
#'   process of order one: temporal correlation.
#'
#' @slot nTime (numeric) number of time periods
#'
#' @rdname correlation
#' @export
corAR1(nTime ~ numeric | integer) %type% {
  .Object
}

#' @name corAR1
#' @usage corAR1(nTime, ...)
#'
#' @param nTime (numeric) number of time periods
#'
#' @rdname correlation
#' @export corAR1
corAR1

#' @rdname correlation
#' @export
corAR1 : corSAR1 : corSAR1AR1() %type% .Object

#' @details \code{corSAR1AR1} can be used to model to model spatial and temporal
#'   correlation
#'
#' @name corSAR1AR1
#' @usage corSAR1AR1(...)
#'
#' @rdname correlation
#' @export corSAR1AR1
corSAR1AR1
