## Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
#' Root Mean Square Estimation
#'
#' @param y vector
#' @param ypred vecotr
#'
#' @return numeric number. A measure of difference between y and ypred.
#' @export
#' @example 
#' ypred = smooth_spline(cars$dist, cars$speed)
#' RMSE(cars$speed, ypred)
RMSE = function(y, ypred){
  idx = !is.na(y) & !is.na(ypred)
  sqrt(sum((y[idx]-ypred[idx])^2)/sum(idx))
}

#' Normalized Mean Square Estimation
#'
#' @param y vector
#' @param ypred vector
#'
#' @return numeric number. A measure of difference between y and ypred.
#' @export
#' @example 
#' ypred = smooth_spline(cars$dist, cars$speed)
#' NMSE(cars$speed, ypred)
NMSE = function(y, ypred){
  idx = !is.na(y) & !is.na(ypred)
  sum((y[idx]-ypred[idx])^2)/sum(y[idx]^2)
}

#' Absolute relative error
#'
#' @param y vector
#' @param ypred vector
#'
#' @return numeric number. A measure of difference between y and ypred.
#' @export
#' @example 
#' ypred = smooth_spline(cars$dist, cars$speed)
#' ARE(cars$speed, ypred)
ARE = function(y, ypred){
  idx = !is.na(y) & !is.na(ypred)
  mean(abs(y[idx] - ypred[idx])/y[idx])
}
