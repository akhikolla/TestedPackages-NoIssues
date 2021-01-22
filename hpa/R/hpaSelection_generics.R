#' Predict outcome and selection equation values from hpaSelection model
#' @description This function predicts outcome and selection equation values from hpaSelection model.
#' @param object Object of class "hpaSelection"
#' @param method string value indicating prediction method based on hermite polynomial approximation "HPA" or Newey method "Newey".
#' @template newdata_Template
#' @param is_cond logical; if \code{TRUE} (default) then conditional predictions will be estimated. Otherwise unconditional predictions will be returned.
#' @param is_outcome logical; if \code{TRUE} (default) then predictions for selection equation will be estimated using "HPA" method.
#' Otherwise selection equation predictions (probabilities) will be returned.
#' @template elipsis_Template
#' @details Note that Newey method can't predict conditional outcomes for zero selection equation value. Conditional probabilities for selection equation
#' could be estimated only when dependent variable from outcome equation is observable.
#' @return This function returns the list which structure depends on \code{method}, \code{is_probit} and \code{is_outcome} values.
predict.hpaSelection <- function (object,  ...,
                              newdata = NULL, method = "HPA", 
                              is_cond = TRUE, 
                              is_outcome = TRUE)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  return(predict_hpaSelection(object, newdata, 
                       method, is_cond, 
                       is_outcome))
}
###
#' Summarizing hpaSelection Fits
#' @description This function summarizing hpaSelection Fits
#' @param object Object of class "hpaSelection"
#' @template elipsis_Template
#' @return This function returns the same list as \code{\link[hpa]{hpaSelection}} 
#' function changing its class to "summary.hpaSelection".
summary.hpaSelection <- function (object, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  return(summary_hpaSelection(object))
}
###
#' Summary for "hpaSelection" object
#' @param x Object of class "hpaSelection"
#' @template elipsis_Template
print.summary.hpaSelection <- function (x, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  return(print_summary_hpaSelection(x))
}
###
#' Plot hpaSelection random errors approximated density
#' @param x Object of class "hpaSelection"
#' @param y this parameter currently ignored
#' @param is_outcome logical; if TRUE then function plots the graph for outcome equation random errors. 
#' Otherwise plot for selection equation random errors will be plotted.
#' @template elipsis_Template
#' @return This function returns the list containing random error expected value \code{errors_exp}
#' and variance \code{errors_var} estimates for selection (if \code{is_outcome = TRUE}) or outcome
#' (if \code{is_outcome = FALSE}) equation.
plot.hpaSelection <- function (x, y = NULL, ..., is_outcome = TRUE) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  if (!is.null(y))
  {
    warning("Note that y parameter currently ignored.")   
  }
  return(plot_hpaSelection(x, is_outcome))
}

###
#' Calculates log-likelihood for "hpaSelection" object
#' @description This function calculates log-likelihood for "hpaSelection" object
#' @usage \method{logLik}{hpaSelection}(object, ...)
#' @param object Object of class "hpaSelection"
#' @template elipsis_Template
logLik.hpaSelection <- function (object, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  lnL <- logLik_hpaSelection(object)
  attr(lnL, "class") <- "logLik"
  attr(lnL, "df") <- length(as.vector(object$x1))
  
  return(lnL)
}
###
#' Print method for "hpaSelection" object
#' @param x Object of class "hpaSelection"
#' @template elipsis_Template
print.hpaSelection <- function (x, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  cat(paste("It is the object of class",class(x),"\n"))
  cat("It contains the following elements:\n")
  cat(names(x), sep = ", ")
  cat("\n")
  cat("---\n")
  cat("Estimation results:\n")
  print(x$results)
  cat("---\n")
  cat(paste("Log-likelihood function value is:", round(x$'log-likelihood', 3), "\n"))
  cat("---\n")
  cat("Please, use summary() function to get additional information\n")
}
