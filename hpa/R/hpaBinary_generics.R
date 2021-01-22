#' Predict method for hpaBinary
#' @param object Object of class "hpaBinary"
#' @template newdata_Template
#' @param is_prob logical; if TRUE (default) then function returns 
#' predicted probabilities. Otherwise latent variable
#' (single index) estimates will be returned.
#' @template elipsis_Template
#' @return This function returns predicted probabilities based on 
#' \code{\link[hpa]{hpaBinary}} estimation results.
predict.hpaBinary <- function (object, ..., 
                              newdata = NULL, 
                              is_prob = TRUE)
{
  if (length(list(...)) > 0)
  {
      warning("Additional arguments passed through ... are ignored.")   
  }
  return(predict_hpaBinary(object, newdata, is_prob))
}
###
#' Summarizing hpaBinary Fits
#' @param object Object of class "hpaBinary"
#' @template elipsis_Template
#' @return This function returns the same list as \code{\link[hpa]{hpaBinary}} 
#' function changing its class to "summary.hpaBinary".
summary.hpaBinary <- function (object, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  return(summary_hpaBinary(object))
}
###
#' Summary for "hpaBinary" object
#' @param x Object of class "hpaBinary"
#' @template elipsis_Template
print.summary.hpaBinary <- function (x, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  return(print_summary_hpaBinary(x))
}
###
#' Plot hpaBinary random errors approximated density
#' @param x Object of class "hpaBinary"
#' @param y this parameter currently ignored
#' @template elipsis_Template
plot.hpaBinary <- function (x, y = NULL, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  if (!is.null(y))
  {
    warning("Note that y parameter currently ignored.")   
  }

  return(plot_hpaBinary(x))
}

###
#' Calculates log-likelihood for "hpaBinary" object
#' @description This function calculates log-likelihood for "hpaBinary" object
#' @usage \method{logLik}{hpaBinary}(object, ...)
#' @param object Object of class "hpaBinary"
#' @template elipsis_Template
logLik.hpaBinary <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  lnL <- logLik_hpaBinary(object)
  attr(lnL, "class") <- "logLik"
  attr(lnL, "df") <- length(as.vector(object$x1))
  
  return(lnL)
}
###
#' Print method for "hpaBinary" object
#' @param x Object of class "hpaBinary"
#' @template elipsis_Template
print.hpaBinary <- function (x, ...) 
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
