if(!require(parallel)){
  install.packages("parallel", repos = "http://cran.us.r-project.org")
  library(parallel)
}
if(!require(minerva)){
  install.packages("minerva", repos = "http://cran.us.r-project.org")
  library(minerva)
}
if(!require(Hmisc)){
  install.packages("Hmisc", repos = "http://cran.us.r-project.org")
  library(Hmisc)
}

#'@importFrom Rcpp evalCpp
#'@useDynLib testforDEP

#' @title Test dependence for two data
#' @description This function computes test statistic, p value, and confidence interval for dependence based on classic methods: Pearson, Kendall, Spearman, and modern methods: Vexler, Kallenberg, MIC, Hoeffding, and Empirical Likelihood tests.
#' @param x a numeric vector stores first variable.
#' @param y numeric vector stores second variable.
#' @param data (Optional) a data frame stores data to be tested.
#' @param test a character indicating which test to implement.. Must be one of \{"PEARSON", "KENDALL", "SPEARMAN", "VEXLER", "TS2", "V", "MIC", "HOEFFD", "EL"\}
#' @param p.opt a character specifying p value to be obtained by distribution or by Monte Carlo simulation. Must be "dist", "MC" or "table".
#' @param num.MC a numeric for number of Monte Carlo simulations.
#' @param BS.CI a numeric specifying alpha for Bootstrap confidence interval. When equal 0, confidence interval won't be computed.
#' @param rm.na a TRUE/ FALSE flag indicating whether remove missing data (NA) in input.
#' @param set.seed a TRUE/ FALSE flag indicating whether set seed for Monte Carlo simulation and bootstrap sampling.
#'
#'
#' @author Jeffrey C. Miecznikowski, En-shuo Hsu, Yanhua Chen, Albert Vexler
#' @details
#'Argument "x, y" and "data" are two different ways to input data. When x or y is missing, data will be taken as input; while x, y and data all exist leads to error.
#'Argument data is a two-column numeric data frame. The order of columns does not affect results.
#'Since modern test methods: "VEXLER", "TS2", "V", "MIC", "HOEFFD", and "EL" have no continuous probability density function, argument p.opt = "dist" does not apply. For classic methods, when p.opt is "dist", argument num.MC will be ignored.
#'p.opt = "table" use interpolation from pre stored simulated tables. Current version only supports "VEXLER", "MIC", "HOEFFD" and "EL" tests.
#'For Vexler, MIC and EL, since computation is more time-consuming, a warning with estimated execution time will be returned when input size > 100. Input size <= 100 is recommanded for Monte Carlo p-value. For input size > 100 use table.
#'num.MC should be a integer between 100 and 10,000 for acceptable computation times.
#'NA in input is not acceptable. Set rm.na = TRUE to remove.
#'More details see \link{Pearson}, \link{Kendall}, \link{Spearman}, \link{Vexler}, \link{Kallenberg}, \link{MIC}, \link{Hoeffding}, \link{EL}.
#'
#' @export
#' @import parallel
#' @import minerva
#' @import Hmisc
#' @import stats
#' @import graphics
#' @import methods
#'
#' @return
#' an S4 object of class "testforDEP_result", having attributes: test statistics (TS), p value (p_value) and confidence interval (CI) if apply.
#'
#' @examples
#' set.seed(123)
#' x = runif(100, 0, 1)
#' y = runif(100, 0, 1)
#'
#' testforDEP(x, y, test = "SPEARMAN", p.opt = "MC",
#'            num.MC = 10000, BS.CI = 0, set.seed = TRUE)
#'
#'
#' #An object of class "testforDEP_result"
#' #Slot "TS":
#' #[1] 59.54311
#'
#' #Slot "p_value":
#' #[1] 0.6735326
#'
#' #Slot "CI":
#' #list()
#'
#' @seealso
#' Technical report:
#' http://sphhp.buffalo.edu/content/dam/sphhp/biostatistics/Documents/techreports/UB-Biostatistics-TR1701.pdf


testforDEP = function(x = NA, y = NA, data = NA, test, p.opt = "MC", num.MC = 10000, BS.CI = 0, rm.na = FALSE, set.seed = FALSE){

  #Check test is valid
  if(!(test == "PEARSON" || test == "KENDALL" || test == "SPEARMAN" ||
    test == "VEXLER" || test == "TS2" || test == "V" || test == "MIC" || test == "HOEFFD"
   || test == "EL")){
    stop("Please specify a test.")
  }

  #Check input is valid
  if(!is.numeric(x) || !is.numeric(y)){
    if(!is.data.frame(data))
      stop("Please input paired vectors x, y. Or input data frame data.")
    else
      warning("Using data as input. Ignoring x, y.")
  }
  else if(is.numeric(x) && is.numeric(y) && is.data.frame(data)){
      stop("Cannot take both x, y pairs and data.")
  }
  else{
    if(length(x) != length(y))
      stop("x and y must have same length.")
    data = data.frame(cbind(x, y))
  }

  #Check NA and Remove NA
  if(anyNA(data)){
    if(rm.na)
      data = data[!is.na(data[,1]) & !is.na(data[,2]),]
     else
       stop("Missing data (NA) detected. Please handle or use rm.na = True to remove.")
  }


  pData = pointer(data, environment())
  obj = new(test, pdata = pData, p.opt = p.opt, num.MC = num.MC, BS.CI = BS.CI, set.seed = set.seed)
  delete("pData", environment())

  return(test(obj))

}
