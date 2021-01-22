#' @title
#' Trimmed L Moments
#' @description
#' Calculates empirical Trimmed L-moments of specific order(s) and trimming.
#' @param x numeric vector of data.
#' @param order integer, order of Trimmed Moment, has to be greater than 1.
#' @param leftrim,rightrim integer indicating trimming parameters, have to be greater than 0.
#' @param na.rm logical, indicates if NAs should be removed.
#' @param computation.method character, indicating if the computation is performed via
#' PWMs, direct, recursive, or recurrence (see References Hosking & Balakrishnan, 2015).
#' Possible values are \code{auto} (default, automatically choose appropriate method), \code{pwm},
#' \code{direct}, \code{recursive}, or \code{recurrence}. Only if empirical moments are calculated.
#'
#' @return numeric vector, empirical TL(\code{leftrim},\code{rightrim})-moments of orders \code{order} of \code{x}
#'
#' @references Elamir, E. A., & Seheult, A. H. (2003). Trimmed L-moments. Computational Statistics & Data Analysis, 43(3), 299-314.
#' @references Hosking, J. R. (1990). L-moments: analysis and estimation of distributions using linear combinations of order statistics. Journal of the Royal Statistical Society. Series B (Methodological), 105-124.
#' @references Hosking, J. R. M. (2007). Some theory and practical uses of trimmed L-moments. Journal of Statistical Planning and Inference, 137(9), 3024-3039.
#' @references Hosking, J. R. M., & Balakrishnan, N. (2015). A uniqueness result for L-estimators, with applications to L-moments. Statistical Methodology, 24, 69-80.
#' @seealso \code{\link{TLMoments}}
#'
#' @examples
#' x <- rnorm(100)
#' TLMoment(x, order = 1)
#' TLMoment(x, order = 2, leftrim = 0, rightrim = 1)
#' TLMoment(x, order = c(1, 2, 3), leftrim = 2, rightrim = 2)
#' TLMoment(x, order = c(1, 3, 2), leftrim = 2, rightrim = 2)
#' @export
TLMoment <- function(x, order = 1L, leftrim = 0L, rightrim = 0L, na.rm = FALSE, computation.method = "auto") {
  # What should happen with NAs?
  if (na.rm) {
    x <- na.exclude(x)
  }
  # At least as much data points as maximal order
  if (max(order) > length(x)) stop("order higher than number of data points!")
  # Error if invalid computation.method is chosen
  if (!(computation.method %in% c("auto", "pwm", "direct", "recursive", "recurrence"))) stop("computation.method must be either auto, pwm, direct, recursive, or recurrence")

  # Check if computation.method is possible
  if (computation.method == "auto")
    computation.method <- select_computation(leftrim, rightrim)
  if (computation.method %in% c("pwm", "recurrence") & (leftrim != as.integer(leftrim) | rightrim != as.integer(rightrim))) {
    stop("Non-integer trimmings cannot be calculated with computation.method \"pwm\" or \"recurrence\"! Use \"direct\" or \"recursive\" as computation.method")
  }

  ### Forwarding to C++-functions
  if (computation.method %in% c("pwm", "direct")) {
    out <- vapply(order, function(r) {
      switch(computation.method,
             pwm = TLMoment_PWM(x, r, leftrim, rightrim),
             direct = TLMoment_direct(x, r, leftrim, rightrim))
   }, numeric(1))
  } else if (computation.method %in% c("recursive", "recurrence")) {
    out <- switch(computation.method,
                  recursive = TLMoments_recursive(x, max(order), leftrim, rightrim),
                  recurrence = TLMoments_recurrence(x, max(order), leftrim, rightrim))
    out <- out[order]
  }

  # Return results
  setNames(out, paste0("L", order))
}
