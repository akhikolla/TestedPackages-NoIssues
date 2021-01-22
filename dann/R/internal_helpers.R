###################################
# Helper functions
###################################
#' Computes the distance between x0 and x1 using the DANN metric
#'
#' @param x1 A numeric matrix with training predictors as columns.
#' @param x2 A numeric matrix with training predictors as columns.
#' @param sigma A numeric matrix defined in Hastie's DANN publication.
#' @keywords internal
DANN_distance <- function(x0, x1, sigma) {
  difference <- x0 - x1
  distance <- difference %*% sigma %*% t(difference)
  return(distance)
}

#' Computes mode.
#' Code found at \href{https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode}{Stack Overflow}
#'
#' @param x A numeric vector.
#' @keywords internal
MODE <- function(x) {
  ux <- sort(unique(x))
  return(ux[which.max(tabulate(match(x, ux)))])
}

#' Computes class proportions
#'
#' @param x A numeric vector.
#' @param possibleValues A vector of all possible values x can contain
#' @keywords internal
class_proportions <- function(x, possibleValues) {
  counts <- purrr::map_dbl(possibleValues, function(i) {
    temp <- x[which(x == i)]
    length(temp)
  })

  out <- counts / sum(counts)
  return(out)
}

#' @keywords internal
.onUnload <- function(libpath) {
  library.dynam.unload("dann", libpath)
}
