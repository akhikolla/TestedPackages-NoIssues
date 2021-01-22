#' A helper for sub_dann
#'
#' @param xTrain Train features. Something easily converted to a numeric matrix.
#' @param yTrain Train classes. Something easily converted to a numeric vector.
#' @param neighborhood_size The number of data points used to calculate between and within class covariance.
#' @param weighted weighted argument to ncoord. See \code{\link[fpc]{ncoord}} for details.
#' @param sphere sphere argument to ncoord. See \code{\link[fpc]{ncoord}} for details.
#' @return  A ggplot graph.
#' @details This function plots the eigenvalues found by \code{\link[fpc]{ncoord}}. The user
#' should make a judgement call on how many eigenvalues are large and set sub_dann's
#' numDim to that number.
#' @importFrom rlang .data
#' @examples
#' library(dann)
#' library(mlbench)
#' library(magrittr)
#' library(dplyr)
#' 
#' ######################
#' # Circle data with 2 related variables and 5 unrelated variables
#' ######################
#' set.seed(1)
#' train <- mlbench.circle(300, 2) %>%
#'   tibble::as_tibble()
#' colnames(train)[1:3] <- c("X1", "X2", "Y")
#' 
#' # Add 5 unrelated variables
#' train <- train %>%
#'   mutate(
#'     U1 = runif(300, -1, 1),
#'     U2 = runif(300, -1, 1),
#'     U3 = runif(300, -1, 1),
#'     U4 = runif(300, -1, 1),
#'     U5 = runif(300, -1, 1)
#'   )
#' 
#' xTrain <- train %>%
#'   select(X1, X2, U1, U2, U3, U4, U5) %>%
#'   as.matrix()
#' 
#' yTrain <- train %>%
#'   pull(Y) %>%
#'   as.numeric() %>%
#'   as.vector()
#' 
#' # Data suggests a subspace with 2 dimentions. The correct answer.
#' graph_eigenvalues(
#'   xTrain = xTrain, yTrain = yTrain,
#'   neighborhood_size = 50, weighted = FALSE, sphere = "mcd"
#' )
#' 
#' 
#' rm(train)
#' rm(xTrain, yTrain)
#' @export
graph_eigenvalues <- function(xTrain, yTrain,
                              neighborhood_size = max(floor(nrow(xTrain) / 5), 50),
                              weighted = FALSE, sphere = "mcd") {
  ###################################
  # Input checking
  ###################################
  # Convert to matrices
  if (!is.matrix(xTrain)) {
    xTrain <- as.matrix(xTrain)
  }
  if (!is.vector(yTrain)) {
    yTrain <- as.vector(yTrain)
  }

  # Confirm converstion worked
  if (!is.matrix(xTrain)) {
    stop("Was not able to convert argment xTrain to a matrix.")
  }
  if (!is.vector(yTrain)) {
    stop("Was not able to convert argment yTrain to a vector.")
  }

  # Confirm numeric
  if (!is.numeric(xTrain)) {
    stop("Argument xTrain should be numeric.")
  }
  if (!is.numeric(yTrain)) {
    stop("Argument yTrain should be numeric.")
  }

  # Missing values.
  if (any(is.na(xTrain))) {
    stop("Argument xTrain should not have any NA values.")
  }
  if (any(is.na(yTrain))) {
    stop("Argument yTrain should not have any NA values.")
  }

  # Confirm structure looks right
  if (nrow(xTrain) != length(yTrain)) {
    stop("nrow(xTrain) should match length(yTrain).")
  }
  if (ncol(xTrain) < 1) {
    stop("Argument xTrain should have at least one column.")
  }
  if (nrow(xTrain) < 1) {
    stop("Argument xTrain should have at least one row.")
  }
  if (length(yTrain) < 1) {
    stop("Argument yTrain should have positive length.")
  }

  # neighborhood_size is valid
  if (length(neighborhood_size) != 1) {
    stop("Argument neighborhood_size should be a length 1 vector.")
  }
  if (!is.numeric(neighborhood_size)) {
    stop("Argument neighborhood_size should be numeric.")
  }
  if (neighborhood_size > nrow(xTrain)) {
    stop("Argument neighborhood_size should be less than or equal to the numer of rows in xTrain.")
  }
  if (neighborhood_size <= 1) {
    stop("Argument neighborhood_size should be at least 2.")
  }

  # weighted is valid
  if (length(weighted) != 1) {
    stop("Argument weighted should be a length 1 vector.")
  }
  if (!is.logical(weighted)) {
    stop("Argument weighted should be logical.")
  }

  # sphere is valid
  if (length(sphere) != 1) {
    stop("Argument sphere should be a length 1 vector.")
  }
  if (!is.character(sphere)) {
    stop("Argument sphere should be a character vector.")
  }
  if (!(sphere %in% c("mve", "mcd", "classical"))) {
    stop("Argument sphere should be a one mve, mcd, or classical.")
  }

  # Find subspace
  subspace <- fpc::ncoord(
    xd = xTrain, clvecd = yTrain,
    nn = neighborhood_size, weighted = weighted,
    sphere = "mcd", countmode = 999999999999999
  )

  eigen <- tibble::enframe(subspace$ev, value = "eigenValues", name = "order")

  graph <- ggplot2::ggplot(eigen, ggplot2::aes(x = .data$order, y = .data$eigenValues)) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = 1:nrow(eigen)) +
    ggplot2::labs(x = "Rank Order", y = "Eigenvalues")

  return(graph)
}
