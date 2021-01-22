#' Discriminant Adaptive Nearest Neighbor With Subspace Reduction
#'
#' @param xTrain Train features. Something easily converted to a numeric matrix.
#' @param yTrain Train classes. Something easily converted to a numeric vector.
#' @param xTest Test features. Something easily converted to a numeric matrix.
#' @param k The number of data points used for final classification.
#' @param neighborhood_size The number of data points used to calculate between and within class covariance.
#' @param epsilon Diagonal elements of a diagonal matrix. 1 is the identity matrix.
#' @param probability Should probabilities instead of classes be returned?
#' @param weighted weighted argument to ncoord. See \code{\link[fpc]{ncoord}} for details.
#' @param sphere weighted argument to ncoord. See \code{\link[fpc]{ncoord}} for details.
#' @param numDim Dimension of subspace used by dann. See \code{\link[fpc]{ncoord}} for details.
#' @return  A numeric vector containing predicted class or a numeric matrix containing class probabilities.
#' @keywords internal
sub_dann_source <- function(xTrain, yTrain, xTest,
                            k = 5, neighborhood_size = max(floor(nrow(xTrain) / 5), 50),
                            epsilon = 1, probability = FALSE,
                            weighted = FALSE, sphere = "mcd", numDim = ncol(xTrain) / 2) {
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
  if (!is.matrix(xTest)) {
    xTest <- as.matrix(xTest)
  }

  # Confirm converstion worked
  if (!is.matrix(xTrain)) {
    stop("Was not able to convert argment xTrain to a matrix.")
  }
  if (!is.vector(yTrain)) {
    stop("Was not able to convert argment yTrain to a vector.")
  }
  if (!is.matrix(xTest)) {
    stop("Was not able to convert argment xTest to a matrix.")
  }

  # Confirm numeric
  if (!is.numeric(xTrain)) {
    stop("Argument xTrain should be numeric.")
  }
  if (!is.numeric(yTrain)) {
    stop("Argument yTrain should be numeric.")
  }
  if (!is.numeric(xTest)) {
    stop("Argument xTest should be numeric.")
  }

  # Missing values.
  if (any(is.na(xTrain))) {
    stop("Argument xTrain should not have any NA values.")
  }
  if (any(is.na(yTrain))) {
    stop("Argument yTrain should not have any NA values.")
  }
  if (any(is.na(xTest))) {
    stop("Argument xTest should not have any NA values.")
  }

  # Confirm structure looks right
  if (ncol(xTrain) != ncol(xTest)) {
    stop("Argument xTrain and xTest should have the same number of columns.")
  }
  if (nrow(xTrain) != length(yTrain)) {
    stop("nrow(xTrain) should match length(yTrain).")
  }
  if (ncol(xTrain) < 1) {
    stop("Argument xTrain should have at least one column.")
  }
  if (ncol(xTest) < 1) {
    stop("Argument xTest should have at least one column.")
  }
  if (nrow(xTrain) < 1) {
    stop("Argument xTrain should have at least one row.")
  }
  if (length(yTrain) < 1) {
    stop("Argument yTrain should have positive length.")
  }
  if (nrow(xTest) < 1) {
    stop("Argument xTest should have at least one row.")
  }

  # k is valid
  if (length(k) != 1) {
    stop("Argument k should be a length 1 vector.")
  }
  if (!is.numeric(k)) {
    stop("Argument k should be numeric.")
  }
  if (k > nrow(xTrain)) {
    stop("Argument k should be less than or equal to the numer of rows in xTrain.")
  }
  if (k <= 0) {
    stop("Argument k should be at least 1.")
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
  if (k > neighborhood_size) {
    stop("Argument k should be less than argument neighborhood_size.")
  }

  # epsilon is valid
  if (length(epsilon) != 1) {
    stop("Argument epsilon should be a length 1 vector.")
  }
  if (!is.numeric(epsilon)) {
    stop("Argument epsilon should be numeric.")
  }
  if (epsilon < 0) {
    stop("Argument epsilon should be at least 0.")
  }

  # probability is valid
  if (length(probability) != 1) {
    stop("Argument probability should be a length 1 vector.")
  }
  if (!is.logical(probability)) {
    stop("Argument probability should be logical.")
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

  # numDim is valid
  if (length(numDim) != 1) {
    stop("Argument numDim should be a length 1 vector.")
  }
  if (!is.numeric(numDim)) {
    stop("Argument numDim should be numeric.")
  }
  if (numDim < 1) {
    stop("Argument numDim should be at least 1.")
  }

  # Find subspace
  subspace <- fpc::ncoord(
    xd = xTrain, clvecd = yTrain,
    nn = neighborhood_size, weighted = weighted,
    sphere = "mcd", countmode = 999999999999999
  )

  xTrain2 <- subspace$proj[, 1:numDim, drop = FALSE]
  xTest2 <- xTest %*% subspace$units[, 1:numDim, drop = FALSE]

  # Get predictions
  predictions <- dann(xTrain2, yTrain, xTest2, k, neighborhood_size, epsilon, probability)

  return(predictions)
}

#' Discriminant Adaptive Nearest Neighbor With Subspace Reduction
#'
#' @param xTrain Train features. Something easily converted to a numeric matrix.
#'               Generally columns should have mean zero and standard deviation one beforehand.
#' @param yTrain Train classes. Something easily converted to a numeric vector.
#' @param xTest Test features. Something easily converted to a numeric matrix.
#'              Generally columns should be centered and scaled according to xTrain beforehand.
#' @param k The number of data points used for final classification.
#' @param neighborhood_size The number of data points used to calculate between and within class covariance.
#' @param epsilon Diagonal elements of a diagonal matrix. 1 is the identity matrix.
#' @param probability Should probabilities instead of classes be returned?
#' @param weighted weighted argument to ncoord. See \code{\link[fpc]{ncoord}} for details.
#' @param sphere sphere argument to ncoord. See \code{\link[fpc]{ncoord}} for details.
#' @param numDim Dimension of subspace used by dann. See \code{\link[fpc]{ncoord}} for details.
#' @return  A numeric vector containing predicted class or a numeric matrix containing class probabilities.
#' @details
#' This is an implementation of Hastie and Tibshirani's sub-dann in section 4.1 of
#' \href{https://web.stanford.edu/~hastie/Papers/dann_IEEE.pdf}{Discriminant Adaptive Nearest
#' Neighbor Classification publication.}. It uses package fpc's ncoord to find the subspace. Then calls
#' dann.
#'
#' dann's performance suffers when noise variables are included in the model. Simulations show sub_dann
#' will generally be more performant in this scenario. However there is no replacement for good feature
#' selection.
#' @examples
#' library(dann)
#' library(mlbench)
#' library(magrittr)
#' library(dplyr)
#' library(ggplot2)
#' 
#' ######################
#' # Circle data with unrelated variables
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
#' test <- mlbench.circle(100, 2) %>%
#'   tibble::as_tibble()
#' colnames(test)[1:3] <- c("X1", "X2", "Y")
#' 
#' # Add 5 unrelated variables
#' test <- test %>%
#'   mutate(
#'     U1 = runif(100, -1, 1),
#'     U2 = runif(100, -1, 1),
#'     U3 = runif(100, -1, 1),
#'     U4 = runif(100, -1, 1),
#'     U5 = runif(100, -1, 1)
#'   )
#' 
#' xTest <- test %>%
#'   select(X1, X2, U1, U2, U3, U4, U5) %>%
#'   as.matrix()
#' 
#' yTest <- test %>%
#'   pull(Y) %>%
#'   as.numeric() %>%
#'   as.vector()
#' 
#' dannPreds <- dann(
#'   xTrain = xTrain, yTrain = yTrain, xTest = xTest,
#'   k = 3, neighborhood_size = 50, epsilon = 1,
#'   probability = FALSE
#' )
#' mean(dannPreds == yTest) # Not a good model
#' 
#' # Data suggests a subspace with 2 dimentions. The correct answer.
#' graph_eigenvalues(
#'   xTrain = xTrain, yTrain = yTrain, neighborhood_size = 50,
#'   weighted = FALSE, sphere = "mcd"
#' )
#' 
#' subDannPreds <- sub_dann(
#'   xTrain = xTrain, yTrain = yTrain, xTest = xTest,
#'   k = 3, neighborhood_size = 50, epsilon = 1,
#'   probability = FALSE,
#'   weighted = FALSE, sphere = "classical", numDim = 2
#' )
#' # sub_dan does much better when unrelated variables are present.
#' mean(subDannPreds == yTest)
#' 
#' rm(train, test)
#' rm(xTrain, yTrain)
#' rm(xTest, yTest)
#' rm(dannPreds, subDannPreds)
#' @export
sub_dann <- compiler::cmpfun(f = sub_dann_source, options = list(optimize = 3))
