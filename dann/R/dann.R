#' Discriminant Adaptive Nearest Neighbor Classification
#'
#' @param xTrain Train features. Something easily converted to a numeric matrix.
#' @param yTrain Train classes. Something easily converted to a numeric vector.
#' @param xTest Test features. Something easily converted to a numeric matrix.
#' @param k The number of data points used for final classification.
#' @param neighborhood_size The number of data points used to calculate between and within class covariance.
#' @param epsilon Diagonal elements of a diagonal matrix. 1 is the identity matrix.
#' @param probability Should probabilities instead of classes be returned?
#' @return  A numeric vector containing predicted class or a numeric matrix containing class probabilities.
#' @keywords internal
dann_source <- function(xTrain, yTrain, xTest, k = 5, neighborhood_size = max(floor(nrow(xTrain) / 5), 50), epsilon = 1, probability = FALSE) {
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
    stop("Argument k be at length 1 vector.")
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
    stop("Argument neighborhood_size be at length 1 vector.")
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
    stop("Argument epsilon be at length 1 vector.")
  }
  if (!is.numeric(epsilon)) {
    stop("Argument epsilon should be numeric.")
  }
  if (epsilon < 0) {
    stop("Argument epsilon should be at least 0.")
  }

  # probability is valid
  if (length(probability) != 1) {
    stop("Argument probability be at length 1 vector.")
  }
  if (!is.logical(probability)) {
    stop("Argument probability should be logical.")
  }

  ###################################
  # Shift classes if needed. Need min(yTrain) > 0
  ###################################
  if (min(yTrain) <= 0) {
    shiftedBy <- abs(min(yTrain)) + 1
    yTrain <- yTrain + shiftedBy
    shifted <- TRUE
  } else {
    shifted <- FALSE
  }

  ###################################
  # Calculate predictions
  ###################################

  if (!probability) {
    predictions <- rep(-1, nrow(xTest))
  } else {
    predictions <- matrix(0, nrow = nrow(xTest), ncol = length(unique(yTrain)))
    colnames(predictions) <- stringr::str_c("Class", as.character(sort(unique(yTrain))))
  }

  NCOLX <- ncol(xTrain)

  ###################################
  # Count number of rows per class
  ###################################
  # Used in dann distance sorting
  # If there is a tie in distance, break tie with most common class.
  Y_counts <- vector(mode = "numeric", length = length(unique(yTrain)))
  names(Y_counts) <- sort(unique(yTrain))
  for (i in seq_along(1:length(Y_counts))) {
    Y_counts[i] <- length(yTrain[which(yTrain == names(Y_counts)[i])])
  }
  Y_counts <- sort(Y_counts, decreasing = TRUE)

  Y_class_presidence <- vector(mode = "numeric", length = length(yTrain))
  for (i in seq_along(1:length(Y_counts))) {
    Y_class_presidence[which(yTrain == names(Y_counts)[i])] <- i
  }

  for (i in seq_along(1:nrow(xTest))) {

    ###########
    # Find neighborhood for x[i,]
    ###########
    distances <- calc_distance_C(xTrain, xTest[i, ])

    nearest_neighbors <- order(distances, Y_class_presidence)[1:neighborhood_size]
    neighborhood_xTrain <- xTrain[nearest_neighbors, 1:NCOLX, drop = FALSE]
    neighborhood_X_mean <- colMeans(neighborhood_xTrain)
    neighborhood_y <- yTrain[nearest_neighbors]
    neighborhood_classes <- unique(neighborhood_y)

    ###########
    # Between and within matrices
    ###########
    class_frequencies <- vector(mode = "numeric", length = length(neighborhood_classes))
    within_class_cov <- matrix(0, nrow = NCOLX, ncol = NCOLX)
    between_class_cov <- matrix(0, nrow = NCOLX, ncol = NCOLX)

    for (kth in seq_along(1:length(neighborhood_classes))) {
      target_class <- neighborhood_classes[kth]
      class_indices <- which(neighborhood_y == target_class)
      class_frequencies[target_class] <- sum(neighborhood_y == target_class) / neighborhood_size

      class_covariance <- stats::var(neighborhood_xTrain[class_indices, 1:ncol(neighborhood_xTrain), drop = FALSE])
      # Deal with 1 row in class edge case
      if (all(is.na(class_covariance))) {
        class_covariance <- matrix(0, nrow = nrow(class_covariance), ncol = ncol(class_covariance))
      }

      within_class_cov <- class_covariance * class_frequencies[target_class] + within_class_cov
      class_mean <- colMeans(neighborhood_xTrain[class_indices, 1:ncol(neighborhood_xTrain), drop = FALSE])
      between_class_cov <- outer(class_mean - neighborhood_X_mean, class_mean - neighborhood_X_mean) *
        class_frequencies[target_class] + between_class_cov
    }

    # W* = W^-.5
    # B* = W*BW*
    W_star <- within_class_cov^.5
    # Deal with NA case
    for (kth in seq_along(1:ncol(W_star))) {
      W_star[which(is.na(W_star[, kth])), kth] <- 0
    }
    W_star <- MASS::ginv(W_star)
    B_star <- W_star %*% between_class_cov %*% W_star
    I <- diag(NCOLX)

    sigma <- W_star %*% (B_star + epsilon * I) %*% W_star

    ###########
    # DANN distance using sigma
    ###########
    distances <- vector(mode = "numeric", length = nrow(xTrain))
    for (kth in seq_along(1:length(distances))) {
      distances[kth] <- DANN_distance_C(xTest[i, 1:NCOLX, drop = FALSE], xTrain[kth, 1:NCOLX, drop = FALSE ], sigma)
    }
    nearest <- order(distances, Y_class_presidence)[1:k]
    if (!probability) {
      predictions[i] <- MODE(yTrain[nearest])
    } else {
      predictions[i, ] <- class_proportions(yTrain[nearest], sort(unique(yTrain)))
    }
  }

  ###################################
  # Shift classes back if needed.
  ###################################
  if (shifted & probability) {
    yTrain <- yTrain - shiftedBy
    colnames(predictions) <- stringr::str_c("Class", as.character(sort(unique(yTrain))))
  } else if (shifted & !probability) {
    predictions <- predictions - shiftedBy
  }

  return(predictions)
}

#' Discriminant Adaptive Nearest Neighbor Classification
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
#' @return  A numeric vector containing predicted class or a numeric matrix containing class probabilities.
#' @details
#' This is an implementation of Hastie and Tibshirani's
#' \href{https://web.stanford.edu/~hastie/Papers/dann_IEEE.pdf}{Discriminant Adaptive Nearest
#' Neighbor Classification publication.}.
#' The code is a port of Christopher Jenness's
#' python \href{https://github.com/christopherjenness/ML-lib}{implementation.}
#' @examples
#' library(dann)
#' library(mlbench)
#' library(magrittr)
#' library(dplyr)
#' library(ggplot2)
#' 
#' ######################
#' # Circle Data
#' ######################
#' set.seed(1)
#' train <- mlbench.circle(300, 2) %>%
#'   tibble::as_tibble()
#' colnames(train) <- c("X1", "X2", "Y")
#' 
#' ggplot(train, aes(x = X1, y = X2, colour = Y)) +
#'   geom_point() +
#'   labs(title = "Train Data")
#' 
#' xTrain <- train %>%
#'   select(X1, X2) %>%
#'   as.matrix()
#' 
#' yTrain <- train %>%
#'   pull(Y) %>%
#'   as.numeric() %>%
#'   as.vector()
#' 
#' test <- mlbench.circle(100, 2) %>%
#'   tibble::as_tibble()
#' colnames(test) <- c("X1", "X2", "Y")
#' 
#' ggplot(test, aes(x = X1, y = X2, colour = Y)) +
#'   geom_point() +
#'   labs(title = "Test Data")
#' 
#' xTest <- test %>%
#'   select(X1, X2) %>%
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
#' mean(dannPreds == yTest) # An accurate model.
#' 
#' rm(train, test)
#' rm(xTrain, yTrain)
#' rm(xTest, yTest)
#' rm(dannPreds)
#' @export
dann <- compiler::cmpfun(f = dann_source, options = list(optimize = 3))
