context("testing dann")
library(dann)
library(mlbench)
library(magrittr)
library(dplyr)

###############################################
# Same Resutls as python version?
###############################################
######################
# Problem 1
######################
xTest <- matrix(0, nrow = 5, ncol = 2)
xTrain <- matrix(0, nrow = 5, ncol = 2)

xTrain[, 1] <- c(1, 2, 3, 4, 5)
xTrain[, 2] <- c(6, 7, 8, 9, 10)
yTrain <- c(rep(1, 2), rep(2, 3))

xTest[, 1] <- c(5, 4, 3, 2, 1)
xTest[, 2] <- c(10, 9, 8, 7, 6)
dannPreds <- dann(xTrain, yTrain, xTest, 3, 5, 1, FALSE)

test_that("Validate structure", {
  expect_true(is.vector(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(length(dannPreds) == nrow(xTest))
})

test_that("Compare results to python version. Problem #1", {
  expect_true(all(dannPreds == c(2, 2, 2, 1, 1)))
})

rm(xTest, xTrain, yTrain, dannPreds)

######################
# Problem 2
######################
xTest <- matrix(0, nrow = 10, ncol = 3)
xTrain <- matrix(0, nrow = 10, ncol = 3)

xTrain[, 1] <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
xTrain[, 2] <- c(6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
xTrain[, 3] <- c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
yTrain <- c(rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2), rep(5, 2))

xTest[, 1] <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
xTest[, 2] <- c(7, 7, 8, 8, 10, 10, 12, 12, 14, 14)
xTest[, 3] <- c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
dannPreds <- dann(xTrain, yTrain, xTest, 1, 4, 1, FALSE)

test_that("Validate structure", {
  expect_true(is.vector(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(length(dannPreds) == nrow(xTest))
})

test_that("Compare results to python version. Problem #2", {
  expect_true(all(dannPreds == c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)))
})

rm(xTest, xTrain, yTrain, dannPreds)

###############################################
# Easy problems
###############################################
######################
# Problem 1
######################
set.seed(1)
train <- mlbench.2dnormals(1000, cl = 2, r = sqrt(2), sd = .2) %>%
  tibble::as_tibble()
colnames(train) <- c("X1", "X2", "Y")

xTrain <- train %>%
  select(X1, X2) %>%
  as.matrix()

yTrain <- train %>%
  pull(Y) %>%
  as.numeric() %>%
  as.vector()

test <- mlbench.2dnormals(1000, cl = 2, r = sqrt(2), sd = .2) %>%
  tibble::as_tibble()
colnames(test) <- c("X1", "X2", "Y")

xTest <- test %>%
  select(X1, X2) %>%
  as.matrix()

yTest <- test %>%
  pull(Y) %>%
  as.numeric() %>%
  as.vector()

dannPreds <- dann(xTrain, yTrain, xTest, 5, 50, 1, FALSE)

test_that("Validate structure", {
  expect_true(is.vector(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(length(dannPreds) == nrow(xTest))
})

test_that("Compare predictions to observed #1", {
  expect_true(all(dannPreds == yTest))
})

rm(train, test)
rm(xTrain, yTrain)
rm(xTest, yTest)
rm(dannPreds)

######################
# Problem 2
######################
set.seed(1)
train <- mlbench.hypercube(n = 1000, d = 3, sides = rep(1, 3), sd = 0.1) %>%
  tibble::as_tibble()
colnames(train) <- c("X1", "X2", "X3", "Y")

xTrain <- train %>%
  select(X1, X2, X3) %>%
  as.matrix()

yTrain <- train %>%
  pull(Y) %>%
  as.numeric() %>%
  as.vector()

test <- mlbench.hypercube(n = 1000, d = 3, sides = rep(1, 3), sd = 0.1) %>%
  tibble::as_tibble()
colnames(test) <- c("X1", "X2", "X3", "Y")

xTest <- test %>%
  select(X1, X2, X3) %>%
  as.matrix()

yTest <- test %>%
  pull(Y) %>%
  as.numeric() %>%
  as.vector()

dannPreds <- dann(xTrain, yTrain, xTest, 7, 100, 1, FALSE)

test_that("Validate structure", {
  expect_true(is.vector(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(length(dannPreds) == nrow(xTest))
})

test_that("Compare predictions to observed #2", {
  expect_true(all(dannPreds == yTest))
})

rm(train, test)
rm(xTrain, yTrain)
rm(xTest, yTest)
rm(dannPreds)

######################
# Problem 3
######################
set.seed(1)
train <- mlbench.circle(500, 2) %>%
  tibble::as_tibble()
colnames(train) <- c("X1", "X2", "Y")

xTrain <- train %>%
  select(X1, X2) %>%
  as.matrix()

yTrain <- train %>%
  pull(Y) %>%
  as.numeric() %>%
  as.vector()

test <- mlbench.circle(500, 2) %>%
  tibble::as_tibble()
colnames(test) <- c("X1", "X2", "Y")

xTest <- test %>%
  select(X1, X2) %>%
  as.matrix()

yTest <- test %>%
  pull(Y) %>%
  as.numeric() %>%
  as.vector()

dannPreds <- dann(xTrain, yTrain, xTest, 7, 50, 1, FALSE)

test_that("Validate structure", {
  expect_true(is.vector(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(length(dannPreds) == nrow(xTest))
})

test_that("Compare predictions to observed #3", {
  expect_true(mean(dannPreds == yTest) > .95)
})

rm(train, test)
rm(xTrain, yTrain)
rm(xTest, yTest)
rm(dannPreds)

######################
# Ties are broken by most common class
######################
xTrain <- matrix(c(0, 0, 1), nrow = 3, ncol = 1)

yTrain <- matrix(c(0, 0, 1), nrow = 3, ncol = 1)

# Splits 0 and 1 evenly.
xTest <- matrix(.5, nrow = 1, ncol = 1)

test_that("Run multiple times to confirm results are consitant", {
  expect_true(dann(xTrain, yTrain, xTest, 1, 3, 1, FALSE) == 0)
  expect_true(dann(xTrain, yTrain, xTest, 1, 3, 1, FALSE) == 0)
  expect_true(dann(xTrain, yTrain, xTest, 1, 3, 1, FALSE) == 0)
  expect_true(dann(xTrain, yTrain, xTest, 1, 3, 1, FALSE) == 0)
  expect_true(dann(xTrain, yTrain, xTest, 1, 3, 1, FALSE) == 0)
})

rm(xTrain, yTrain, xTest)

###############################################
# Confirm class numeric value shifting works as expected
###############################################
######################
# Class
######################
xTest <- matrix(0, nrow = 5, ncol = 2)
xTrain <- matrix(0, nrow = 5, ncol = 2)

xTrain[, 1] <- c(1, 2, 3, 4, 5)
xTrain[, 2] <- c(6, 7, 8, 9, 10)
yTrain <- c(rep(-3, 2), rep(-2, 3))

xTest[, 1] <- c(5, 4, 3, 2, 1)
xTest[, 2] <- c(10, 9, 8, 7, 6)
dannPreds <- dann(xTrain, yTrain, xTest, 3, 5, 1, FALSE)

test_that("Validate structure", {
  expect_true(is.vector(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(length(dannPreds) == nrow(xTest))
})

test_that("Shift logic test #1", {
  expect_true(all(dannPreds == c(-2, -2, -2, -3, -3)))
})

# Repeat with non sequental classes
yTrain <- c(rep(-3, 2), rep(2, 3))
dannPreds <- dann(xTrain, yTrain, xTest, 3, 5, 1, FALSE)
test_that("Validate structure", {
  expect_true(is.vector(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(length(dannPreds) == nrow(xTest))
})

test_that("Shift logic test #2", {
  expect_true(all(dannPreds == c(2, 2, 2, -3, -3)))
})

# Repeat normal use case. 0 and 1.
yTrain <- c(rep(0, 2), rep(1, 3))
dannPreds <- dann(xTrain, yTrain, xTest, 3, 5, 1, FALSE)
test_that("Validate structure", {
  expect_true(is.vector(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(length(dannPreds) == nrow(xTest))
})

test_that("Shift logic test #3", {
  expect_true(all(dannPreds == c(1, 1, 1, 0, 0)))
})

######################
# probabilities
######################

yTrain <- c(rep(0, 2), rep(1, 3))
dannPreds <- dann(xTrain, yTrain, xTest, 3, 5, 1, TRUE)
test_that("Validate structure", {
  expect_true(is.matrix(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(nrow(dannPreds) == nrow(xTest))
  expect_true(ncol(dannPreds) == 2)
  expect_true(all(colnames(dannPreds) == c("Class0", "Class1")))
})

yTrain <- c(rep(-4, 2), rep(-1, 3))
dannPreds <- dann(xTrain, yTrain, xTest, 3, 5, 1, TRUE)
test_that("Validate structure", {
  expect_true(is.matrix(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(nrow(dannPreds) == nrow(xTest))
  expect_true(ncol(dannPreds) == 2)
  expect_true(all(colnames(dannPreds) == c("Class-4", "Class-1")))
})

rm(xTest, xTrain, yTrain, dannPreds)

###############################################
# Confirm class probabilities look correct.
###############################################
######################
# Generate problem
######################
set.seed(1)
train <- mlbench.circle(20, 2) %>%
  tibble::as_tibble()
colnames(train) <- c("X1", "X2", "Y")

xTrain <- train %>%
  select(X1, X2) %>%
  as.matrix()

yTrain <- train %>%
  pull(Y) %>%
  as.numeric() %>%
  as.vector()

test <- mlbench.circle(20, 2) %>%
  tibble::as_tibble()
colnames(test) <- c("X1", "X2", "Y")

xTest <- test %>%
  select(X1, X2) %>%
  as.matrix()

yTest <- test %>%
  pull(Y) %>%
  as.numeric() %>%
  as.vector()

######################
# K equal 1
######################
K <- 1
dannPreds <- dann(xTrain, yTrain, xTest, K, 15, 1, TRUE)
test_that("Validate structure", {
  expect_true(is.matrix(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(nrow(dannPreds) == nrow(xTest))
  expect_true(ncol(dannPreds) == 2)
  expect_true(all(colnames(dannPreds) == c("Class1", "Class2")))
})

test_that("Confirm class probabilities sum to 1", {
  expect_true(all(rowSums(dannPreds) == 1))
})

possibleProb <- 0:K / K
test_that("Confirm class probabilities are divisible by k", {
  expect_true(all(dannPreds %in% possibleProb))
})

######################
# K equal 2
######################
K <- 2
dannPreds <- dann(xTrain, yTrain, xTest, K, 15, 1, TRUE)
test_that("Validate structure", {
  expect_true(is.matrix(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(nrow(dannPreds) == nrow(xTest))
  expect_true(ncol(dannPreds) == 2)
  expect_true(all(colnames(dannPreds) == c("Class1", "Class2")))
})

test_that("Confirm class probabilities sum to 1", {
  expect_true(all(rowSums(dannPreds) == 1))
})

possibleProb <- 0:K / K
test_that("Confirm class probabilities are divisible by k", {
  expect_true(all(dannPreds %in% possibleProb))
})

######################
# K equal 3
######################
K <- 3
dannPreds <- dann(xTrain, yTrain, xTest, K, 15, 1, TRUE)
test_that("Validate structure", {
  expect_true(is.matrix(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(nrow(dannPreds) == nrow(xTest))
  expect_true(ncol(dannPreds) == 2)
  expect_true(all(colnames(dannPreds) == c("Class1", "Class2")))
})

test_that("Confirm class probabilities sum to 1", {
  expect_true(all(rowSums(dannPreds) == 1))
})

possibleProb <- 0:K / K
test_that("Confirm class probabilities are divisible by k", {
  expect_true(all(dannPreds %in% possibleProb))
})

######################
# K equal 5
######################
K <- 5
dannPreds <- dann(xTrain, yTrain, xTest, K, 15, 1, TRUE)
test_that("Validate structure", {
  expect_true(is.matrix(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(nrow(dannPreds) == nrow(xTest))
  expect_true(ncol(dannPreds) == 2)
  expect_true(all(colnames(dannPreds) == c("Class1", "Class2")))
})

test_that("Confirm class probabilities sum to 1", {
  expect_true(all(rowSums(dannPreds) == 1))
})

possibleProb <- 0:K / K
test_that("Confirm class probabilities are divisible by k", {
  expect_true(all(dannPreds %in% possibleProb))
})

######################
# K equal 10
######################
K <- 10
dannPreds <- dann(xTrain, yTrain, xTest, K, 15, 1, TRUE)
test_that("Validate structure", {
  expect_true(is.matrix(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(nrow(dannPreds) == nrow(xTest))
  expect_true(ncol(dannPreds) == 2)
  expect_true(all(colnames(dannPreds) == c("Class1", "Class2")))
})

test_that("Confirm class probabilities sum to 1", {
  expect_true(all(rowSums(dannPreds) == 1))
})

possibleProb <- 0:K / K
test_that("Confirm class probabilities are divisible by k", {
  expect_true(all(dannPreds %in% possibleProb))
})

rm(train, test)
rm(xTrain, yTrain)
rm(xTest, yTest)
rm(K)

###############################################
# Create data for checking
###############################################
set.seed(1)
xTest <- matrix(0, nrow = 100, ncol = 2)
xTrain <- matrix(0, nrow = 100, ncol = 2)

xTrain[, 1] <- runif(100, -10, 1)
xTrain[, 2] <- runif(100, -1, 1)
yTrain <- c(rep(1, 50), rep(2, 50))

xTest[, 1] <- runif(100, -1, 1)
xTest[, 2] <- runif(100, -1, 1)

###############################################
# All legitimate values of probability work
###############################################
dannPreds <- dann(xTrain, yTrain, xTest, 2, 5, 1, FALSE)
test_that("Validate structure", {
  expect_true(is.vector(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(length(dannPreds) == nrow(xTest))
  expect_true(all(colnames(dannPreds) == "Class"))
})

dannPreds <- dann(xTrain, yTrain, xTest, 2, 5, 1, TRUE)
test_that("Validate structure", {
  expect_true(is.matrix(dannPreds))
  expect_true(is.numeric(dannPreds))
  expect_true(nrow(dannPreds) == nrow(xTest))
  expect_true(ncol(dannPreds) == 2)
  expect_true(all(colnames(dannPreds) == c("Class1", "Class2")))
})

###############################################
# Input checking
###############################################
#######
# Data checks
#######
chars <- matrix("A", nrow = 5, ncol = 2)
test_that("Nonnumeric inputs error", {
  expect_error(dann(chars, yTrain, xTest, 3, 10, 1), NULL)
  expect_error(dann(xTrain, chars, xTest, 3, 10, 1), NULL)
  expect_error(dann(xTrain, yTrain, chars, 3, 10, 1), NULL)
})
rm(chars)

missingValues <- yTrain
missingValues[1] <- NA
test_that("Missing values in inputs error", {
  expect_error(dann(missingValues, yTrain, xTest), NULL)
  expect_error(dann(xTrain, missingValues, xTest), NULL)
  expect_error(dann(xTrain, yTrain, missingValues), NULL)
})
rm(missingValues)

xTrainrowMissing <- xTrain[1:(nrow(xTrain) - 1), ]
yTrainrowMissing <- yTrain[1:(length(yTrain) - 1)]
test_that("Differnet number of rows in xTrain and yTrain error.", {
  expect_error(dann(xTrainrowMissing, yTrain, xTest), NULL)
  expect_error(dann(xTrain, yTrainrowMissing, xTest), NULL)
})
rm(xTrainrowMissing, yTrainrowMissing)

noDataxTrain <- xTrain[0, ]
noDatayTrain <- yTrain[0]
noDataxTest <- xTest[0, ]
test_that("No rows in inputs error", {
  expect_error(dann(noDataxTrain, noDatayTrain, xTest), NULL)
  expect_error(dann(xTrain, noDatayTrain, noDataxTest), NULL)
})
rm(noDataxTrain, noDatayTrain, noDataxTest)

WrongVarTrain <- xTrain[0, 1:(ncol(xTrain) - 1)]
WrongVarTest <- xTest[0, 1:(ncol(xTest) - 1)]
test_that("Different number of columns in xTrain and xTest error.", {
  expect_error(dann(WrongVarTrain, yTrain, xTest), NULL)
  expect_error(dann(xTrain, yTrain, WrongVarTest), NULL)
})
rm(WrongVarTrain, WrongVarTest)

TooManyyTrain <- cbind(yTrain, yTrain)
test_that("Too many columns in yTrain error.", {
  expect_error(dann(xTrain, TooManyyTrain, xTest), NULL)
})
rm(TooManyyTrain)

#######
# non data checks
#######
test_that("k checks works", {
  expect_error(dann(xTrain, yTrain, xTest, c(3, 2), 3, 1), NULL)
  expect_error(dann(xTrain, yTrain, xTest, "3", 3, 1), NULL)
  expect_error(dann(xTrain, yTrain, xTest, 100000, 3, 1), NULL)
  expect_error(dann(xTrain, yTrain, xTest, 0, 3, 1), NULL)
})

test_that("neighborhood_size checks works", {
  expect_error(dann(xTrain, yTrain, xTest, 2, c(2, 3), 1), NULL)
  expect_error(dann(xTrain, yTrain, xTest, 2, "3", 1), NULL)
  expect_error(dann(xTrain, yTrain, xTest, 2, 100000, 1), NULL)
  expect_error(dann(xTrain, yTrain, xTest, 2, 0, 1), NULL)
  expect_error(dann(xTrain, yTrain, xTest, 4, 3, 1), NULL)
})

test_that("epsilon checks works", {
  expect_error(dann(xTrain, yTrain, xTest, 2, 2, c(2, 3)), NULL)
  expect_error(dann(xTrain, yTrain, xTest, 2, 2, "1"), NULL)
  expect_error(dann(xTrain, yTrain, xTest, 2, 2, -1), NULL)
})

test_that("probability checks works", {
  expect_error(dann(xTrain, yTrain, xTest, 2, 2, 1, c(TRUE, FALSE)), NULL)
  expect_error(dann(xTrain, yTrain, xTest, 2, 2, 1, "TRUE"), NULL)
})
