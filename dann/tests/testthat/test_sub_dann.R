context("testing sub_dann")
library(dann)
library(mlbench)
library(magrittr)
library(dplyr)

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

subDannPreds <- sub_dann(xTrain, yTrain, xTest, 5, 50, 1, FALSE, FALSE, "mcd", 2)

test_that("Validate structure", {
  expect_true(is.vector(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(length(subDannPreds) == nrow(xTest))
})

test_that("Compare predictions to observed #1", {
  expect_true(all(subDannPreds == yTest))
})

rm(train, test)
rm(xTrain, yTrain)
rm(xTest, yTest)
rm(subDannPreds)

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

subDannPreds <- sub_dann(xTrain, yTrain, xTest, 1, 50, 1, FALSE, FALSE, "mcd", 3)

test_that("Validate structure", {
  expect_true(is.vector(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(length(subDannPreds) == nrow(xTest))
})

test_that("Compare predictions to observed #2", {
  expect_true(mean(subDannPreds == yTest) > .98)
})

rm(train, test)
rm(xTrain, yTrain)
rm(xTest, yTest)
rm(subDannPreds)

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
subDannPreds <- sub_dann(xTrain, yTrain, xTest, K, 15, 1, TRUE, FALSE, "mcd", 2)
test_that("Validate structure", {
  expect_true(is.matrix(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(nrow(subDannPreds) == nrow(xTest))
  expect_true(ncol(subDannPreds) == 2)
  expect_true(all(colnames(subDannPreds) == c("Class1", "Class2")))
})

test_that("Confirm class probabilities sum to 1", {
  expect_true(all(rowSums(subDannPreds) == 1))
})

possibleProb <- 0:K / K
test_that("Confirm class probabilities are divisible by k", {
  expect_true(all(subDannPreds %in% possibleProb))
})

######################
# K equal 2
######################
K <- 2
subDannPreds <- sub_dann(xTrain, yTrain, xTest, K, 15, 1, TRUE, FALSE, "mcd", 2)
test_that("Validate structure", {
  expect_true(is.matrix(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(nrow(subDannPreds) == nrow(xTest))
  expect_true(ncol(subDannPreds) == 2)
  expect_true(all(colnames(subDannPreds) == c("Class1", "Class2")))
})

test_that("Confirm class probabilities sum to 1", {
  expect_true(all(rowSums(subDannPreds) == 1))
})

possibleProb <- 0:K / K
test_that("Confirm class probabilities are divisible by k", {
  expect_true(all(subDannPreds %in% possibleProb))
})

######################
# K equal 3
######################
K <- 3
subDannPreds <- sub_dann(xTrain, yTrain, xTest, K, 15, 1, TRUE, FALSE, "mcd", 2)
test_that("Validate structure", {
  expect_true(is.matrix(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(nrow(subDannPreds) == nrow(xTest))
  expect_true(ncol(subDannPreds) == 2)
  expect_true(all(colnames(subDannPreds) == c("Class1", "Class2")))
})

test_that("Confirm class probabilities sum to 1", {
  expect_true(all(rowSums(subDannPreds) == 1))
})

possibleProb <- 0:K / K
test_that("Confirm class probabilities are divisible by k", {
  expect_true(all(subDannPreds %in% possibleProb))
})

######################
# K equal 5
######################
K <- 5
subDannPreds <- sub_dann(xTrain, yTrain, xTest, K, 15, 1, TRUE, FALSE, "mcd", 2)
test_that("Validate structure", {
  expect_true(is.matrix(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(nrow(subDannPreds) == nrow(xTest))
  expect_true(ncol(subDannPreds) == 2)
  expect_true(all(colnames(subDannPreds) == c("Class1", "Class2")))
})

test_that("Confirm class probabilities sum to 1", {
  expect_true(all(rowSums(subDannPreds) == 1))
})

possibleProb <- 0:K / K
test_that("Confirm class probabilities are divisible by k", {
  expect_true(all(subDannPreds %in% possibleProb))
})

######################
# K equal 10
######################
K <- 10
subDannPreds <- sub_dann(xTrain, yTrain, xTest, K, 15, 1, TRUE, FALSE, "mcd", 2)
test_that("Validate structure", {
  expect_true(is.matrix(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(nrow(subDannPreds) == nrow(xTest))
  expect_true(ncol(subDannPreds) == 2)
  expect_true(all(colnames(subDannPreds) == c("Class1", "Class2")))
})

test_that("Confirm class probabilities sum to 1", {
  expect_true(all(rowSums(subDannPreds) == 1))
})

possibleProb <- 0:K / K
test_that("Confirm class probabilities are divisible by k", {
  expect_true(all(subDannPreds %in% possibleProb))
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
subDannPreds <- sub_dann(xTrain, yTrain, xTest, 2, 5, 1, FALSE, FALSE, "mcd", 2)
test_that("Validate structure", {
  expect_true(is.vector(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(length(subDannPreds) == nrow(xTest))
  expect_true(all(colnames(subDannPreds) == "Class"))
})

subDannPreds <- sub_dann(xTrain, yTrain, xTest, 2, 5, 1, TRUE, FALSE, "mcd", 2)
test_that("Validate structure", {
  expect_true(is.matrix(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(nrow(subDannPreds) == nrow(xTest))
  expect_true(ncol(subDannPreds) == 2)
  expect_true(all(colnames(subDannPreds) == c("Class1", "Class2")))
})

###############################################
# All legitimate values of weighted work
###############################################
subDannPreds <- sub_dann(xTrain, yTrain, xTest, 2, 5, 1, FALSE, FALSE, "mcd", 2)
test_that("Validate structure", {
  expect_true(is.vector(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(length(subDannPreds) == nrow(xTest))
  expect_true(all(colnames(subDannPreds) == "Class"))
})

subDannPreds <- sub_dann(xTrain, yTrain, xTest, 2, 5, 1, FALSE, TRUE, "mcd", 2)
test_that("Validate structure", {
  expect_true(is.vector(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(length(subDannPreds) == nrow(xTest))
  expect_true(all(colnames(subDannPreds) == "Class"))
})

###############################################
# All legitimate values of sphere work
###############################################
subDannPreds <- sub_dann(xTrain, yTrain, xTest, 2, 50, 1, FALSE, FALSE, "mve", 2)
test_that("Validate structure", {
  expect_true(is.vector(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(length(subDannPreds) == nrow(xTest))
  expect_true(all(colnames(subDannPreds) == "Class"))
})

subDannPreds <- sub_dann(xTrain, yTrain, xTest, 2, 50, 1, FALSE, FALSE, "mcd", 2)
test_that("Validate structure", {
  expect_true(is.vector(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(length(subDannPreds) == nrow(xTest))
  expect_true(all(colnames(subDannPreds) == "Class"))
})

subDannPreds <- sub_dann(xTrain, yTrain, xTest, 2, 50, 1, FALSE, FALSE, "classical", 2)
test_that("Validate structure", {
  expect_true(is.vector(subDannPreds))
  expect_true(is.numeric(subDannPreds))
  expect_true(length(subDannPreds) == nrow(xTest))
  expect_true(all(colnames(subDannPreds) == "Class"))
})

###############################################
# Input checking
###############################################
#######
# Data checks
#######
chars <- matrix("A", nrow = 5, ncol = 2)
test_that("Nonnumeric inputs error", {
  expect_error(sub_dann(chars, yTrain, xTest, 3, 10, 1), NULL)
  expect_error(sub_dann(xTrain, chars, xTest, 3, 10, 1), NULL)
  expect_error(sub_dann(xTrain, yTrain, chars, 3, 10, 1), NULL)
})
rm(chars)

missingValues <- yTrain
missingValues[1] <- NA
test_that("Missing values in inputs error", {
  expect_error(sub_dann(missingValues, yTrain, xTest), NULL)
  expect_error(sub_dann(xTrain, missingValues, xTest), NULL)
  expect_error(sub_dann(xTrain, yTrain, missingValues), NULL)
})
rm(missingValues)

xTrainrowMissing <- xTrain[1:(nrow(xTrain) - 1), ]
yTrainrowMissing <- yTrain[1:(length(yTrain) - 1)]
test_that("Differnet number of rows in xTrain and yTrain error.", {
  expect_error(sub_dann(xTrainrowMissing, yTrain, xTest), NULL)
  expect_error(sub_dann(xTrain, yTrainrowMissing, xTest), NULL)
})
rm(xTrainrowMissing, yTrainrowMissing)

noDataxTrain <- xTrain[0, ]
noDatayTrain <- yTrain[0]
noDataxTest <- xTest[0, ]
test_that("No rows in inputs error", {
  expect_error(sub_dann(noDataxTrain, noDatayTrain, xTest), NULL)
  expect_error(sub_dann(xTrain, yTrain, noDataxTest), NULL)
})
rm(noDataxTrain, noDatayTrain, noDataxTest)

WrongVarTrain <- xTrain[0, 1:(ncol(xTrain) - 1)]
WrongVarTest <- xTest[0, 1:(ncol(xTest) - 1)]
test_that("Different number of columns in xTrain and xTest error.", {
  expect_error(sub_dann(WrongVarTrain, yTrain, xTest), NULL)
  expect_error(sub_dann(xTrain, yTrain, WrongVarTest), NULL)
})
rm(WrongVarTrain, WrongVarTest)

TooManyyTrain <- cbind(yTrain, yTrain)
test_that("Too many columns in yTrain error.", {
  expect_error(sub_dann(xTrain, TooManyyTrain, xTest), NULL)
})
rm(TooManyyTrain)

#######
# non data checks
#######
test_that("k checks works", {
  expect_error(sub_dann(xTrain, yTrain, xTest, c(3, 2), 3, 1, FALSE, FALSE, "mcd", 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, "3", 3, 1, FALSE, FALSE, "mcd", 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 100000, 3, 1, FALSE, FALSE, "mcd", 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 0, 3, 1, FALSE, FALSE, "mcd", 2), NULL)
})

test_that("neighborhood_size checks works", {
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, c(2, 3), 1, FALSE, FALSE, "mcd", 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, "3", 1, FALSE, FALSE, "mcd", 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 100000, 1, FALSE, FALSE, "mcd", 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 0, 1, FALSE, FALSE, "mcd", 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 4, 3, 1, FALSE, FALSE, "mcd", 2), NULL)
})

test_that("epsilon checks works", {
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, c(2, 3), FALSE, FALSE, "mcd", 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, "1", FALSE, FALSE, "mcd", 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, -1, FALSE, FALSE, "mcd", 2), NULL)
})

test_that("probability checks works", {
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, 1, c(TRUE, FALSE), FALSE, "mcd", 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, 1, "TRUE", FALSE, "mcd", 2), NULL)
})

test_that("weighted checks works", {
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, 1, FALSE, c(TRUE, FALSE), "mcd", 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, 1, FALSE, "TRUE", "mcd", 2), NULL)
})

test_that("sphere checks works", {
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, 1, FALSE, FALSE, c("mcd", "mcd"), 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, 1, FALSE, FALSE, FALSE, 2), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, 1, FALSE, FALSE, "foo", 2), NULL)
})

test_that("numDim checks works", {
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, 1, FALSE, FALSE, "mcd", c(1, 2)), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, 1, FALSE, FALSE, "mcd", "2"), NULL)
  expect_error(sub_dann(xTrain, yTrain, xTest, 2, 2, 1, FALSE, FALSE, "mcd", 0), NULL)
})
