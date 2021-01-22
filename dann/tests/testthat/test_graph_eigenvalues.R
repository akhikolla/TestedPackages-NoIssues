context("testing graph_eigenvalues")
library(dann)
library(mlbench)
library(magrittr)
library(dplyr)

######################
# Circle data with 2 related variables and 5 unrelated variables
######################
set.seed(1)
train <- mlbench.circle(500, 2) %>%
  tibble::as_tibble()
colnames(train)[1:3] <- c("X1", "X2", "Y")

# Add 5 unrelated variables
train <- train %>%
  mutate(
    U1 = runif(500, -1, 1),
    U2 = runif(500, -1, 1),
    U3 = runif(500, -1, 1),
    U4 = runif(500, -1, 1),
    U5 = runif(500, -1, 1)
  )

xTrain <- train %>%
  select(X1, X2, U1, U2, U3, U4, U5) %>%
  as.matrix()

yTrain <- train %>%
  pull(Y) %>%
  as.numeric() %>%
  as.vector()

# Data suggests a subspace with 2 dimentions. The correct answer.
graph <- graph_eigenvalues(xTrain, yTrain, 50, FALSE, "mcd")

test_that("Validate structure", {
  expect_true(all(class(graph) == c("gg", "ggplot")))
})

rm(graph)

###############################################
# All legitimate values of weighted work
###############################################

test_that("Validate structure", {
  expect_true(all(class(graph_eigenvalues(xTrain, yTrain, 50, FALSE, "mcd")) == c("gg", "ggplot")))
})

test_that("Validate structure", {
  expect_true(all(class(graph_eigenvalues(xTrain, yTrain, 50, TRUE, "mcd")) == c("gg", "ggplot")))
})


###############################################
# All legitimate values of sphere work
###############################################
test_that("Validate structure", {
  expect_true(all(class(graph_eigenvalues(xTrain, yTrain, 50, FALSE, "mve")) == c("gg", "ggplot")))
})

test_that("Validate structure", {
  expect_true(all(class(graph_eigenvalues(xTrain, yTrain, 50, TRUE, "mcd")) == c("gg", "ggplot")))
})

test_that("Validate structure", {
  expect_true(all(class(graph_eigenvalues(xTrain, yTrain, 50, TRUE, "classical")) == c("gg", "ggplot")))
})

###############################################
# Input checking
###############################################
#######
# Data checks
#######
chars <- matrix("A", nrow = 5, ncol = 2)
test_that("Nonnumeric inputs error", {
  expect_error(graph_eigenvalues(chars, yTrain, 50, FALSE, "mcd"), NULL)
  expect_error(graph_eigenvalues(xTrain, chars, 50, FALSE, "mcd"), NULL)
})
rm(chars)

missingValues <- yTrain
missingValues[1] <- NA
test_that("Nonnumeric inputs error", {
  expect_error(graph_eigenvalues(missingValues, yTrain, 50, FALSE, "mcd"), NULL)
  expect_error(graph_eigenvalues(xTrain, missingValues, 50, FALSE, "mcd"), NULL)
})
rm(missingValues)

xTrainrowMissing <- xTrain[1:(nrow(xTrain) - 1), ]
yTrainrowMissing <- yTrain[1:(length(yTrain) - 1)]
test_that("Differnet number of rows in xTrain and yTrain error.", {
  expect_error(graph_eigenvalues(xTrainrowMissing, yTrain), NULL)
  expect_error(graph_eigenvalues(xTrain, yTrainrowMissing), NULL)
})
rm(xTrainrowMissing, yTrainrowMissing)

noDataxTrain <- xTrain[0, ]
noDatayTrain <- yTrain[0]
test_that("No rows in inputs error", {
  expect_error(graph_eigenvalues(noDataxTrain, noDatayTrain), NULL)
})
rm(noDataxTrain, noDatayTrain)

#######
# non data checks
#######
test_that("neighborhood_size checks works", {
  expect_error(graph_eigenvalues(xTrain, yTrain, c(2, 3), FALSE, "mcd"), NULL)
  expect_error(graph_eigenvalues(xTrain, yTrain, "3", FALSE, "mcd"), NULL)
  expect_error(graph_eigenvalues(xTrain, yTrain, 100000, FALSE, "mcd"), NULL)
  expect_error(graph_eigenvalues(xTrain, yTrain, 0, FALSE, "mcd"), NULL)
})

test_that("weighted checks works", {
  expect_error(graph_eigenvalues(xTrain, yTrain, 3, c(TRUE, FALSE), "mcd"), NULL)
  expect_error(graph_eigenvalues(xTrain, yTrain, 3, "FALSE", "mcd"), NULL)
})

test_that("sphere checks works", {
  expect_error(graph_eigenvalues(xTrain, yTrain, 3, FALSE, c("mcd", "mcd")), NULL)
  expect_error(graph_eigenvalues(xTrain, yTrain, 3, FALSE, FALSE), NULL)
  expect_error(graph_eigenvalues(xTrain, yTrain, 3, FALSE, "foo"), NULL)
})

rm(train)
rm(xTrain, yTrain)
