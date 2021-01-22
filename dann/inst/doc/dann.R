## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----simData-------------------------------------------------------------
library(dann)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(mlbench)

set.seed(1)
train <- mlbench.2dnormals(600, cl = 6, r = sqrt(2), sd = .5) %>%
  tibble::as_tibble()
colnames(train) <- c("X1", "X2", "Y")

ggplot(train, aes(x = X1, y = X2, colour = Y)) + 
  geom_point() + 
  labs(title = "Train Data")


test <- mlbench.2dnormals(600, cl = 6, r = sqrt(2), sd = .5) %>%
  tibble::as_tibble()
colnames(test) <- c("X1", "X2", "Y")
ggplot(test, aes(x = X1, y = X2, colour = Y)) + 
  geom_point() + 
  labs(title = "Test Data")

## ----Cluster-------------------------------------------------------------
xTrain <- train %>%
  select(X1, X2) %>%
  as.matrix()
yTrain <- train %>%
  pull(Y) %>%
  as.numeric() %>%
  as.vector()

xTest <- test %>%
  select(X1, X2) %>%
  as.matrix()
yTest <- test %>%
  pull(Y) %>%
  as.numeric() %>%
  as.vector()

## ----model---------------------------------------------------------------
dannPreds <- dann(xTrain = xTrain, yTrain = yTrain, xTest = xTest,
                  k = 7, neighborhood_size = 150, epsilon = 1)
round(mean(dannPreds == yTest), 2)

