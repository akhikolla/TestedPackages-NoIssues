## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----exampleP1-----------------------------------------------------------
 library(dann)
 library(mlbench)
 library(magrittr)
 library(dplyr, warn.conflicts = FALSE)
 library(ggplot2)

 ######################
 # Circle data with unrelated variables
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
 
 test <- mlbench.circle(500, 2) %>%
   tibble::as_tibble()
 colnames(test)[1:3] <- c("X1", "X2", "Y")

 # Add 5 unrelated variables
 test <- test %>%
   mutate(
     U1 = runif(500, -1, 1),
     U2 = runif(500, -1, 1),
     U3 = runif(500, -1, 1),
     U4 = runif(500, -1, 1),
     U5 = runif(500, -1, 1)
   )
 

## ----shape---------------------------------------------------------------
 xTrain <- train %>%
   select(X1, X2, U1, U2, U3, U4, U5) %>%
   as.matrix()

 yTrain <- train %>%
   pull(Y) %>%
   as.numeric() %>%
   as.vector()

 xTest <- test %>%
   select(X1, X2, U1, U2, U3, U4, U5) %>%
   as.matrix()

 yTest <- test %>%
   pull(Y) %>%
   as.numeric() %>%
   as.vector()

## ----dann----------------------------------------------------------------
 dannPreds <- dann(xTrain = xTrain, yTrain = yTrain, xTest = xTest, 
                   k = 3, neighborhood_size = 50, epsilon = 1, probability = FALSE)
 mean(dannPreds == yTest) # Not a good model

## ----graph---------------------------------------------------------------
 
 graph_eigenvalues(xTrain = xTrain, yTrain = yTrain, 
                   neighborhood_size = 50, weighted = FALSE, sphere = "mcd")

## ----subDann-------------------------------------------------------------
 subDannPreds <- sub_dann(xTrain = xTrain, yTrain = yTrain, xTest = xTest, 
                          k = 3, neighborhood_size = 50, epsilon = 1, 
                          probability = FALSE, 
                          weighted = FALSE, sphere = "mcd", numDim = 2)
 mean(subDannPreds == yTest) # sub_dan does much better when unrelated variables are present.

## ----dann2---------------------------------------------------------------
 variableSelectionDann <- dann(xTrain = xTrain[, 1:2], yTrain = yTrain, xTest = xTest[, 1:2],
                               k = 3, neighborhood_size = 50, epsilon = 1, probability = FALSE)
 
 mean(variableSelectionDann == yTest) # Best model found when only true predictors are used.

