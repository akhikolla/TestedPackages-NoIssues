\name{eive.cga}
\alias{eive.cga}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Errors-in-variable estimation in linear regression with compact genetic algorithms
}
\usage{
eive.cga(dirtyx, otherx = NULL, y, numdummies = 10, popsize = 20)
}
\arguments{
  \item{dirtyx}{
Vector of values of independent variable measured with error
}
  \item{otherx}{
Matrix of other independent variables.
}
  \item{y}{
Vector of values of dependent variable
}
  \item{numdummies}{
Number of dummy variables used in algorithm. By default, it is 10.
}
  \item{popsize}{
Population size parameter used in CGA. By default it is 20.
}
}
\details{
Algorithm performs a genetic search to separate mismeasured independent variable into clean and error parts.
}
\value{
 \item{ols}{lm object calculated using original values}
 \item{eive}{lm object calculated using the predicted variable by eive}
 \item{proxy}{lm object of proxy regression obtained by genetic search.}
}
\author{
Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>
Erkin Diyarbakirlioglu <ediyarbakirlioglu@gmail.com>
}
\description{
Method performs a genetic search to find dummy variables that used in a two stage linear regression to reduce errors-in-variables
bias in linear regression.
}

\examples{
# Creating an artificial data

# Loading required package
require ("eive")

# Setting random number generator seed to 12345
# so each time the script runs, same numbers will
# be generated
set.seed(12345)

# Number of observations is set to 30
n<-30

# Unobserved X values are drawn from a Normal distribution
# with mean 10 and variance 7
clean.x <- rnorm(n, mean=10, sd=sqrt(7))

# Measurement error values are dranw from a Normal distribution
# with mean 0 and variance 3
delta.x <- rnorm(n, mean=0,  sd=sqrt(3))

# Error term of regression. Normally distributed with mean 0 and 
# variance 5
e <- rnorm(n, mean=0, sd=sqrt(5))

# Generating Y values using the linear model
# In this model, intercept is 20 and slope is 10. 
y<- 20 + 10* clean.x + e

# Generating observed X values by adding measurement errors
# to unobserved X
dirty.x <- clean.x + delta.x

# Performs a genetic search to find dummy variables that
# used in two stage least squares.
# Please un-comment the line below
# result <- eive.cga (dirtyx=dirty.x, y=y, numdummies=10)

# Print the result
# Please un-comment the line below
# print(result)

########################################### OUTPUT #############################################
# $ols
#
# Call:
# lm(formula = y ~ dirtyx)
# 
# Coefficients:
# (Intercept)       dirtyx  
#     63.590        5.533  
#
#
# $eive
# 
# Call:
# lm(formula = y ~ ols.proxy$fitted.values)
#
# Coefficients:
#            (Intercept)  ols.proxy$fitted.values  
#                 23.863                    9.229  
#
#
# $proxy
#
# Call:
# lm(formula = dirtyx ~ matrix(best, nrow = n))
# 
# Coefficients:
#              (Intercept)   matrix(best, nrow = n)1   matrix(best, nrow = n)2  
#                 12.9321                   -0.6252                   -1.9923  
# matrix(best, nrow = n)3   matrix(best, nrow = n)4   matrix(best, nrow = n)5  
#                  0.7537                   -0.7076                   -0.5247  
# matrix(best, nrow = n)6   matrix(best, nrow = n)7   matrix(best, nrow = n)8  
#                 -0.9196                   -2.0802                   -0.9246  
# matrix(best, nrow = n)9  matrix(best, nrow = n)10  
#                 -0.6164                    1.9694  
######################################### END OF OUTPUT ##########################################

}
