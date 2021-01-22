\name{generate.eive.data}
\alias{generate.eive.data}
\title{
Generates data for errors-in-variables model
}
\description{
This function generates data using a linear regression model with size of n. Then one of the independent variables is contaminated
by adding measurement errors. Another independent variable can be included in model. 
}
\usage{
generate.eive.data(n, e.sd, delta.sd, seed = 12345, useotherx = FALSE)
}
\arguments{
  \item{n}{
Number of observations.
}
  \item{e.sd}{
Standard deviation of error term of regression.
}
  \item{delta.sd}{
Standard deviation of measurement error. 
}
  \item{seed}{
Random number seed. By default, it is 12345.
}
  \item{useotherx}{
Boolean variable. If it is TRUE, another variable will be created with no errors. By default, it is FALSE.
}
}
\value{
Returns a matrix of contaminated variable, other variable (if exists) and independent variable in its columns.
}
\author{
Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>
Erkin Diyarbakirlioglu <ediyarbakirlioglu@gmail.com>
}


