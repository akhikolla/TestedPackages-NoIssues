\name{predict}
\alias{predict}
\title{Inference of classes based on received parameters of large-margin classifier.}
\usage{
predict(model, X)
}
\arguments{
  \item{model}{Large-margin classifier parameters computed using GGClassification::model function.}
  
  \item{X}{Data matrix to be classified.}
  
}
\value{

A vector of labels returned by the classifier.

}
\description{

The function receives a model containing the parameters of the classifier and a data matrix to be classified
Using the parameters and data, the classes/labels of input data is infered and returned as a vector.

}
\examples{

X <- matrix(0.3 * rnorm(100) + 2.5, ncol = 2)
y <- c(rep(0, times = 25), rep(1, times = 25))

mdl <- model(X, y)

x_prd <- matrix(0.3 * rnorm(100) + 4.5, ncol=2)

y_hat <- predict(mdl, x_prd)

}
