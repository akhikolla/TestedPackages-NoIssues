\name{GGClassification-package}
\alias{GGClassification-package}
\alias{GGClassification}
\docType{package}
\title{
  Gabriel Graph Based Large-Margin Classifiers.
}
\description{
  Contains the implementation of a binary large margin classifier based on Gabriel Graph.
}
\details{
  Functions were implemented in C++.
  
  The first function, "GabrielGraph", generates a graph based on Gabriel Graph's construction rule. Edges are linked according to the Euclidean distance among data samples.
  
  The second function, called "model", is used to calculate classifier parameters. It returns a named list with four parameters, the midpoints between opposite-class edges, a vector w containing all weights, a vector representing classifier's bias terms and the corresponding labels of input data.
  
  The third function, named "predict", takes the parameters returned by function "model" in order to classify input data.
  
}
\author{
Thiago Malta Coutinho <thiagomcoutinho@ufmg.br>, Felipe Velloso Campos <felipevellosoc@ufmg.br>, Gustavo Rodrigues Lacerda Silva <gustavolacerdas@gmail.com>, Luiz Carlos Bambirra <luizlitc@gmail.com>, Antonio de Padua Braga <apbraga@ufmg.br> 
}
\references{
Gabriel, K. R., & Sokal, R. R. (1969). A New Statistical Approach to Geographic Variation Analysis. Systematic Zoology, 18(3), 259. doi:10.2307/2412323 

Matula, David & Sokal, Robert. (1980). Properties of Gabriel Graphs Relevant to Geographic Variation Research and the Clustering of Points in the Plane. Geographical Analysis. 12. 205 - 222. 10.1111/j.1538-4632.1980.tb00031.x. 

L. C. B. Torres, C. L. Castro, F. Coelho, F. Sill Torres and A. P. Braga, "Distance-based large margin classifier suitable for integrated circuit implementation," in Electronics Letters, vol. 51, no. 24, pp. 1967-1969, 19 11 2015.
  
} 
\keyword{ package } 
\seealso{
  For more related papers, please check our research group website: 
  http://litc.cpdee.ufmg.br/
}
\examples{
  
    # The example shows a binary classification problem, characterized by two Gaussian classes, 
    # in order to demonstrate how to use the classification functions.
     
    nc = 100
    
    xc1 <- matrix(0.3 * rnorm(nc) + 2.5, ncol = 2)
    xc2 <- matrix(0.3 * rnorm(nc) + 3.5, ncol = 2)
  
    xc1 <- cbind(xc1, rep(0, times = nc/2))
    xc2 <- cbind(xc2, rep(1, times = nc/2))
  
    X <- rbind(xc1, xc2)
    suffled_indexes <- sample(nc)
  
    train_size = nc * 0.7
  
    X_train <- X[suffled_indexes[1:train_size], cbind(1,2)]
    y_train <- X[suffled_indexes[1:train_size], 3]
  
    X_test <- X[suffled_indexes[(71:100)], cbind(1,2)]
    y_test <- X[suffled_indexes[(71:100)], 3]
  
    mdl <- model(X_train, y_train)
    prd <- predict(mdl, X_test)

}
