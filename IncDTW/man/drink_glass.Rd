\name{drink_glass}
\alias{drink_glass}
\alias{walk}
\alias{brush_teeth}
\docType{data}
\title{
Accelerometer: drink a glass, walk, brush teeth.
}
\description{
3-dimensional acceleration time series recorded during the activities of walking, drinking a glass or brushing teeth.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("drink_glass")}

\format{A list of matrices, where each matrix has 3 columns (x, y, and z axis of the accelerometer). The number of rows differ. 
  
}
\details{
list of 3-dimensional time series stored as matrix. The data is recorded with 32Hz. The data is z-scaled (z-normalized).

}
\source{
UCI Machine Learning Repository
\url{https://archive.ics.uci.edu/ml/datasets/Dataset+for+ADL+Recognition+with+Wrist-worn+Accelerometer}
}

\examples{
\dontrun{
data(drink_glass)
class(drink_glass)
length(drink_glass)
dim(drink_glass[[1]])
matplot(drink_glass[[1]], type="l")

data(walk)
class(walk)
length(walk)
dim(walk[[1]])
matplot(walk[[1]], type="l")

data(brush_teeth)
class(brush_teeth)
length(brush_teeth)
dim(brush_teeth[[1]])
matplot(brush_teeth[[1]], type="l")
}
}
\keyword{datasets}
