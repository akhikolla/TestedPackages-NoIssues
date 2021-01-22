\name{lowerbound}
\alias{lowerbound}
\alias{lowerbound_tube}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
lowerbound
}
\description{
Calculate the lowerbound for the DTW distance measure in linear time.}
\usage{
lowerbound(C, ws, scale = c("z", "01", "none"),
           dist_method = c("norm1", "norm2", "norm2_square"), 
           Q = NULL, tube = NULL)
                     
lowerbound_tube(Q, ws, scale = c("z", "01", "none"))
                     
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{
vector or matrix, the query time series
}
  \item{C}{
vector or matrix, the query time series
}
  \item{dist_method}{
distance method, one of ("norm1", "norm2", "norm2_square")
}

  \item{scale}{
either "none", so no scaling is performed, or one of ("z", "01") to scale both Q and C. Also see \code{\link{dtw}}
}
  \item{ws}{
see \code{\link{dtw}}
}
\item{tube}{tube for lower bounding. "tube" can be the output from lowerbound_tube(). If tube = NULL, then Q must not be NULL, so that tube can be defined. If the tube is passed as argument to lowerbound(), then it is necessary that the scale parameter in the lowerbound() call is identical to the scaling method applied on Q before calculating the tube.
}

}


\value{
lowerbound distance measure that is proven to be smaller than the DTW distance measure
}


\details{

Lower Bounding: The following methods are implemented:
\itemize{
\item \code{LB_Keogh} for univariate time series (Keogh et al. 2005)
\item \code{LB_MV} for multivariate time series with the \code{dist_method = "norm2_square"}, (Rath et al. 2002)
\item Adjusted for different distance methods "norm1" and "norm2", inspired by (Rath et al. 2002).
}
}


\references{
\itemize{
\item Keogh, Eamonn, and Chotirat Ann Ratanamahatana. "Exact indexing of dynamic time warping." Knowledge and information systems 7.3 (2005): 358-386.
\item Rath, Toni M., and R. Manmatha. "Lower-bounding of dynamic time warping distances for multivariate time series." University of Massachusetts Amherst Technical Report MM 40 (2002).
\item Sakurai, Yasushi, Christos Faloutsos, and Masashi Yamamuro. "Stream monitoring under the time warping distance." Data Engineering, 2007. ICDE 2007. IEEE 23rd International Conference on. IEEE, 2007.
}

}


\examples{

\dontrun{

#--- Univariate time series Q and C
ws <- sample(2:40, size = 1)
dist_method <- "norm1"
N <- 50
N <- 50
Q <- cumsum(rnorm(N))
C <- cumsum(rnorm(N))
Q.z <- IncDTW::scale(Q, "z")
C.z <- IncDTW::scale(C, "z")

lb.z <- lowerbound(C = C.z, ws = ws, scale ="none", dist_method = dist_method, Q = Q.z)
lb <- lowerbound(C = C, ws = ws, scale ="z", dist_method = dist_method, Q = Q)
d1 <- dtw2vec(Q = Q.z, C = C.z, step_pattern = "symmetric1", 
  dist_method = dist_method, ws = ws)$distance
d2 <- dtw2vec(Q = Q.z, C = C.z, step_pattern = "symmetric2", 
  dist_method = dist_method, ws = ws)$distance

c(lb, lb.z, d1, d2)


#--- with pre-calculated tube
ws <- sample(2:40, size = 1)
dist_method <- "norm1"
N <- 50
N <- 50
Q <- cumsum(rnorm(N))
C <- cumsum(rnorm(N))
Q.z <- IncDTW::scale(Q, "z")
C.z <- IncDTW::scale(C, "z")

tube <- lowerbound_tube(Q, ws, scale = "z")

lb.z <- lowerbound(C = C.z, ws = ws, scale ="none", dist_method = dist_method, tube = tube)
lb <- lowerbound(C = C, ws = ws, scale ="z", dist_method = dist_method, tube = tube)
d1 <- dtw2vec(Q = Q.z, C = C.z, step_pattern = "symmetric1", 
  dist_method = dist_method, ws = ws)$distance
d2 <- dtw2vec(Q = Q.z, C = C.z, step_pattern = "symmetric2", 
  dist_method = dist_method, ws = ws)$distance

c(lb, lb.z, d1, d2)



#--- Multivariate time series Q and C
ws <- sample(2:40, size = 1)
dist_method <- sample(c("norm1", "norm2", "norm2_square"), size = 1)
N <- 50
Q <- matrix(cumsum(rnorm(N * 3)), ncol = 3)
C <- matrix(cumsum(rnorm(N * 3)), ncol = 3)
Q.z <- IncDTW::scale(Q, "z")
C.z <- IncDTW::scale(C, "z")

lb.z <- lowerbound(C = C.z, ws = ws, scale ="none", dist_method = dist_method, Q = Q.z)
lb <- lowerbound(C = C, ws = ws, scale ="z", dist_method = dist_method, Q = Q)
d1 <- dtw2vec(Q = Q.z, C = C.z, step_pattern = "symmetric1", 
  dist_method = dist_method, ws = ws)$distance
d2 <- dtw2vec(Q = Q.z, C = C.z, step_pattern = "symmetric2", 
  dist_method = dist_method, ws = ws)$distance

c(lb, lb.z, d1, d2)

}
}
