\name{rundtw}
\alias{rundtw}
\alias{print.rundtw}
\alias{summary.rundtw}
\alias{is.rundtw}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rundtw
}
\description{
Detect recurring patterns similar to given query pattern by measuring the distance with DTW. A window of the length of the query pattern slides along the longer time series and calculates computation-time-efficiently the DTW distance for each point of time. The function incrementally updates the scaling of the sliding window, incrementally updates the cost matrix, applies the vector-based implementation of the DTW algorithm, early abandons and applies lower bounding methods to decrease the calculation time.
}
\usage{
rundtw(Q, C, dist_method = c("norm1", "norm2", "norm2_square"), 
      step_pattern = c("symmetric1", "symmetric2"), k = NULL, 
      scale = c("01", "z", "none"), ws = NULL, threshold = NULL, 
      lower_bound = TRUE, overlap_tol = 0, return_QC = FALSE,
      normalize = c("01", "z", "none"))
      
\method{print}{rundtw}(x, ...)

\method{summary}{rundtw}(object, ...)

is.rundtw(x) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{
vector or matrix, the query time series
}
  \item{C}{
vector or matrix (equal number of columns as Q), the longer time series which is scanned for multiple fits of the query time series. C can also be a list of time series (either all vectors, or all matrices of equal number of columns) with varying lengths. See Details.
}
  \item{dist_method}{
see \code{\link{dtw}}
}
  \item{step_pattern}{
see \code{\link{dtw}}, only for "symmetric1" the lower bounding is implemented yet
}
  \item{k}{
integer >= 0. If \code{k > 0}, then the k-nearest neighbors to the query pattern that are found in all possible sub-sequences of the long time series C are returned. Per default the found fits don't overlap, except the \code{overlap_tol} parameter is adjusted (this should be done with care!). If \code{k > 0} then \code{lowerbound} is set to \code{TRUE}.
}
  \item{scale}{
character, one of c("01", "z", "none") (default = "01"), if not "none" then \code{Q} (once at the start) and \code{C} (running scaling) are scaled. Either min-max ("01") or the z-scaling ("z")  is applied. TRUE (identical to '01') and FALSE (identical to 'none') are deprecated and will be dropped in the next package version.
}
  \item{ws}{
see \code{\link{dtw}}
}
  \item{threshold}{
numeric >= 0, global threshold for early abandoning DTW calculation if this threshold is hit. (also see \code{\link{dtw}}). If NULL (default) no early abandoning is applied.  
}
\item{lower_bound}{
logical, (default = TRUE) If TRUE (default) then lower bounding is applied (see Details).
}
\item{overlap_tol}{
integer between 0 and length of Q, (default = 0) gives the number of observations that two consecutive fits are accepted to overlap. 
}
\item{return_QC}{logical, default = FALSE. If TRUE then \code{Q} and \code{C} are appended to the return list.
}
\item{normalize}{deprecated, use \code{scale} instead. If \code{normalize} is set, then \code{scale} is overwritten by \code{normalize} for compatibility.
}
\item{x}{the output object from \code{\link{rundtw}}.
}
\item{object}{any R object
}
\item{...}{further arguments passed to print or summary.
}

}


\value{
\item{dist}{vector of DTW distances}
\item{counter}{named vector of counters. Gives information how the algorithm proceeded. see Details}
\item{knn_indices}{indices of the found kNN}
\item{knn_values}{DTW distances of the found kNN}
\item{knn_list_indices}{indices of list entries of C, where to find the kNN. Only returned if C is a list of time series. See examples.}
\item{Q, C}{input time series}
}


\details{
This function and algorithm was inspired by the work of Sakurai et al. (2007) and refined for running min-max scaling and lower bounding.

Lower Bounding: The following methods are implemented:
\itemize{
\item \code{LB_Keogh} for univariate time series (Keogh et al. 2005)
\item \code{LB_MV} for multivariate time series with the \code{dist_method = "norm2_square"}, (Rath et al. 2002)
\item Adjusted for different distance methods "norm1" and "norm2", inspired by (Rath et al. 2002).
}


Counter vector:
\itemize{
\item "scale_reset" counts how many times the min and max of the sliding window and the scaling need to be reset completely
\item "scale_new_extreme" how many times the min or max of the sliding window are adjusted incrementally and the scaling need to be reset completely
\item "scale_1step" how many times only the new observation in the sliding window needs to be scaled based on the current min and max
\item "cm_reset" how many times the cost matrix for the sliding window needs to be recalculated completely
\item "cm_1step" how many times only the front running column of the cost matrix is calculated
\item "early_abandon" how many times the early abandon method aborts the DTW calculation before finishing
\item "lower_bound" how many times the lower bounding stops the initialization of the DTW calculation
\item "completed" for how many subsequences the DTW calculation finished
}



C is a list of time series:
If C is a list of time series, the algorithm concatenates the list to one long time series to apply the logic of early abandoning, lower bounding, and finding the kNN. Finally the results are split to match the input. The 
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
#--- Simulate a query pattern Q and a longer time series C with
# distorted instances of Q within C. Apply rundtw() do detect 
# these instances of C. 

rw <- function(nr) cumsum(rnorm(nr))
noise <- function(nr) rnorm(nr)
set.seed(1234)
nC <- 500
nQ <- 40
nfits <- 5

nn <- nC - nfits * nQ # length of noise
nn <- nn/nfits + 1

Q <- sin(seq(from = 1, to = 4 * pi, length.out = nQ))
Qscale <- IncDTW::scale(Q, type = "01")
C <- rw(0)
for(i in 1:nfits){
   C <- c(C, rw(nn)  , 
          Q * abs(rnorm(1, 10, 10)) + 
            rnorm(1, 0, 10) + noise(nQ))
}

# Apply running min-max scaling and allow lower 
# bounding to find the 3 NN
x <- rundtw(Q, C, scale = '01', ws = 10, k = 3, 
            lower_bound = TRUE, return_QC = TRUE)

# Have a look at the result and get the best indices 
# with lowest distance.
x
summary(x)
find_peaks(x$dist, nQ)
plot(x, scale = "01")

# The fourth and fifth simuated fits are not returned, 
# since the DTW distances are higher than the other found 3 NN. 
# The algorithm early abandons and returns NA for these 
# indices. Get all distances by the following command:
x_all <- rundtw(Q, C, scale = '01', ws = 10, 
                k = 0, lower_bound = FALSE)
plot(x_all)

# Do min-max-scaling and lower bound
rundtw(Q, C, scale = '01', ws = 10, lower_bound = TRUE)

# Do z-scaling and lower bound
rundtw(Q, C, scale = 'z', ws = 10, lower_bound = TRUE)

# Don't scaling and don't lower bound
rundtw(Q, C, scale = 'none', ws = 10, lower_bound = FALSE)

# kNN: Do z-scaling and lower bound 
rundtw(Q, C, scale = 'z', ws = 10, k = 3)



#--- For multivariate time series
rw <- function(nr, nco) {
        matrix(cumsum(rnorm(nr * nco)), nrow = nr, ncol = nco)
}

nC <- 500
nQ <- 50
nco <- 2
nfits <- 5

nn <- nC - nfits * nQ# length of noise
nn <- nn/nfits

Q <- rw(nQ, nco)
Qscale <- IncDTW::scale(Q, type="01")
C <- matrix(numeric(), ncol=nco)
for(i in 1:nfits){
   C <- rbind(C, rw(nn, nco), Q)
}

# Do min-max-scaling and lower bound
rundtw(Q, C, scale = '01', ws = 10, threshold = Inf, 
       lower_bound = TRUE)

# Do z-scaling and lower bound
rundtw(Q, C, scale = 'z', ws = 10, threshold = NULL, 
       lower_bound = TRUE)

# Don't scale and don't lower bound
rundtw(Q, C, scale = 'none', ws = 10, threshold = NULL, 
       lower_bound = FALSE)



#--- C can also be a list of (multivariate) time series. 
# So rundtw() detects the closest fits of a query pattern
# across all time series in C.
nC <- 500
nQ <- 50
nco <- 2
rw <- function(nr, nco){
        matrix(cumsum(rnorm(nr * nco)), nrow = nr, ncol = nco)
}

Q <- rw(nQ, nco)
C1 <- rbind(rw(100, nco), Q, rw(20, nco))
C2 <- rbind(rw(10, nco), Q, rw(50, nco))
C3 <- rbind(rw(200, nco), Q, rw(30, nco))
C_list <- list(C1, C2, C3)

# Do min-max-scaling and lower bound
x <- rundtw(Q, C_list, scale = '01', ws = 10, threshold = Inf, 
            lower_bound = TRUE, k = 3, return_QC = TRUE)
x
# Plot the kNN fit of the 2nd or 3rd list entry of C
plot(x, lix = 2)
plot(x, lix = 3)

# Do z-scaling and lower bound
rundtw(Q, C_list, scale = 'z', ws = 10, threshold = Inf, 
       lower_bound = TRUE, k = 3)

# Don't scale and don't lower bound
x <- rundtw(Q, C_list, scale = 'none', ws = 10, 
            lower_bound = FALSE, k = 0, return_QC = TRUE)
x
plot(x)
}
}
