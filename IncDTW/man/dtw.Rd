\name{dtw}
\alias{dtw}
\alias{cm}
\alias{print.idtw}
\alias{summary.idtw}
\alias{is.idtw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dynamic Time Warping
}
\description{
Calculate the DTW distance, cost matrices and direction matrices including the warping path two multivariate time series.
}
\usage{
dtw(Q, C, dist_method = c("norm1", "norm2", "norm2_square"), 
    step_pattern = c("symmetric2", "symmetric1"), ws = NULL,
    return_cm = FALSE,
    return_diffM = FALSE,
    return_wp = FALSE,
    return_diffp = FALSE,
    return_QC = FALSE)
                     
cm(Q, C, dist_method = c("norm1", "norm2", "norm2_square"), 
   ws = NULL, ...)
   
\method{print}{idtw}(x, digits = getOption("digits"), ...)

\method{summary}{idtw}(object, ...)

is.idtw(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{Query time series. Q needs to be one of the following: (1) a one dimensional vector, (2) a matrix where each row is one observations and each column is one dimension of the time series, or (3) a matrix of differences/ costs (diffM, cm). If Q and C are matrices they need to have the same number of columns.
}
\item{C}{
Candidate time series. C needs to be one of the following: (1) a one dimensional vector, (2) a matrix where each row is one observations and each column is one dimension of the time series, or (3) if Q is a matrix of differences or costs C needs to be the respective character string 'diffM' or 'cm'.
}
\item{dist_method}{character, describes the method of distance measure for multivariate time series (this parameter is ignored for univariate time series). Currently supported methods are 'norm1' (default, is the Manhattan distance), 'norm2' (is the Euclidean distance) and 'norm2_square'. For the function \code{cm()} the parameter \code{dist_method} can also be a user defined distance function (see details and examples).
}
\item{step_pattern}{character, describes the step pattern. Currently implemented are the patterns \code{symmetric1} and \code{symmetric2}, see details.
}
\item{ws}{
integer, describes the window size for the sakoe chiba window. If NULL, then no window is applied. (default = NULL)
}
\item{return_cm}{logical, if TRUE then the Matrix of costs (the absolute value) is returned. (default = FALSE)
}
\item{return_diffM}{
logical, if TRUE then the Matrix of differences (not the absolute value) is returned. (default = FALSE)
}
\item{return_wp}{
logical, if TRUE then the warping path is returned. (default = FALSE) If return_diffp == TRUE, then return_wp is set to TRUE as well.
}
\item{return_diffp}{
logical, if TRUE then the path of differences (not the absolute value) is returned. (default = FALSE)
}
\item{return_QC}{
logical, if TRUE then the input vectors Q and C are appended to the returned list. This is useful for the \code{\link{plot.idtw}} function.  (default = FALSE)
}
\item{x}{
output from \code{dtw} or \code{\link{idtw}}.
}
\item{object}{
any R object
}
\item{...}{
additional arguments, e.g. passed to \code{print}, \code{summary}, or a user defined distance function for \code{cm()} 
}
\item{digits}{
passed to \code{round} and \code{print}
}

}


\details{
The dynamic time warping distance is the element in the last row and last column of the global cost matrix. 

For the multivariate case where Q is a matrix of n rows and k columns and C is a matrix of m rows and k columns
the dist_method parameter defines the following distance measures: 

 norm1: \deqn{dist(Q_{i,.}, C_{j,.}) = \sum {l = 1:k}  |Q_{i,l} - C_{j,l}|}

 norm2: \deqn{dist(Q_{i,.}, C_{j,.}) = (\sum {l = 1:k} (Q_{i,l} - C_{j,l})^2)^0.5}
 
 norm2_square: \deqn{dist(Q_{i,.}, C_{j,.}) = \sum{l = 1:k} (Q_{i,l} - C_{j,l})^2}
 

The parameter \code{step_pattern} describes how the two time series are aligned. 
If \code{step_pattern == "symmetric1"} then  \deqn{gcm_{i,j} = cm{i,j} + min(gcm_{i-1,j}, gcm{i-1, j-1}, gcm{i, i-1} }.

If \code{step_pattern == "symmetric2"} then  \deqn{gcm_{i,j} = cm{i,j} + min(gcm_{i-1,j}, cm{i,j}+ gcm{i-1, j-1}, gcm{i, i-1} }.
 

To make DTW distances comparable for many time series of different lengths use the \code{normlized_distance} with the setting \code{step_method = 'symmetric2'}. Please find a more detailed discussion and further references here: \url{http://dtw.r-forge.r-project.org/}.


User defined distance function:
To calculate the DTW distance measure of two time series a distance function for the local distance of two observations \code{Q[i, ]} and \code{C[j, ]} of the time series \code{Q} and \code{C} has to be selected. The predefined distance function are 'norm1', 'norm2' and 'norm2-square'. It is also possible to define a customized distance function and use the cost matrix \code{cm} as input for the DTW algorithm, also for the incremental functions. In the following experiment we apply the cosine distance as local distance function:

\deqn{
 d_cos(C_i, Q_j) = 1 - (\sum{o=1:O} Q_{io} * C_{jo})/ ((\sum{o=1:O} Q_{io}^2)^0.5 * (\sum{o=1:O} C_{jo}^2)^0.5).
}





}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{distance}{the DTW distance, that is the element of the last row and last column of gcm}
\item{normalized_distance}{the normalized DTW distance, that is the distance divided by \code{N+M}, where N and M are the lengths of the time series Q and C, respectively. If \code{step_pattern == 'symmetric1'} no normalization is performed and \code{NA} is returned (see details).}
\item{gcm}{global cost matrix}
\item{dm}{direction matrix (3=up, 1=diagonal, 2=left)}
\item{wp}{warping path}
\item{ii}{indices of Q of the optimal path}
\item{jj}{indices of C of the optimal path}
\item{cm}{Matrix of costs}
\item{diffM}{Matrix of differences}
\item{diffp}{path of differences}
\item{Q}{input Q}
\item{C}{input C}
}
\references{
Sakoe, H.; Chiba, S., Dynamic programming algorithm optimization for spoken word recognition, Acoustics, Speech, and Signal Processing [see also IEEE Transactions on Signal Processing], IEEE Transactions on , vol.26, no.1, pp. 43-49, Feb 1978. http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055
}


\examples{

#--- univariate
Q <- cumsum(rnorm(100))
C <- Q[11:100] + rnorm(90, 0, 0.5)
tmp <- dtw(Q = Q, C = C, ws = 15, return_diffM = FALSE, 
           return_QC = TRUE, return_wp = TRUE)
names(tmp)
print(tmp, digits = 3)
plot(tmp)
plot(tmp, type = "warp")


#--- compare different input variations
dtw_base  <- dtw(Q = Q, C = C, ws = 15, return_diffM = TRUE)
dtw_diffM <- dtw(Q = dtw_base$diffM, C = "diffM", ws = 15, 
                 return_diffM = TRUE)
dtw_cm    <- dtw(Q = abs(dtw_base$diffM), C = "cm", ws = 15, 
                 return_diffM = TRUE)

identical(dtw_base$gcm, dtw_cm$gcm)
identical(dtw_base$gcm, dtw_diffM$gcm)

# of course no diffM is returned in the 'cm'-case
dtw_cm$diffM


#--- multivariate case
Q <- matrix(rnorm(100), ncol=2)
C <- matrix(rnorm(80), ncol=2)
dtw(Q = Q, C = C, ws = 30, dist_method = "norm2")


#--- user defined distance function
# We define the distance function d_cos and use it as input for the cost matrix function cm. 
# We can pass the output from cm() to dtw2vec(), and also to idtw2vec() for the 
# incrermental calculation:

d_cos <- function(x, y){ 
  1 - sum(x * y)/(sqrt(sum(x^2)) * sqrt(sum(y^2)))
}

Q <- matrix(rnorm(100), ncol=5, nrow=20)
C <- matrix(rnorm(150), ncol=5, nrow=30)
cm1 <- cm(Q, C, dist_method = d_cos)
dtw2vec(Q = cm1, C = "cm")$distance

res0 <- idtw2vec(Q = cm1[ ,1:20], newObs =  "cm")
idtw2vec(Q = cm1[ ,21:30], newObs =  "cm", gcm_lc = res0$gcm_lc_new)$distance

# The DTW distances -- based on the customized distance function -- of the 
# incremental calculation and the one from scratch are identical.


}
\keyword{ cluster }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line

