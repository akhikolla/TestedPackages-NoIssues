\name{find_peaks}
\alias{find_peaks}
\title{find_peaks}
\usage{
find_peaks(x, w, get_min = TRUE, strict = TRUE)
}
\arguments{
\item{x}{vector}

\item{w}{window, at least w-many values need to be in-between two consecutive
peaks to find both, otherwise only the bigger one is returned}

\item{get_min}{logical (default TRUE) if TRUE, then minima are returned, else maxima}

\item{strict}{logical, if TRUE (default) then a local minimum needs to be smaller then all neighbors. 
If FALSE, then a local minimum needs to be smaller or equal all neighbors.}
}
\value{
integer vector of indices where x has local extreme values
}
\description{
Find negative or positive peaks of a vector in a predefined neighborhood w
}
\examples{

#--- Find the peaks (local minima and maxima), 
# and also the border peak at index 29. First the local maxima:
x <- c(1:10, 9:1, 2:11)
peak_indices <- find_peaks(x, w=3, get_min=FALSE)
peak_indices
x[peak_indices]

# and now the local minima
peak_indices <- find_peaks(x, w=3, get_min=TRUE)
peak_indices
x[peak_indices]


#--- What exactly does the neigbohood parameter 'w' mean?
# At least w-many values need to be inbetween two consecutive peaks:
x <- -c(1:10, 9, 9, 11, 9:8, 7)
peak_indices <- find_peaks(x, w=3)
peak_indices
x[peak_indices]

x <- -c(1:10, 9, 9,9, 11, 9:8, 7)
peak_indices <- find_peaks(x, w=3)
peak_indices
x[peak_indices]


#--- What does the parameter 'strict' mean?
# If strict = TRUE, then the peak must be '<' (or '>') 
# then the neighbors, other wise '<=' (or '>=')
x <- c(10:1,  1:10)
peak_indices <-  find_peaks(x, w=3, strict = TRUE)
peak_indices
x[peak_indices]

peak_indices <-  find_peaks(x, w=3, strict = FALSE)
peak_indices
x[peak_indices]
   
}

