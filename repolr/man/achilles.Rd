\name{achilles}
\alias{achilles}
\docType{data}
\title{Achilles Tendon Rupture}
\description{
EuroQol activity scores, measured on a three point ordinal score scale, post-surgery for 48 patients at three occasions. 
}
\usage{data(achilles)}
\format{
A data frame with 144 observations on the following four variables.
  \describe{
    \item{\code{Patient}}{a patient identifier variable.}
    \item{\code{Treat}}{post-surgery treatments are either immediate mobilisation in a carbon-fibre orthosis with three 1.5cm heel raises (1) or traditional plaster cast immobilisation (2).}
    \item{\code{Time}}{recorded at baseline (1), six months (2) and one year (3) post-surgery.}
    \item{\code{Activity}}{ability to undertake usual activities post-surgery; this was scored by each patient as either no problem (1), some problem (2) or an inability (3) to perform usual activity (e.g. work, leisure, housework etc).}
  }
}
\source{
Costa, M.L., MacMillan, K., Halliday, D., Chester, R., Shepstone, L., Robinson, A.H.N., Donell, S.T. (2006). Randomised controlled trials of immediate weight-bearing mobilisation for rupture of the tendon Achillis. \emph{Journal of Bone and Joint Surgery (British)} 88-B, 69-77.
}
\examples{
data(achilles)
head(achilles)
}
\keyword{datasets}
