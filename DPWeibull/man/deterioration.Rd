\name{deterioration}
\alias{deterioration}
\docType{data}
\title{Time to Cosmetic Deterioration of Breast Cancer Patients}
\description{
  We reused the dataset "deterioration" of the time to cosmetic deterioration of the breast for breast
  cancer patients provided by the package DPpackage, which is orphaned in CRAN now. The two treatments
  are radiation alone and radiation coupled with chemotherapy. The event time is either right censored
  or interval censored.
}
\usage{data(deterioration)}
\format{
  A data frame with 94 observations on the following 3 variables.
  \describe{
    \item{\code{left}}{a numeric vector giving the left limit of the interval}
    \item{\code{right}}{a numeric vector giving the right limit of the interval, -999 stands for right censoring}
    \item{\code{trt}}{a numeric vector giving the treatment (0 = radiation only, 
                                                             1 = radiation plus chemotherapy)}
  }
}

\source{
  Beadle, G., Come, S., Henderson, C., Silver, B., and Hellman, S. (1984a). The effect of
  adjuvant chemotherapy on the cosmetic results after primary radiation treatment for early
  stage breast cancer. International Journal of Radiation Oncology, Biology and Physics, 
  10: 2131-2137.
  
  Beadle, G., Harris, J., Silver, B., Botnick, L., and Hellman, S. (1984b). Cosmetic results
  following primary radiation therapy for early breast cancer. Cancer, 54: 2911-2918.
  
  Finkelstein, D.M. and Wolfe, R.A. (1985). A semiparametric model for regression analysis of
  interval-censored failure time data. Biometrics, 41: 933-945.
}

\examples{
  data(deterioration)
  ## maybe str(deterioration) ; plot(deterioration) ...
}
\keyword{datasets}
