\name{unemployment}
\alias{unemployment}
\docType{data}
\title{Unemployment data}
\description{
  Unemployment data about some European countries in 2011.
}
\usage{data(unemployment)}
\format{
 A data.frame with 32 rows and 3 columns.
}
\details{
The source is  Eurostat news-release 104/2012 - 4 July 2012.
The 32 observations are European countries: BELGIUM, BULGARIA, CZECHREPUBLIC, DENMARK, GERMANY, ESTONIA, IRELAND, GREECE, SPAIN, FRANCE, ITALY, CYPRUS, LATVIA, LITHUANIA, LUXEMBOURG, HUNGARY, MALTA, NETHERLANDS, AUSTRIA, POLAND, PORTUGAL, ROMANIA, SLOVENIA, SLOVAKIA, FINLAND, SWEDEN, UNITEDKINGDOM, ICELAND, NORWAY, SWITZERLAND, CROATIA, TURKEY.
The 3 variables are: the total unemployment rate, defined as the percentage of unemployed persons aged 15-74 in the economically active population (Variable 1); the youth unemployment rate, defined as the unemployment rate for young people aged between 15 and 24 (Variable 2); the long-term unemployment share, defined as the Percentage of unemployed persons who have been unemployed for 12 months or more (Variable 3).
Non-spherical clusters seem to be present in the data. The Gustafson and Kessel-like fuzzy k-means should be used for finding them.}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{Fclust}}, \code{\link{FKM}}, \code{\link{FKM.gk}}}
\examples{
## unemployment data
data(unemployment)
## fuzzy k-means (only spherical clusters)
unempFKM=FKM(unemployment,k=3)
## Gustafson and Kessel-like fuzzy k-means (non-spherical clusters)
unempFKM.gk=FKM.gk(unemployment,k=3,RS=10)}

\keyword{data}
\keyword{multivariate}
