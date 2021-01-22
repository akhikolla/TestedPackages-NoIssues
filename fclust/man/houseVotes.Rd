\name{houseVotes}
\alias{houseVotes}
\docType{data}
\title{
Congressional Voting Records Data}
\description{
1984 United Stated Congressional Voting Records for each of the U.S. House of
      Representatives Congressmen on the 16 key votes identified by the
      Congressional Quarterly Almanac.
}
\usage{data(houseVotes)}
\format{
  A data.frame with 435 rows on 17 columns (16 qualitative variables and 1 classification variable).
}
\details{
The data collect 1984 United Stated Congressional Voting Records for each of the 435 U.S. House of Representatives Congressmen on the 16 key votes identified by the Congressional Quarterly Almanac (CQA). The variable \code{class} splits the observations in \code{democrat} and \code{republican}. The qualitative variables refer to the votes on \code{handicapped-infants}, \code{water-project-cost-sharing}, \code{adoption-of-the-budget-resolution}, \code{physician-fee-freeze}, \code{el-salvador-aid}, \code{religious-groups-in-schools}, \code{anti-satellite-test-ban}, \code{aid-to-nicaraguan-contras}, \code{mx-missile}, \code{immigration}, \code{synfuels-corporation-cutback}, \code{education-spending}, \code{superfund-right-to-sue}, \code{crime}, \code{duty-free-exports}, and \code{export-administration-act-south-africa}. All these 16 variables are objects of class \code{factor} with three levels according to the CQA scheme: \code{y} refers to the types of votes ''voted for'', ''paired for'' and ''announced for''; \code{n} to ''voted against'', ''paired against'' and ''announced against''; code{yn} to ''voted present'', ''voted present to avoid conflict of interest'' and ''did not vote or otherwise make a position known''.}
\source{
https://archive.ics.uci.edu/ml/datasets/congressional+voting+records
}
\references{
Schlimmer, J.C., 1987. Concept acquisition through representational adjustment. Doctoral dissertation, Department of Information and Computer Science, University of California, Irvine, CA. 
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{NEFRC}}, \code{\link{NEFRC.noise}}}
\examples{
data(houseVotes)
X=houseVotes[,-1]
class=houseVotes[,1]

}
\keyword{data}
\keyword{multivariate}