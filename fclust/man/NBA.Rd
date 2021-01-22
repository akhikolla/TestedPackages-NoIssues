\name{NBA}
\alias{NBA}
\docType{data}
\title{
NBA teams data }
\description{
NBA team statistics from the 2017-2018 regular season.}
\usage{data(NBA)}
\format{
 A data.frame with 30 rows and 22 columns.
}
\details{
Data refer to some statistics of the NBA teams for the regular season 2017-2018. The teams are distinguished according to two classification variables.

The statistics are: number of wins (\code{W}), field goals made (\code{FGM}), field goals attempted (\code{FGA}), field goals percentage (\code{FGP}), 3 point field goals made (\code{3PM}), 3 point field goals attempted (\code{3PA}), 3 point field goals percentage (\code{3PP}), free throws made (\code{FTM}), free throws attempted (\code{FTA}), free throws percentage (\code{FTP}), offensive rebounds (\code{OREB}), defensive rebounds (\code{DREB}), assists (\code{AST}), turnovers (\code{TOV}), steals (\code{STL}), blocks (\code{BLK}), blocked field goal attempts (\code{BLKA}), personal fouls (\code{PF}), personal fouls drawn (\code{PFD}) and points (\code{PTS}). Moreover, reported are the conference (\code{Conference}) and the playoff appearance (\code{Playoff}).
}
\source{
https://stats.nba.com/teams/traditional/
}
\seealso{\code{\link{FKM}}}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\examples{
\dontrun{

data(NBA)
## A subset of variables is considered
X <- NBA[,c(4,7,10,11,12,13,14,15,16,17,20)]
clust.FKM=FKM(X=X,k=2:6,m=1.5,RS=50,stand=1,index="SIL.F",alpha=1)
summary(clust.FKM)
}
}
\keyword{data}
\keyword{multivariate}
