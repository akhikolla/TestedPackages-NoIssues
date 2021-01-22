\name{markerinsert}
\alias{markerinsert}
\title{To insert marker in genotype.}
\description{
a method that can insert marker in genotype.}
\usage{
markerinsert(mp,geno,map,cl,gg1,gg2,gg0,flagRIL)
}
\arguments{
      \item{mp}{linkage map matrix after insert.}
	\item{geno}{genotype matrix.}
	\item{map}{linkage map matrix.}
	\item{cl}{walk speed.}
	\item{gg1}{raw covariate matrix.}
	\item{gg2}{code for type 1.}
       \item{gg0}{code for missing.}
       \item{flagRIL}{RIL population or not.}
       }
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
\dontrun{
mp<-matrix(c(197.9196,198.7536,199.5876,200.4216,201.2453,
202.0691,202.8928,203.7521,204.6113,205.4706,206.3298,207.1891,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,3,3,3,3,3,3,
1,1,1,2,2,2,3,3,3,3,3,3,1,2,3,4,5,6,7,8,9,10,11,12),12,5)
map<-matrix(c(1,1,1,1,197.9196,200.4216,202.8928,207.1891),4,2)
geno<-matrix(c(1,99,99,99),1,4)
mark_insert<-QTL.gCIMapping::markerinsert(mp,geno,map,cl=1,gg1=1,gg2=-1,
gg0=99,flagRIL=1)
}
}


