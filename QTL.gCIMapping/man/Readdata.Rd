\name{Readdata}
\alias{Readdata}
\title{Read raw data}
\description{
 Read raw data which have not been transformed
}
\usage{
Readdata(file,fileFormat,fileICIMcov)
}
\arguments{
	\item{file}{Dataset input}
	\item{fileFormat}{Format of dataset .}
	\item{fileICIMcov}{Format of covariate file of QTLIciMapping.}
      	}
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, and Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
\dontrun{
data(F2data)
Readdata(file=F2data,fileFormat="GCIM",fileICIMcov=NULL)
}
}
