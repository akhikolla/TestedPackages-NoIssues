\name{Dodata}
\alias{Dodata}
\title{Process raw data}
\description{
 Process raw data for later use
}
\usage{
Dodata(fileFormat,Population,Model,readraw)
}
\arguments{
	\item{fileFormat}{Dataset format}
	\item{Population}{Population type.}
	\item{Model}{Random or fixed model.}
	\item{readraw}{Raw data.}
      	}
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, and Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
\dontrun{
data(DHdata)
readraw<-Readdata(file=DHdata,fileFormat="GCIM",fileICIMcov=NULL)
doda<-Dodata(fileFormat="GCIM",Population="DH",Model="Random",readraw)
}
}
