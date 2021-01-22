\name{WangF}
\alias{WangF}
\title{To perform QTL mapping with wang method}
\description{
 Genome-wide Composite Interval Mapping
}
\usage{
WangF(pheRaw,genRaw,mapRaw1,yygg1,flagRIL,cov_en,Population,WalkSpeed,CriLOD)
}
\arguments{
	\item{pheRaw}{phenotype matrix.}
	\item{genRaw}{genotype matrix.}
	\item{mapRaw1}{linkage map matrix.}
	\item{yygg1}{the transformed covariate matrix .}
      \item{flagRIL}{if RIL or not.}
      \item{cov_en}{raw covariate matrix.}
	\item{Population}{population flag.}
	\item{WalkSpeed}{Walk speed for Genome-wide Scanning.(WalkSpeed=1).}
       \item{CriLOD}{Critical LOD scores for significant QTL (CriLOD=2.5).}
      	}
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, and Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
\dontrun{
data(DHdata)
readraw<-Readdata(file=DHdata,fileFormat="GCIM",fileICIMcov=NULL)
###
DoResult<-Dodata(fileFormat="GCIM",Population="DH",Model="Random",readraw)
###
ws<-WangF(pheRaw=DoResult$pheRaw,genRaw=DoResult$genRaw,
mapRaw1=DoResult$mapRaw1,yygg1=DoResult$yygg1,
flagRIL=DoResult$flagRIL,cov_en=DoResult$cov_en,
Population="DH",WalkSpeed=1,CriLOD=2.5)
}
}
