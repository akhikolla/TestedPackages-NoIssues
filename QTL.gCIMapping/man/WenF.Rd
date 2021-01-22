\name{WenF}
\alias{WenF}
\title{To perform QTL mapping with Wen method}
\description{
An efficient multi-locus mixed model framework for the detection of small and linked QTLs in F2}
\usage{
WenF(pheRaw,genRaw,mapRaw1,yygg1,cov_en,WalkSpeed,CriLOD,dir)
}
\arguments{
      \item{pheRaw}{phenotype matrix.}
	\item{genRaw}{genotype matrix.}
	\item{mapRaw1}{linkage map matrix.}
	\item{yygg1}{the transformed covariate matrix .}
	\item{cov_en}{raw covariate matrix.}
	\item{WalkSpeed}{Walk speed for Genome-wide Scanning.(WalkSpeed=1).}
       \item{CriLOD}{Critical LOD scores for significant QTL (CriLOD=2.5).}
       \item{dir}{file path in your computer.}
       }
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, and Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
\dontrun{
data(F2data)
###
readraw<-Readdata(file=F2data,fileFormat="GCIM",fileICIMcov=NULL)
###
DoResult<-Dodata(fileFormat="GCIM",Population="F2",Model="Random",readraw)
###
wf<-WenF(pheRaw=DoResult$pheRaw,genRaw=DoResult$genRaw,
mapRaw1=DoResult$mapRaw1,yygg1=DoResult$yygg1,
cov_en=DoResult$cov_en,WalkSpeed=1,CriLOD=2.5,
dir=tempdir())
}
}



