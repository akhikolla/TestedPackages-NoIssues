\name{WenS}
\alias{WenS}
\title{The second step of Wen method}
\description{
An efficient multi-locus mixed model framework for the detection of small and linked QTLs in F2}
\usage{
WenS(flag,CriLOD,NUM,pheRaw,Likelihood,setseed,flagrqtl,yygg,mx,phe,
chr_name,v.map,gen.raw,a.gen.orig,d.gen.orig,n,names.insert2,X.ad.tran.data,X.ad.t4,dir)
}
\arguments{
      \item{flag}{random or fix model.}
	\item{CriLOD}{LOD score.}
	\item{NUM}{the number of trait.}
	\item{pheRaw}{raw phenotype matrix .}
	\item{Likelihood}{likelihood function.}
      \item{setseed}{random seed set in which,
            the cross validation is needed.}
	\item{flagrqtl}{do CIM or not.}
      \item{yygg}{covariate matrix.}
      \item{mx}{raw genotype matrix.}
      \item{phe}{phenotype matrix.}
	\item{chr_name}{chromosome name.}
	\item{v.map}{linkage map matrix.}
	\item{gen.raw}{raw genotype matrix.}
	\item{a.gen.orig}{additive genotype matrix.}
	\item{d.gen.orig}{dominant genotype matrix.}
      \item{n}{number of individual.}
      \item{names.insert2}{linkage map after insert.}
      \item{X.ad.tran.data}{genotype matrix after insert.}
	\item{X.ad.t4}{genotype matrix.}
      \item{dir}{file storage path.}
       }
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
\dontrun{
data(genf2)
data(phef2)
data(mapf2)
WEN1re<-WenF(pheRaw=phef2,genRaw=genf2,mapRaw1=mapf2,
yygg1=NULL,cov_en=NULL,WalkSpeed=1,CriLOD=2.5,dir=tempdir())
###
ws<-WenS(flag=1,CriLOD=2.5,NUM=1,pheRaw=phef2,
Likelihood="REML",setseed=11001,flagrqtl=FALSE,
yygg=WEN1re$yygg,mx=WEN1re$mx,phe=WEN1re$phe,
chr_name=WEN1re$chr_name,v.map=WEN1re$v.map,
gen.raw=WEN1re$gen.raw,a.gen.orig=WEN1re$a.gen.orig,
d.gen.orig=WEN1re$d.gen.orig,n=WEN1re$n,
names.insert2=WEN1re$names.insert2,
X.ad.tran.data=WEN1re$X.ad.tran.data,
X.ad.t4=WEN1re$X.ad.t4,dir=tempdir())
}
}

