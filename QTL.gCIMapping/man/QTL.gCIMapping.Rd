\name{QTL.gCIMapping}
\alias{QTL.gCIMapping}
\alias{QTL.gCIMapping}
\docType{package}
\title{QTL Genome-Wide Composite Interval Mapping}
\description{
Conduct multiple quantitative trait loci (QTL) mapping under the framework of random-QTL-effect mixed linear model. First, each position on the genome is detected in order to construct a negative logarithm P-value curve against genome position. Then, all the peaks on each effect (additive or dominant) curve are viewed as potential QTL, all the effects of the potential QTL are included in a multi-QTL model, their effects are estimated by empirical Bayes in doubled haploid or by adaptive lasso in F2, and true QTL are identified by likelihood radio test.}
\details{
\tabular{ll}{
Package: \tab QTL.gCIMapping\cr
Type: \tab Package\cr
Version: \tab 3.3.1\cr
Date: \tab 2020-10-8\cr
Depends: \tab MASS,qtl\cr
Imports: \tab methods,openxlsx,stringr,Rcpp\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
}
\usage{
QTL.gCIMapping(file,fileFormat,fileICIMcov,Population,Model,WalkSpeed,
CriLOD,Likelihood,SetSeed,flagrqtl,DrawPlot,PlotFormat,Resolution,Trait,dir) 
 }
\arguments{
	\item{file}{File path and name in your computer,
       i.e.,"D:/Users/GCIM_Format_DH.csv".}
	\item{fileFormat}{Format for input file (GCIM, ICIM, Cart).}
	\item{fileICIMcov}{File path and name in your computer,
       i.e.,"D:/Users/ICIM_Cov.csv".}
	\item{Population}{BC1, BC2, DH, RIL, F2.}
	\item{Model}{Random (random model) or Fixed (fixed model) for QTL 
       effects.}
	\item{WalkSpeed}{Walk speed for Genome-wide Scanning.(WalkSpeed=1)}
       \item{CriLOD}{Critical LOD scores for significant QTL (CriLOD=2.5).}
       \item{Likelihood}{This parameter is only for F2 population, 
       including REML (restricted maximum likelihood ) and ML(maximum  
       likelihood).}
       \item{SetSeed}{In which the cross validation experiment is needed. Generally speaking, the random seed in the cross-validation experiment was set as 11001. If some known genes are not identified by the seed, users may try to use some new random seeds. At this case, one better result may be obtained.}
       \item{flagrqtl}{This parameter is only for F2 population,   flagrqtl="FALSE" in the first running. If the other software detects only one QTL in a neighborhood but this software finds two linked QTLs (one with additive effect and another with dominant effect) in the region, let flagrqtl="TRUE"}
      \item{DrawPlot}{This parameter is for all the populations, including FALSE and TRUE, DrawPlot=FALSE indicates no figure output, DrawPlot=TRUE indicates the output of the figure against genome position.}
	\item{PlotFormat}{This parameter is for all the figure files, including *.jpeg, *.png, *.tiff and *.pdf,i.e.,PlotFormat="tiff".}
	\item{Resolution}{This parameter is for all the figure files, including Low and High.}
       \item{Trait}{Trait=1:3 indicates the analysis from the first trait to the third trait.}
       \item{dir}{This parameter is for the save path,i.e.,"D:/Users".}

}
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, and Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\note{
     The User Manual file (name: Instruction.pdf) in the folder of 
     ".../QTL.gCIMapping/inst/doc".
}
\references{
1.  Wang Shi-Bo, Wen Yang-Jun, Ren Wen-Long, Ni Yuan-Li, Zhang Jin, Feng Jian-Ying, Zhang Yuan-Ming*. Mapping small-effect and linked quantitative trait loci for complex traits in backcross or DH populations via a multi-locus GWAS methodology. Scientific Reports 2016, 6: 29951.
2.  Wen Yang-Jun, Zhang Ya-Wen, Zhang Jin, Feng Jian-Ying, Jim M. Dunwell, Zhang Yuan-Ming*. An efficient multi-locus mixed model framework for the detection of small and linked QTLs in F2. Briefings in Bioinformatics 2019, 20(5): 1913-1924.
3.  Zhang Ya-Wen, Wen Yang-Jun, Jim M. Dunwell, Zhang Yuan-Ming*. QTL.gCIMapping.GUI v2.0: An R software for detecting small-effect and linked QTLs for quantitative traits in bi-parental segregation populations. Computational and Structural Biotechnology Journal 2020, 18: 59-65.
4.  Wen Yang-Jun, Zhang Ya-Wen, Zhang Jin, Feng Jian-Ying, Zhang Yuan-Ming*. The improved FASTmrEMMA and GCIM algorithms for genome-wide association and linkage studies in large mapping populations. The Crop Journal 2020, 8(5): 723-732.
}
\examples{
\dontrun{
data(F2data)
QTL.gCIMapping(file=F2data,Population="F2",WalkSpeed=1,CriLOD=2.5,
Trait=1,dir=tempdir()) 
}
}