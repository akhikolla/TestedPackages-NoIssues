\name{estimaterates_f}
\alias{estimaterates_f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate substitution rate matrix.
}
\description{
Beta version of code; still in development. More or less identical to \code{\link{estimaterates}}. Here, some entries of the substitution rate matrix can be specified by the user to be held fixed during the parameter estimation.}
\usage{
estimaterates_f(usertree = NULL, userphyl = NULL, matchtipstodata = FALSE,
              unobserved = NULL, alphabet = NULL, modelmat = NULL, 
              bgtype = "listofnodes", bg = NULL, partition = NULL, 
              ratevar = FALSE, nocat = 4, reversible = FALSE, numhessian = TRUE, 
              rootprob = NULL, rpvec = NULL, init = 0.9, lowli = 0.001, upli = 100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{usertree}{
Rooted binary tree of class "phylo". Read in Newick tree using read.tree() from package ape before passing this argument. The branch lengths must be in expected substitutions per site (but see lowli and upli arguments). Trees estimated using MrBayes and BEAST, for instance, yield branch lengths using that scale. If an unrooted tree is available, look into the ``root" and ``midpoint.root" functions from the APE and phytools packages, respectively.
}
  \item{userphyl}{
A matrix or data frame of phyletic patterns. Rows represent the discrete character patterns and columns represent the taxa. These data can be numeric or character (but not both). 
}
  \item{matchtipstodata}{
The default is FALSE, which means that the user must ensure that the ordering of the taxa in the data matrix must match the internal ordering of tip labels of the tree. Set to TRUE, if the column names of the data matrix, i.e., the taxa names, are all present, with the same spelling and notation in the tip labels of the tree provided, in which case, the restriction on ordering is not necessary.
}
  \item{unobserved}{
A matrix of unobserved phyletic patterns, representing possible sampling or acquisition bias. Each row should be a unique phyletic pattern.
}
  \item{alphabet}{
The set of discrete characters used. May be integer only or character only.
}
  \item{modelmat}{
Accepts the same arguments as \code{\link{estimaterates}}.  Here, some entries of the substitution rate matrix can be specified by the user to be held fixed during the parameter estimation. A square matrix can be input that consists of integer values denoting the indices of the substitution rates to be estimated and entries with decimal values (that are not equivalent to an integer in R, i.e., something other than, for example, 1.0, 2.0, etc.) to be held fixed. See example.
}
  \item{bgtype}{
Use this to group branches hypothesized to follow the same rates (but differ from other branches). If clade-specific insertion and deletion rates are required to be estimated, use argument option "ancestornodes". If, on the other hand, a group of branches (not in a clade) are hypothesized to follow the same rates, use argument option "listofnodes".}
  \item{bg}{
A vector of nodes should be provided if the "ancestornodes" option was chosen for argument "bgtype". If, on the other hand, "listofnodes" was chosen, a list should be provided with each element of the list being a vector of nodes that limit the branches that follow the same rates. See examples and vignette.
}
  \item{partition}{
A list of vectors (of sites) subject to different evolutionary constraints. For example, supplying list(c(1:2500), c(2501:5000)) means that sites 1 through 2500 follow their own substitution rates distinct from sites 2501 through 5000. These sites correspond to the rows of the data supplied in userphyl. Partition models can be fitted with a common (or unique) gamma distribution for rate variation over all partitions, and/or common root probabilities can be estimated over all partitions.
}
  \item{ratevar}{
Default option is FALSE. Specifying "discgamma" implements the discrete gamma approximation model of Yang, 1994. Even if a partition of sites is specified, this option uses the same \eqn{\alpha} parameter for the gamma distribution over all partitions. Specifying "partitionspecificgamma" implements the discrete gamma approximation model in each partition separately.
The number of categories can be specified using the "nocat" argument. See examples and vignette.
}
  \item{nocat}{
The number of categories for the discrete gamma approximation. 
}
  \item{reversible}{
This option forces a model to be reversible, i.e., the flow from state 'a' to 'b' is the same as the flow from state 'b' to 'a'. Only symmetric transition matrices can be specified in modelmat with this option. Inspired by the DiscML package.
}
  \item{numhessian}{
Set to FALSE if standard errors are not required. This speeds up the algorithm. Default is TRUE. Although the function being used to calculate these errors is reliable, it is rare but possible that the errors are not calculated (due to approximation-associated issues while calculating the Hessian). In this case, try bootstrapping with numhessian=FALSE.
}
  \item{rootprob}{
Four options are available: "equal", "stationary", "maxlik", and "user". 
 \itemize{
 \item Option "equal" means that all the discrete characters are given equal weight at the root. 
 \item Using "stationary" means that the discrete characters are weighted at the root by the stationary frequencies implied by the substitution rate matrix. Note that these can differ based on the partition. In the case that the "bgtype" argument is also provided, an average of the stationary frequencies of all branch groupings is used.
 \item If "maxlik" is supplied, then the root frequencies are also estimated. These do not differ based on the partition.
 \item If "user" is supplied here, a vector of root frequency parameters can be provided to argument "rpvec".
 }
}
  \item{rpvec}{
If option "user" is specified for argument "rootprob", supply a vector of the same length as that provided to argument "alphabet", representing the root frequency parameters.
}
\item{init     }{	
Initial value for the rates. The default value is 0.9.}
\item{lowli}{
  For finer control of the boundaries of the optimization problem. The default value is 0.001. This usually suffices if the branch lengths are in expected substitutions per site. However, if branch lengths are in different units, this should be changed accordingly. 
}
  \item{upli}{
  For finer control of the boundaries of the optimization problem. The default value is 100. This usually suffices if the branch lengths are in expected substitutions per site. However, if branch lengths are in different units, this should be changed accordingly. 
}
  \item{\dots}{
Passing other arguments to the optimization algorithm "nlminb". For example, control = list(trace = 5) will print progress at every 5th iteration.}
}
\details{
See vignette for detailed examples.
}
\value{
All arguments used while calling the estimaterates function are attached in a list. The following components are also returned in the same list:
\item{call        }{Function call used.}
\item{conv        }{A vector of convergence indicators for the model run. 0 denotes successful convergence.}
\item{time        }{Time taken in seconds.}
\item{tree        }{The phylogenetic tree used.}
\item{bg          }{List of group of nodes that capture branches that follow unique substitution rates.}
\item{results     }{List of results including the parameter estimates for each unique entry of modelmat (excluding zeros), standard errors, number of parameters fit, and AIC and BIC values. Furthermore, details from the optimization routine applied are also available. Use ...$results$wop$parsep.}
\item{data_red    }{Unique phyletic patterns observed.}
\item{w           }{List of number of times each gene phyletic pattern was observed.}
}

\author{
Utkarsh J. Dang and G. Brian Golding

\email{udang@binghamton.edu}
}
\note{
Experimental function
}
\references{
Paradis, E., Claude, J. & Strimmer, K. 2004. APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290. R package version 3.2.

Tane Kim and Weilong Hao 2014. DiscML: An R package for estimating evolutionary rates of discrete characters using maximum likelihood. R package version 1.0.1. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
##############
\dontrun{
rm(list=ls())
set.seed(1) #set seed
library(markophylo)
tree1 <- ape::rtree(4, br = runif(6, 0.01, 0.3)) #generate a random 4-taxa tree (6 branches) 
bf <- c(0.33, 0.33, 0.34) #vector of equal root frequencies
gf <- 5000 #number of character patterns to be simulated
rootseq <- sample(1:3, gf, replace = TRUE, prob = bf) #simulate the MRCA sequence (1s and 2s).
par1 <- 1.5 
par2 <- 3   
par3 <- 5
Q_sim <- matrix(c(NA, 0.6, par2, par1, NA, 0.8, par2, par3, NA), 3, 3) 
#generated substitution rate matrix
diag(Q_sim) <- -rowSums(Q_sim, na.rm = TRUE)
print(Q_sim)
#simulate data using the geiger package and Q_sim
data_sim <- sapply(1:gf, FUN = function(X) geiger::sim.char(tree1, 
Q_sim, model = "discrete", root = rootseq[X])) 
simulateddat <- list(data = t(data_sim), Q = Q_sim, tree = tree1)

model1 <- markophylo::estimaterates_f(usertree = simulateddat$tree, 
userphyl = simulateddat$data, alphabet = c(1, 2, 3), rootprob = "equal", 
modelmat = matrix(c(NA, 0.6, 2, 1, NA, 0.8, 2, 3, NA), 3, 3))
}
}