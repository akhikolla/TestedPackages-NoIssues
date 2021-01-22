\name{bosclust}
\alias{bosclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function to perform a clustering
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function performs a clustering algorithm on ordinal data by using the multiple latent block model (see references for further  details). It allows the user to define D groups of variables that have different numbers of levels. The BOS distribution is used, and the parameters inference is obtained using the SEM-Gibbs algorithm.
}
\usage{
bosclust(x, idx_list=c(1), kr, init, nbSEM, nbSEMburn, 
        nbindmini, m=0, percentRandomB=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  %%     ~~Describe \code{x} here~~
  Matrix made of ordinal data of dimension N*Jtot. The features with the same numbers of levels must be placed side by side. The missing values should be coded as NA.
  }
  \item{idx_list}{
  %%     ~~Describe \code{x} here~~
  Vector of length D. This argument is useful when variables have different numbers of levels. Element d should indicate where the variables with number of levels m[d] begin in matrix x.
  }
  \item{kr}{
  %%     ~~Describe \code{kr} here~~
  Number of row clusters.
  }
  \item{m}{
  %%     ~~Describe \code{m} here~~
  Vector of length D. The d^{th} element defines the number of levels of the ordinal data.
  }
  \item{nbSEM}{
  %%     ~~Describe \code{nbSEM} here~~
  Number of SEM-Gibbs iterations realized to estimate the parameters.
  }
  \item{nbSEMburn}{
  %%     ~~Describe \code{nbSEMburn} here~~
  Number of SEM-Gibbs burn-in iterations for estimating parameters. This parameter must be inferior to nbSEM.
  }
  \item{nbindmini}{
  %%     ~~Describe \code{nbindmini} here~~
  Minimum number of cells belonging to a block.
  }
  \item{init}{
  %%     ~~Describe \code{init} here~~
  String that indicates the kind of initialisation. Must be one of the following words :  "kmeans", "random" or "randomBurnin".
  }
  \item{percentRandomB}{
  Vector of length 1. Indicates the percentage of resampling when init is equal to "randomBurnin".
  }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

  \item{@V}{Matrix of dimension N*kr such that V[i,g]=1 if i belongs to cluster g.}
  \item{@zr}{Vector of length N with resulting row partitions.}
  \item{@pi}{Vector of length kr. This corresponds to the row mixing proportions.}
  \item{@m}{Vector of length D. The d^{th} element represents the number of levels of the d^{th} group of variables.}
  \item{@icl}{ICL value for clustering.}       
  \item{@name}{Name of the result.}
  \item{@params}{List of length D. The d^{th} item stores the resulting position and precision parameters mu and pi.}
  \item{@paramschain}{List of length nbSEMburn. The parameters of the blocks are stored for each iteration of the SEM-Gibbs algorithm.}
  \item{@xhat}{List of length D. The d^{th} item represents the dataset of the d^{th} group of variables, with missing values completed.}
  \item{@zrchain}{Matrix of dimension nbSEM*N. Row i represents the row cluster partitions at iteration i.}    
  \item{@pichain}{List of length nbSEM. Item i is a vector of length kr that contains the row mixing proportions at iteration i.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
%%  ~~who you are~~
Margot Selosse, Julien Jacques, Christophe Biernacki.
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{

\donttest{
  library(ordinalClust)
  data("dataqol")
  set.seed(5)

  # loading the ordinal data
  M <- as.matrix(dataqol[,2:29])

  m = 4

  krow = 4

  nbSEM=50
  nbSEMburn=40
  nbindmini=2
  init = "random"


  object <- bosclust(x=M,kr=krow, m=m, nbSEM=nbSEM,
      nbSEMburn=nbSEMburn, nbindmini=nbindmini, init=init)
    }
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"ordinal-data"}% use one of  RShowDoc("KEYWORDS")
\keyword{"clustering"}% __ONLY ONE__ keyword per line
