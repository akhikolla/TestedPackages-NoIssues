\name{boscoclust}
\alias{boscoclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function to perform a co-clustering
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function runs a co-clustering algorithm on ordinal data by using the latent block model (see references for further  details). A BOS distribution is used, and the parameters inference is obtained using the SEM-Gbbs algorithm.
}
\usage{
boscoclust(x=matrix(0,nrow=1,ncol=1), idx_list=c(1), kr, kc, init, nbSEM, nbSEMburn, 
          nbRepeat=1, nbindmini, m=0, percentRandomB=0)
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
  Number of row classes.
  }
    \item{kc}{
  %%     ~~Describe \code{kc} here~~
  Vector of length D. The d^{th} element indicates the number of column clusters. %% Set to 0 to choose a classical multivariate BOS model.
  }
  \item{m}{
  %%     ~~Describe \code{m} here~~
  Vector of length D. The d^{th} element defines the number of levels of the ordinal data.
  }
  \item{nbSEM}{
  %%     ~~Describe \code{nbSEM} here~~
  Number of SEM-Gibbs iterations realized to estimate parameters.
  }
  \item{nbSEMburn}{
  %%     ~~Describe \code{nbSEMburn} here~~
  Number of SEM-Gibbs burn-in iterations for estimating parameters. This parameter must be inferior to nbSEM.
  }
  \item{nbRepeat}{
  %%     ~~Describe \code{nbSEMburn} here~~
  Number of times sampling on rows and columns will be done for each SEM-Gibbs iteration.
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
  Vector of length 2. Indicates the percentage of resampling when init is equal to "randomBurnin".
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
\item{@icl}{ICL value for co-clustering.}
\item{@name}{Name of the result.}
\item{@paramschain}{List of length nbSEMburn. The parameters of the blocks are stored for each iteration of the SEM-Gibbs algorithm.}
\item{@pichain}{List of length nbSEM. Item i is a vector of length kr that contains the row mixing proportions at iteration i.}
\item{@rhochain}{List of length nbSEM. Item i is a list of length D whose d^{th} element contains the column mixing proportions of the group of variables d, for iteration i.}
\item{@zc}{List of length D. The d^{th} item is a vector of length J[d] representing the column partitions for the group of variables d.}
\item{@zr}{Vector of length N with resulting row partitions.}
\item{@W}{List of length D. Item d is a matrix of dimension J*kc[d] such that W[j,h]=1 if j belongs to cluster h.}
\item{@m}{Vector of length D. The d^{th} element represents the number of levels of d^{th} group of variables.}
\item{@params}{List of length D. The d^{th} item represents the blocks parameters for a group of variables d.}
\item{@pi}{Vector of length kr. This corresponds to the row mixing proportions.}
\item{@rho}{List of length D. The d^{th} item represents the column mixing proportion for the d^{th} group of variables.}
\item{@xhat}{List of length D. The d^{th} item represents the dataset of the d^{th} group of variables, with missing values completed.}
\item{@zrchain}{Matrix of dimension nbSEM*N. Row i represents the row cluster partitions at iteration i.}  
\item{@zrchain}{List of length D. Item d is a matrix of dimension nbSEM*J[d]. Row i represents the column cluster partitions at iteration i.}

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
  \dontshow{
    data("Msimulated")

    nbSEM=15
    nbSEMburn=10
    nbindmini=1
    init = "kmeans"

    kr=2
    kc=2
    m=3

    res <- boscoclust(x=Msimulated,kr = kr, kc=kc, m=m, nbSEM=nbSEM,
                      nbindmini=nbindmini, nbSEMburn=nbSEMburn, init=init)
  }
  \donttest{
    
  library(ordinalClust)

  # loading the real dataset
  data("dataqol")
  set.seed(5)

  # loading the ordinal data
  M <- as.matrix(dataqol[,2:29])


  # defining different number of categories:
  m=4


  # defining number of row and column clusters
  krow = 4
  kcol = 4

  # configuration for the inference
  nbSEM=50
  nbSEMburn=40
  nbindmini=2
  init = "randomBurnin"
  percentRandomB=c(20,20)



  # Co-clustering execution
  object <- boscoclust(x=M,kr=krow,kc=kcol,m=m,nbSEM=nbSEM,
            nbSEMburn=nbSEMburn, nbindmini=nbindmini, init=init, percentRandomB=percentRandomB)

  }
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"ordinal-data"}% use one of  RShowDoc("KEYWORDS")
\keyword{"co-clustering"}% __ONLY ONE__ keyword per line
