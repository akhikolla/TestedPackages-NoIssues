#' chickn-package
#' 
#' The R package \code{chickn} implements the Chromatogram Hierarchical Compressive K-means with Nystrom approximation 
#' clustering approach. It is designed to cluster a large collection of high-resolution 
#' mass spectrometry signals (chromatographic profiles) relying on a compressed version of the data (a.k.a. data sketch). 
#' Data compression is achieved following the guidelines for Nystrom approximation provided by \insertCite{wang2019scalable}{chickn}
#' and the sketching operator from \insertCite{DBLP:journals/corr/KerivenBGP16}{chickn}.
#' The Filebacked Big Matrix (FBM) class from the [bigstatsr](https://github.com/privefl/bigstatsr) package 
#' is used to store and to manupulate matrices, which cannot be load in memory.  
#' @references 
#' \itemize{
#' \item Permiakova O, Guibert R, Kraut A, Fortin T, Hesse AM, Burger T (2020) "CHICKN: Extraction of peptide chromatographic 
#' elution profiles from large scale mass spectrometry data by means of Wasserstein compressive hierarchical cluster analysis." 
#' BMC Bioinformatics (under revision).
# \item \insertRef{wang2019scalable}{chickn}
# \item \insertRef{DBLP:journals/corr/KerivenBGP16}{chickn}.
#' }
#' @docType package
#' @author Olga Permiakova, Romain Guibert, Thomas Burger
#' @import bigstatsr RcppParallel mvnfast zipfR MASS pracma nloptr foreach doRNG parallel doParallel
#' @importFrom Rcpp evalCpp
#' @useDynLib chickn, .registration=TRUE
#' @name chickn
NULL  