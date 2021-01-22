# Purpose: Package documentation
# Updated: 2020/10/03

#' @useDynLib RNOmni, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' RNOmni: Rank-Normal Omnibus Association Testing
#'
#' Implementation of genetic association tests that incorporate the rank-based
#' inverse normal transformation (INT) \code{\link{RankNorm}}. In direct-INT
#' \code{\link{DINT}}, the phenotype itself is transformed. In indirect-INT
#' \code{\link{IINT}}, phenotypic residuals are transformed. The omnibus INT
#' \code{\link{OINT}} test adaptively combines the D-INT and I-INT tests into a
#' single robust and statistically powerful procedure. See McCaw ZR, Lane JM,
#' Saxena R, Redline S, Lin X. "Operating characteristics of the rank-based
#' inverse normal transformation for quantitative trait analysis in genome-wide
#' association studies." Biometrics, 2019. <https://doi.org/10.1111/biom.13214>.
#' 
#' @author Zachary R. McCaw
#' @docType package
#' @name RNOmni-help
NULL
