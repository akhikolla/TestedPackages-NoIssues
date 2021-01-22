#' Deconvolute Mixed Genomes with Unknown Proportions
#'
#' Traditional phasing programs are limited to diploid organisms.
#' Our method modifies Li and Stephens algorithm with Markov chain Monte Carlo
#' (MCMC) approaches, and builds a generic framework that allows haplotype
#' searches in a multiple infection setting. This package is primarily developed
#' as part of #' the Pf3k project, which is a global collaboration using the
#' latest sequencing technologies to provide a high-resolution view of natural
#' variation in the malaria parasite Plasmodium falciparum. Parasite DNA are
#' extracted from patient blood sample, which often contains more than one
#' parasite strain, with unknown proportions. This package is used for
#' deconvoluting mixed haplotypes, #' and reporting the mixture proportions from
#' each sample.
#'
#' @author
#' Zhu Sha
#'
#' Maintainer: Joe Zhu \email{sha.joe.zhu@gmail.com}
#'
#' @name DEploid-package
#' @docType package
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib _DEploid_dEploid
#' @useDynLib _DEploid_extractVcf
NULL
