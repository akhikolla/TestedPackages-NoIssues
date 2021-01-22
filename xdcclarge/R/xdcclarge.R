#' Package
#'
#' @description Functions for Estimating a (c)DCC-GARCH Model in large dimensions based on a publication by Engle et,al (2017) and Nakagawa et,al (2018).
#' This estimation method is consist of composite likelihood method by Pakel et al. (2014) and (Non-)linear shrinkage estimation of covariance matrices by Ledoit and Wolf (2004,2015,2016).
#'
#' @details To estimate the covariance matrix in financial time series,  
#'    it is necessary consider two important aspects: the cross section and the time series. 
#'    With regard to the cross section, we have the difficulty of correcting the biases of 
#'    the sample covariance matrix eigenvalues in a large number of time series.
#'    With regard to the time series aspect, we have to account for volatility clustering and time-varying correlations. 
#'    This package is implemented the improved estimation of the covariance matrix based on the following publications:
#'
#'    \itemize{\item Aielli, Gian Piero. (2013).
#'    Dynamic conditional correlation: on properties and estimation. Journal of Business & Economic Statistics 31: 282-99. <doi:10.1080/07350015.2013.771027>
#'    \item Engle, Robert F. (2002). 
#'    Dynamic conditional correlation: A simple class of multivariate generalized autoregressive conditional heteroskedasticity models. 
#'    Journal of Business & Economic Statistics 20: 339-50. <doi:10.1198/073500102288618487> 
#'    \item Engle, Robert F, Olivier Ledoit, and Michael Wolf. (2017).
#'    Large dynamic covariance matrices. Journal of Business & Economic Statistics, 1-13. <doi:10.1080/07350015.2017.1345683>
#'    \item Kei Nakagawa, Mitsuyoshi Imamura and Kenichi Yoshida. (2018). Risk-Based Portfolios with Large Dynamic Covariance Matrices.
#'    International Journal of Financial Studies, 6(2), 1-14. <doi:10.3390/ijfs6020052>
#'    \item Ledoit, O. and Wolf, M. (2004). A well-conditioned
#'    estimator for large-dimensional covariance matrices. Journal of
#'    Multivariate Analysis, 88(2). <doi:10.1016/S0047-259X(03)00096-4>
#'    \item Ledoit, O. and Wolf, M. (2012).
#'    Nonlinear shrinkage estimation of large-dimensional covariance matrices.
#'    Annals of Statistics, 40(2). <doi:10.1214/12-AOS989>
#'    \item Ledoit, O. and Wolf, M. (2015). Spectrum
#'    estimation: a unified framework for covariance matrix estimation and PCA in
#'    large dimensions. Journal of Multivariate Analysis, 139(2). <doi:10.1016/j.jmva.2015.04.006>
#'    \item Pakel, Cavit and Shephard, Neil and Sheppard, Kevin and Engle, Robert F. (2014). Fitting vast dimensional time-varying covariance models.
#'    Technical report <http://paneldataconference2015.ceu.hu/Program/Cavit-Pakel.pdf>
#'    }
#'
#' @name xdcclarge
#' @docType package
#' @import stats Rcpp nlshrink
NULL
