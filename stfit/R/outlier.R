## Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
#' Image Outlier Detection
#'
#' @param mat data matrix. Each row is a row stacked image.
#'
#' @return a list containing the following entries:
#' \itemize{
#'   \item outidx: index of the outlier image
#'   \item outpct: percentage of outlier pixels corresponding to \code{outidx},
#'   \item outlst: a list of the same length as \code{outidx}, with each list the missing pixel index.
#' }
#' @export
#'
#' @examples
#' dfB = landsat106[landsat106$year >= 2000,]
#' matB = as.matrix(dfB[,-c(1:2)])
#' outlier(matB)
outlier <- function(mat){
  .outlier <- function(y){
    whisker = boxplot(y, plot = FALSE)$stats[c(1, 5)]
    which(y < whisker[1] | y > whisker[2])
  }

  lst = lapply(1:ncol(mat), function(i) .outlier(mat[,i])) ## length equals # of pixels
  ## a table of the outlier image index(row index of mat), higher value indicates higher
  ## probability that the corresponding image has problem.
  tbl = table(unlist(lst)) ## table of outlier numbers by image
  idx = as.numeric(names(tbl)) ## image index
  tot = apply(mat[idx, ], 1, function(x) sum(!is.na(x))) ## number of observed pixels
  outpct = tbl/tot
  if(length(idx) > 0){
    outlst = vector("list", length(idx))
    for(i in 1:length(idx)){
      for(j in 1:length(lst)){
        if(idx[i] %in% lst[[j]])
          outlst[[i]] = c(outlst[[i]], j)
      }
    }
  } else
    outlst = NULL
  return(list(outidx = idx, outpct = outpct, outlst = outlst))
}
