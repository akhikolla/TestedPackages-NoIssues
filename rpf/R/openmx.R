##' Convert an OpenMx MxModel object into an IFA group
##'
##' When \dQuote{minItemsPerScore} is passed, EAP scores will be computed
##' from the data and stored. Scores are required for some diagnostic
##' tests. See discussion of \dQuote{minItemsPerScore} in \link{EAPscores}.
##'
##' @template detail-group
##' @param mxModel MxModel object
##' @param data observed data (otherwise the data will be taken from the mxModel)
##' @param container an MxModel in which to search for the latent distribution matrices
##' @template arg-dots
##' @param minItemsPerScore minimum number of items required to compute a score (also see description)
##' @return a groups with item parameters and latent distribution
##' @seealso \href{https://cran.r-project.org/package=ifaTools}{ifaTools}
as.IFAgroup <- function(mxModel, data=NULL, container=NULL, ..., minItemsPerScore=NULL) {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}
	expectation <- mxModel$expectation
  if (!is(expectation, "MxExpectationBA81")) {
    stop(paste("Don't know how to create an IFA group from",
               class(expectation)))
  }
	if (missing(container)) {
		container <- mxModel
	}

  mat <- expectation$item
  if (length(grep("\\.", mat))) {
    stop(paste("Don't know how to obtain the item matrix", mat))
  }
  itemMat <- mxModel[[mat]]
  if (is.null(itemMat)) {
    stop(paste("Item matrix", mat, "not found"))
  }

  ret <- list(spec = expectation$ItemSpec,
              param = itemMat$values,
              free = itemMat$free,
	      labels = itemMat$labels,
	      # TODO maybe should include free variables in latent distribution?
	      uniqueFree = length(unique(itemMat$labels[itemMat$free], incomparables=NA)),
              qpoints = expectation$qpoints,
              qwidth = expectation$qwidth)

  mat <- expectation$mean
  if (length(grep("\\.", mat))) {
    meanMat <- eval(substitute(mxEval(theExpression, container),
			       list(theExpression = parse(text = mat)[[1]])))
    if (is.null(meanMat)) {
      stop(paste("Don't know how to obtain the mean matrix", mat))
    }
    ret$mean <- meanMat
  } else {
	  mxMat <- mxModel[[mat]]
	  if (!is.null(mxMat)) {
		  ret$mean <- mxMat$values
	  }
  }

  mat <- expectation$cov
  if (length(grep("\\.", mat))) {
    covMat <- eval(substitute(mxEval(theExpression, container),
				   list(theExpression = parse(text = mat)[[1]])))
    if (is.null(covMat)) {
      stop(paste("Don't know how to obtain the cov matrix", mat))
    }
    ret$cov <- covMat
  } else {
	  mxMat <- mxModel[[mat]]
	  if (!is.null(mxMat)) {
		  ret$cov <- mxMat$values
	  }
  }

  if (!missing(data)) {
    ret$data <- data
  } else if (!is.null(mxModel$data)) {
    mxData <- mxModel$data
    if (mxData$type != "raw") {
      stop(paste("Not sure how to handle data of type", mxData$type))
    }
    if (!is.null(mxData$frequency) && !is.na(mxData$frequency)) {
	    ret$freqColumn <- mxData$frequency
    }
    if (!is.null(mxData$weight) && !is.na(mxData$weight)) {
	    ret$weightColumn <- mxData$weight
    }
    ret$data <- mxData$observed
  }

  if (!is.na(expectation$weightColumn)) {
    ret$weightColumn <- expectation$weightColumn
    ret$observedStats <- nrow(ret$data) - 1
  } else {
	  if (ncol(ret$param) == ncol(ret$data)) {
		  freq <- tabulateRows(ret$data[orderCompletely(ret$data),])
		  ret$observedStats <- length(freq) - 1L
	  }
  }

  if (max(sapply(ret$spec, function(s) s$factors)) > 0 && !missing(minItemsPerScore)) {
    ret$minItemsPerScore <- minItemsPerScore
    ret$score <- EAPscores(ret, compressed=TRUE)
  }

  ret
}

#' Strip data and scores from an IFA group
#'
#' In addition, the freqColumn and weightColumn are reset to NULL.
#'
#' @template arg-grp
#' @template detail-group
#' @return
#' The same group without associated data.
#' @examples
#' spec <- list()
#' spec[1:3] <- list(rpf.grm(outcomes=3))
#' param <- sapply(spec, rpf.rparam)
#' data <- rpf.sample(5, spec, param)
#' colnames(param) <- colnames(data)
#' grp <- list(spec=spec, param=param, data=data, minItemsPerScore=1L)
#' grp$score <- EAPscores(grp)
#' str(grp)
#' grp <- stripData(grp)
#' str(grp)
stripData <- function(grp) {
	grp$data <- NULL
	grp$score <- NULL
	grp$weightColumn <- NULL
	grp$freqColumn <- NULL
	grp
}
