dist.matrix <- function (M, M2=NULL, method="cosine", p=2, normalized=FALSE, byrow=TRUE, convert=TRUE, as.dist=FALSE, terms=NULL, terms2=terms, skip.missing=FALSE) {
  method <- match.arg(method, c("cosine", "euclidean", "maximum", "manhattan", "minkowski", "canberra", "jaccard", "overlap"))
  similarity <- (method %in% c("cosine", "jaccard", "overlap")) && !convert
  symmetric <- !(method %in% c("overlap")) # FALSE if distance/similarity measure is asymmetric
  cross.distance <- !is.null(M2)  # TRUE if calculating (rectangular) cross-distance matrix
  need.nonneg <- method %in% c("jaccard", "overlap")
  
  if (method == "minkowski") {
    if (p == Inf) {
      method <- "maximum" # same as in rowNorms()
    } else if (p < 0 || !is.finite(p)) {
      stop("Minkowski p-distance can only be computed for 0 <= p < Inf")
    }
  }
  if (as.dist && similarity) stop("cannot create 'dist' object from similarity matrix")
  if (as.dist && cross.distance) stop("cannot create 'dist' object from cross-distance matrix")
  if (as.dist && !symmetric) stop("cannot create 'dist' object for asymmetric distance measure")

  ensure.nonneg <- function (M, name="M") {
    is.nonneg <- dsm.is.canonical(M)$nonneg
    if (is.na(is.nonneg)) is.nonneg <- dsm.is.canonical(M, nonneg.check=TRUE)$nonneg
    if (!is.nonneg) stop(sprintf("%s must be a non-negative matrix", name))
  }
  
  M <- find.canonical.matrix(M) # extract co-occurence matrix from DSM object, ensure canonical format
  sparse.M <- dsm.is.canonical(M)$sparse
  if (need.nonneg) ensure.nonneg(M)
  
  if (cross.distance) {
    M2 <- find.canonical.matrix(M2)
    sparse.M2 <- dsm.is.canonical(M2)$sparse
    if (byrow) {
      if (ncol(M) != ncol(M2)) stop("M and M2 are not conformable (must have same number of columns)")
    } else {
      if (nrow(M) != nrow(M2)) stop("M and M2 are not conformable (must have same number of rows)")        
    }
    if (need.nonneg) ensure.nonneg(M2, name="M2")
  }

  if (!is.null(terms) || !is.null(terms2)) {
    targets.M <- if (byrow) rownames(M) else colnames(M)
    targets.M2 <- if (is.null(M2)) targets.M else if (byrow) rownames(M2) else colnames(M2)

    if (!missing(terms2)) cross.distance <- TRUE # if different filters are applied, we're always dealing with a cross-distance calculation

    ## if cross.distance is FALSE, M2 must be NULL and terms2 unspecified (i.e. M2=M and terms2=terms), so leave M2 set to NULL
    if (cross.distance) {
      if (!is.null(terms2)) {
        terms2 <- as.character(terms2) # in case terms2 is a factor
        found <- terms2 %in% targets.M2
        if (!all(found) && !skip.missing) stop("second term(s) not found in M2: ", paste(terms2[!found], collapse=", "))
        terms2 <- terms2[found]
        if (is.null(M2)) {
          M2 <- if (byrow) M[terms2, , drop=FALSE] else M[ , terms2, drop=FALSE] # need to process terms2 first before overwriting M below
        } else {
          M2 <- if (byrow) M2[terms2, , drop=FALSE] else M2[ , terms2, drop=FALSE]
        }
      } else {
        if (is.null(M2)) M2 <- M # cross-distances with terms2=NULL -> M2 = copy of M before subsetting
      }
      sparse.M2 <- dsm.is.canonical(M2)$sparse # define/update sparse.M2
    }

    if (!is.null(terms)) {
      terms <- as.character(terms) # in case terms is a factor
      found <- terms %in% targets.M
      if (!all(found) && !skip.missing) stop("first term(s) not found in M: ", paste(terms[!found], collapse=", "))
      terms <- terms[found]
      M <- if (byrow) M[terms, , drop=FALSE] else M[ , terms, drop=FALSE]
    }
  }
  
  if (method == "cosine") {
    ## cosine / angular measure is computed as very efficient matrix crossproduct
    
    if (byrow) {
      result <- if (is.null(M2)) tcrossprod(M) else tcrossprod(M, M2)
    } else {
      result <- if (is.null(M2)) crossprod(M) else crossprod(M, M2)
    }
    result <- as.matrix(result) # ensure that cosine similarity matrix is in dense representation
    if (!normalized) {
      norms.M <- if (byrow) rowNorms(M, "euclidean") else colNorms(M, "euclidean")
      if (is.null(M2)) {
        norms.M2 <- norms.M
      } else {
        norms.M2 <- if (byrow) rowNorms(M2, "euclidean") else colNorms(M2, "euclidean")        
      }
      result <- scaleMargins(result, rows=1/norms.M, cols=1/norms.M2, duplicate=FALSE) # transform in-place (newly allocated above)
    }

    if (convert) {
      transform_code <- 0L # cosine -> angle transformation
      tol <- 1e-12
      result <- CPP_similarity_to_distance(result, transform_code, tol, duplicate=FALSE) # can operate inplace on <result>
    }    
    rownames(result) <- if (byrow) rownames(M) else colnames(M)
    colnames(result) <- if (is.null(M2)) rownames(result) else if (byrow) rownames(M2) else colnames(M2)

  } else {
    ## other distance measures are implemented in C code, working on columns (transposed matrix) for efficiency

    .M <- if (byrow) t(M) else M
    if (cross.distance) {
      if (sparse.M != sparse.M2) stop("M and M2 must either be both in dense format or both in sparse format")
      .M2 <- if (byrow) t(M2) else M2
    } else {
      .M2 <- .M
    }

    method.code <- switch(method, euclidean=0, maximum=1, manhattan=2, minkowski=3, canberra=4, jaccard=5, overlap=6) # must be kept in sync with C code
    param1 <- switch(method, euclidean=0, maximum=0, manhattan=0, minkowski=p, canberra=0, jaccard=0, overlap=0)
  
    if (sparse.M) {
      result <- CPP_col_dist_sparse(ncol(.M), .M@p, .M@i, .M@x, ncol(.M2), .M2@p, .M2@i, .M2@x, method.code, param1, symmetric && !cross.distance)
    } else {
      result <- CPP_col_dist_dense(.M, .M2, method.code, param1, symmetric && !cross.distance)
    }
    if (method %in% c("jaccard", "overlap")) {
      if (method == "overlap" && !normalized) {
        ## asymmetric overlap relative to x, so values must be normalised with ||x||_1 = sum x_i
        norms.M <- if (byrow) rowNorms(M, "manhattan") else colNorms(M, "manhattan")
        result <- scaleMargins(result, rows=1/norms.M, duplicate=FALSE) # can operate inplace on <result>
        idx <- norms.M == 0   # special case: o(0, x) = 1
        if (any(idx)) result[idx, ] <- 1 
      }
      if (convert) {
        transform_code <- 1L # d = 1 - sim is a metric (jaccard) or dissimilarity (overlap)
        result <- CPP_similarity_to_distance(result, transform_code, 0, duplicate=FALSE) # can operate inplace on <result>
      }
    }
    dimnames(result) <- list(colnames(.M), colnames(.M2))
  }

  if (as.dist) {
    as.dist(result)
  } else {
    class(result) <- c("dist.matrix", "matrix")
    if (similarity) attr(result, "similarity") <- TRUE
    if (symmetric && !cross.distance) attr(result, "symmetric") <- TRUE
    result
  }
}
