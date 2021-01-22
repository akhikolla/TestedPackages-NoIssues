## internal helper functions for [rc]bind.dsm() and merge.dsm()

# merge tables of marginal frequencies (either columns or rows, as indicated by margin argument)
#  - all tables must have variables "term" and (optionally) "f"; metdata variables taken from first table
#  - marginal frequencies are checked for equality; if there are differences, pmax() is assigned to result table with attribute "adjusted" set
#  - mode="same" requires that all tables have exactly the same terms in the same order and aborts with error message otherwise
#  - mode="intersection" preserves only terms that occur in all tables and reorders them if necessary
#  - mode="union" collects terms from all tables (not implemented yet)
.combine.marginals <- function (tbls, margin=c("column", "row"), mode=c("same", "intersection", "union")) {
  mode <- match.arg(mode)
  if (mode == "union") stop("'union' operation has not been implemented yet")
  margin <- match.arg(margin)
  stopifnot(length(tbls) >= 1)

  have.term.vec <- sapply(tbls, function (.t) "term" %in% colnames(.t))
  if (!all(have.term.vec)) stop("all ", margin, " marginal tables must contain variable 'term'")
  have.f.vec <- sapply(tbls, function (.t) "f" %in% colnames(.t))
  have.f <- all(have.f.vec)
  if (any(have.f.vec) && !have.f) stop("either all ", margin, " marginal tables must contain variable 'f', or none of them may")

  first.tbl <- tbls[[1]]
  other.tbls <- tbls[-1]
  terms <- first.tbl$term
  if (mode == "same") {
    length.ok <- all( sapply(other.tbls, function (.t) nrow(.t) == length(terms)) )
    if (!length.ok) stop("all DSM objects must have the same number of ", margin, "s")
    terms.ok <- all( sapply(other.tbls, function (.t) all(.t$term == terms)) )
    if (!terms.ok) stop("all DSM objects must have the same ", margin, " labels (terms or features)")
  } else {
    for (.t in other.tbls) terms <- intersect(terms, .t$term)
    if (length(terms) < 1) stop("DSM objects have no common ", margin, " labels, cannot merge")
    # intersect() keeps original ordering from first.tbl, so don't re-sort terms
  }

  adjusted <- FALSE
  terms.idx <- match(terms, first.tbl$term)

  res <- data.frame(term=terms, stringsAsFactors=FALSE)
  if (have.f) {
    F <- first.tbl$f[terms.idx] # init marginal frequencies from first table
    for (.t in other.tbls) {
      .idx <- match(terms, .t$term)
      .f <- .t$f[.idx]
      if (any(.f != F)) {
        F <- pmax(F, .f)
        adjusted <- TRUE
      }
    }
    res$f <- F
  }

  other.vars <- setdiff(colnames(first.tbl), c("term", "f")) # metadata variables from first table
  if (length(other.vars) > 0) res <- cbind(res, first.tbl[terms.idx, other.vars, drop=FALSE])

  attr(res, "adjusted") <- adjusted
  return(res)
}

# rbind tables of marginal frequencies (i.e. joining terms from all tables) with some consistency checks
#  - all tables must have exactly the same variables, and terms must be unique across all tables
#  - margin argument is currently used for error messages only
#  - optional term.suffix contains strings to be appended to terms from each table in order to make them unique
.bind.marginals <- function (tbls, margin=c("row", "column"), term.suffix=NULL) {
  margin <- match.arg(margin)
  stopifnot(length(tbls) >= 1)

  have.term.vec <- sapply(tbls, function (.t) "term" %in% colnames(.t))
  if (!all(have.term.vec)) stop("all ", margin, " marginal tables must contain variable 'term'")

  first.tbl <- tbls[[1]]
  other.tbls <- tbls[-1]

  vars <- colnames(first.tbl) # check that all tables are compatible, i.e. have the same variables in the same order
  vars.same <- all( sapply(other.tbls, function (.t) all( colnames(.t) == vars )) )
  if (!vars.same) stop(margin, " information tables of all DSM objects must be compatible (same variables)")

  n.items <- sapply(tbls, function (.t) nrow(.t)) # number of entries from each DSM
  res <- do.call(rbind, tbls)

  if (!is.null(term.suffix)) {
    stopifnot(length(term.suffix) == length(tbls))
    res$orig.term <- res$term
    res$orig.part <- rep(term.suffix, n.items)
    res$term <- paste(res$orig.term, res$orig.part, sep="")
  }
  if (any(duplicated(res$term))) stop(margin, " labels must be unique across all DSM objects")

  rownames(res) <- 1:nrow(res)
  return(res)
}
