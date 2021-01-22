check_ends_at_least_starts <- function(x, stop_if_false=FALSE) {
	result <- TRUE
	if (!is.matrix(x))
		result <- FALSE
	else if (!(ncol(x)==2L))
		result <- FALSE
	else if (!(is.logical(as.vector(x))|is.numeric(as.vector(x))))
		result <- FALSE
	else if (any(x[,1L] > x[,2L]))
		result <- FALSE
	if (stop_if_false && !result) {
		stop("Intervals must be specified by two-column numeric matrices with all end points (col 2) greater than or equal to start points (col 1)")
	} 
	result
}

#' @title Stich together touching intervals and remove empty intervals
#' @description Given an integer matrix specifying disjoint intervals sorted by start position, merge intervals with matching start and ends, and remove intervals of length zero.
#' @param x Integer matrix of two columns, the first column giving the (inclusive) start points of intervals and the second column giving the corresponding (exclusive) end points.
#' @template return-intervals
#' @export
#' @examples 
#' stitch(cbind(1:2, 2:3))
stitch <- function(x) {
	check_ends_at_least_starts(x, TRUE)
	xne <- x[x[,2L] > x[,1L],,drop=FALSE]
	if (nrow(xne) < 2L)
		return(xne)
	touching <- (xne[-1L,1L]-xne[-nrow(xne),2L]) == 0L
	cbind(c(xne[1L,1L], xne[-1L,1L][!touching]), c(xne[-nrow(xne),2L][!touching], xne[nrow(xne),2L]))
}

#' @title Check intervals are detached, sorted and non-empty. 
#' @description Check that \code{x} is an integer matrix specifying intervals, that the specified intervals are detached (i.e. non-overlapping/disjoint and non-touching) and that it is sorted (given that the intervals are detached, sorting by start position gives a unique result), and that the start points are greater than the end points (i.e. that they are non-empty/the lengths of all intervals is greater than zero).
#' @param x Integer matrix of two columns, the first column giving the (inclusive) start points of intervals and the second column giving the corresponding (exclusive) end points.
#' @return Boolean value.
#' @export
#' @examples
#' detached_sorted_nonempty(cbind(1:2, 2:3)) 
#' detached_sorted_nonempty(cbind(c(1, 3), c(2, 4))) 
#' detached_sorted_nonempty(cbind(1, 1)) 
detached_sorted_nonempty <- function(x) {
	if (!check_ends_at_least_starts(x, FALSE))
		return(FALSE)
	if (!all(diff(as.integer(t(x)))>0L))
		return(FALSE)
	TRUE
}

#' @title Compute overlaps of two sets of detached and sorted intervals
#' @rdname nonoverlapping
#' @description Find intervals satisfying particular conditions, including corresponding base R functions \code{intersect} (i.e. find intersections of intervals), \code{union} (i.e. unions of intervals) and \code{setdiff} (i.e. finding intervals which are contained in one set of intervals but not another).
#' @param x Integer matrix of two columns, the first column giving the (inclusive) start points of intervals and the second column giving the corresponding (exclusive) end points.
#' @param y Same as \code{x}.
#' @param check Boolean value determining whether to check that the intervals specified in arguments \code{x} and \code{y} are sorted and non-overlapping (uses function \code{\link{detached_sorted_nonempty}}). Defaults to \code{TRUE}, but setting to \code{FALSE} may allow faster execution.
#' @param in_x Boolean value determining whether to flag \code{TRUE} on intervals contained in \code{x}.
#' @param in_y Boolean value determining whether to flag \code{TRUE} on intervals contained in \code{y}.
#' @param op Character value specifying operator used to combine flags for each interval, either \code{"and"} or \code{"or"}.
#' @param ... Additional arguments to be passed to \code{overlaps}.
#' @template return-intervals
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib IntervalSurgeon
#' @examples
#' intersects(cbind(1, 3), cbind(2, 4))
#' setdiffs(cbind(1, 3), cbind(2, 4))
#' unions(cbind(1, 3), cbind(2, 4))
overlaps <- function(x, y, check=TRUE, in_x=TRUE, in_y=TRUE, op="and") {
	if (check) {
		stopifnot(detached_sorted_nonempty(x))
		stopifnot(detached_sorted_nonempty(y))
	}
	pts <- sort(unique(c(x, y)))
	if (length(pts) == 0L)
		return(matrix(nrow=0L, ncol=2L, data=integer(0L)))
	tick <- dash_set_overlaps(starts1=x[,1L],ends1=x[,2L],starts2=y[,1L],ends2=y[,2L],state1=in_x,state2=in_y,op_is_and=op=="and",pts=pts)
	stitch(cbind(pts[-length(pts)],pts[-1L])[tick,,drop=FALSE])
}

#' @rdname nonoverlapping
#' @export
intersects <- function(x, y, ...) {
	overlaps(x=x, y=y, in_x=TRUE, in_y=TRUE, op="and", ...)
}

#' @rdname nonoverlapping
#' @export
unions <- function(x, y, ...) overlaps(x=x, y=y, in_x=TRUE, in_y=TRUE, op="or", ...)

#' @rdname nonoverlapping
#' @export
setdiffs <- function(x, y, ...) overlaps(x=x, y=y, in_x=TRUE, in_y=FALSE, op="and", ...)

#' @title Get break points for set of intervals
#' @description Get the sorted set start points and end points for a set of intervals specified as an integer matrix.
#' @param x Integer matrix of two columns, the first column giving the (inclusive) start points of intervals and the second column giving the corresponding (exclusive) end points.
#' @return Ordered integer vector of unique interval start/end points.
#' @export
#' @examples
#' breaks(cbind(2*1:5, 3*1:5))
breaks <- function(x) {
	sort(unique(as.integer(x)))
}

#' @title Get the sections from a set of interval breaks
#' @description Given a set of interval breaks (see \code{\link{breaks}}), generate a new set of intervals, the `sections', which partitions the full range of the given set, with an interval between every `break' (i.e. start/end point) in the given set.
#' @param x Sorted integer vector.
#' @template return-intervals
#' @export
#' @examples
#' sections(1:10)
sections <- function(x) {
	cbind(x[-length(x)], x[-1L])
}

#' @title Depth of piled intervals
#' @description Get the depth of piled intervals for each section in the sections of \code{x} (see \code{\link{sections}}).
#' @param x Integer matrix of two columns, the first column giving the (inclusive) start points of intervals and the second column giving the corresponding (exclusive) end points.
#' @param include_intervals Logical value determining whether the function should return a vector of depths at each `section' in the range of \code{x} (see \code{\link{sections}}), or a list with properties \code{intervals} and \code{depths} specifying the intervals of the sections and the corresponding depths respectively.
#' @return Integer vector giving depth of piled intervals from \code{x} (within each sub-interval) or list containing a property \code{"intervals"}, a matrix of sections, and property \code{"depths"}, giving the corresponding pile depths.
#' @export
#' @examples
#' depth(cbind(1:10, 11:20))
depth <- function(x, include_intervals=FALSE) {
	check_ends_at_least_starts(x, stop_if_false=TRUE)
	pts <- breaks(x)
	d <- if (length(pts) == 0L) integer(0) else rcpp_depth(sorted_starts=sort(x[,1L]), sorted_ends=sort(x[,2L]), pts=pts)
	if (include_intervals)
		list(
			intervals=sections(pts),
			depths=d
		)
	else
		d
}

#' @title Flatten a set of intervals
#' @description For a given set of intervals compute the set of intervals where there is overlap with at least one from the given. The resulting intervals are sorted and detached.
#' @param x Integer matrix of two columns, the first column giving the (inclusive) start points of intervals and the second column giving the corresponding (exclusive) end points.
#' @template return-intervals
#' @export
#' @examples
#' flatten(rbind(c(1, 3), c(2, 4), c(5, 6)))
flatten <- function(x) {
	o <- depth(x, include_intervals=TRUE)
	stitch(o$intervals[o$d > 0,,drop=FALSE])
}

#' @title Get IDs of intervals covering each sub-interval
#' @description Get the intervals overlapping each section as a list. 
#' @param x Integer matrix of two columns, the first column giving the (inclusive) start points of intervals and the second column giving the corresponding (exclusive) end points. 
#' @param interval_names Character vector of names for each interval, not necessarily unique. If they are not unique, one might wish to \code{lapply} \code{unique} to the list of members for each sub-interval returned by this function. Defaults to the \code{rownames} of \code{x}.
#' @param output Character value either \code{"list"} or \code{"vector"} determining whether a named list of interval index/name vectors or flat vector of members (corresponding to the output of \code{\link{depth}}) is returned.
#' @return See notes on \code{output} parameter.
#' @export
#' @examples
#' pile(cbind(1:10, 11:20))
pile <- function(x, interval_names=rownames(x), output="list") {
	check_ends_at_least_starts(x, stop_if_false=TRUE)
	if (nrow(x) == 0L) {
		return(if (output == "list") { list() } else { if (is.null(interval_names)) integer(0L) else character(0L) })
	}
	pts <- breaks(x)
	sub_int_starts <- pts[-length(pts)]
	sub_int_ends <- pts[-1L]
	sub_int_names <- if (length(pts) > 1L) paste0("[",sub_int_starts,",",sub_int_ends,")") else character(0)
	if (is.unsorted(x[,1L]))
		stop("Intervals must be sorted by start position!")
	d <- rcpp_depth(sorted_starts=x[,1L], sorted_ends=sort(x[,2L]), pts=pts)
	mem_raw <- rcpp_pile(starts=x[,1L], ends=x[,2L], pts=pts, total_members=sum(d))
	mem <- if (is.null(interval_names)) mem_raw+1L else interval_names[mem_raw+1L]
	if (output=="list")
		split(mem, factor(rep(sub_int_names, d), levels=sub_int_names))
	else if (output=="vector")
		mem
	else
		stop("Invalid 'output' argument")
}

#' @title Get all overlapping tuples of intervals from multiple sets
#' @description Get matrix specifying overlapping tuples of intervals from multiple sets. Each row specifies an overlapping tuple. The \code{n}th element in a row contains the row index of the interval in the \code{n}th set of intervals passed to the function. Depending on the value of the \code{output} argument, there may two additional columns giving the start and end coordinates of the overlap (the default: \code{output="intervals"}, no extra columns (\code{output="indices"}) or one additional column giving the row index of the 'section' of the complete set of intervals (\code{output="sections"}, see \code{\link{sections}}).   
#' @param ... Integer matrices of two columns, the first column giving the (inclusive) start points of intervals and the second column giving the corresponding (exclusive) end points. 
#' @param output Character value, one of \code{"intervals"}, \code{"indices"} and \code{"sections"}. 
#' @return Integer matrix.
#' @export
#' @examples
#' join(rbind(c(1, 100), c(50, 100)), rbind(c(1, 2), c(49, 51), c(50, 200)))
join <- function(..., output="intervals") {
	stopifnot(class(output) == "character")
	valid_output_options <- c("intervals", "indices", "sections")
	if (!any(valid_output_options == output))
		stop(paste0("'output' argument must be one of ", paste0(collapse=", ", "'", valid_output_options, "'")))
	interval_sets <- list(...)
	all_intervals <- do.call(what=rbind, interval_sets)
	check_ends_at_least_starts(all_intervals, stop_if_false=TRUE)
	ord <- order(all_intervals[,1L])
	is_annotation <- factor(rep(seq_along(interval_sets), times=vapply(FUN.VALUE=integer(1), FUN=nrow, X=interval_sets)),levels=seq_along(interval_sets))[ord]
	index_original <- as.integer(do.call(what=c, lapply(interval_sets, function(ints) seq_len(nrow(ints)))))[ord]-1L
	sorted <- all_intervals[ord,,drop=FALSE]
	d <- depth(x=sorted, include_intervals=FALSE)
	tiles <- pile(x=sorted, interval_names=NULL, output="vector")
	set_by_section <- table(set=factor(rep(seq_along(d), times=d), levels=seq_along(d)), section=is_annotation[tiles])
	by_int_set <- split(index_original[tiles], is_annotation[tiles])
	starts <- cumsum(c(0L, vapply(FUN.VALUE=integer(1), FUN=length, X=by_int_set[-length(by_int_set)])))
	vals <- as.integer(do.call(what=c, by_int_set))
	total_overlaps <- sum(apply(set_by_section, 1, prod))
	if (total_overlaps > .Machine$integer.max)
		stop("Too many overlapping tuples of intervals")
	cubes <- rcpp_hyper_cubes(vals, starts, set_by_section) + 1L

	if (output=="sections") {
		cubes
	} else {
		ovrlps <- cubes[,-ncol(cubes),drop=FALSE]
		indices <- ovrlps[!duplicated(ovrlps),,drop=FALSE]
		if (output == "indices") {
			indices	
		} else if (output == "intervals") {
			secs <- sections(breaks(sorted))
			overlap_parts <- unname(split(t(indices), seq_along(interval_sets)))
			bps <- lapply(1:2, function(bp) Map(interval_sets, overlap_parts, f=function(ints, parts) ints[parts,bp]))
			do.call(what=function(...) cbind(indices, ...), Map(list(pmax, pmin), bps, f=function(fun, bp_sets) do.call(what=fun, bp_sets)))
		} else {
			stop("Invalid output option")
		}
	}
}

#' @title Annotate one set of intervals with the names of those which intersect with the other
#' @description Create a list of vectors of indices/names of intervals/points in \code{y} (if \code{y} is a two-column matrix/vector respectively) which intersect with each interval/point in \code{x} (if \code{x} is a two-column matrix/vector respectively). 
#' @param x Integer matrix of two columns, the first column giving the (inclusive) start points of intervals and the second column giving the corresponding (exclusive) end points, or, an integer vector specifying the location of points. 
#' @param annotation Matrix specifying intervals or vector specifying points with which to annotate \code{x}. 
#' @return List of vectors of indices of overlapping intervals/points.
#' @export
#' @examples
#' annotate(rbind(A=c(1, 100), B=c(50, 100)), rbind(a=c(1, 2), b=c(49, 51), c=c(50, 200)))
#' annotate(rbind(A=c(1, 100), B=c(50, 100)), c(a=1, b=49, c=51, d=100))
annotate <- function(x, annotation) {
	x_mat <- if (is.vector(x)) cbind(x, x+1L) else x
	annotation_mat <- if (is.vector(annotation)) cbind(annotation, annotation + 1L) else annotation
	ovlps <- join(x_mat, annotation_mat, output="indices")
	x_rows <- if (is.null(rownames(x_mat))) seq_len(nrow(x_mat)) else rownames(x_mat)
	anno_rows <- if (is.null(rownames(annotation_mat))) seq_len(nrow(annotation_mat)) else rownames(annotation_mat)
	(if (is.null(rownames(x_mat))) unname else identity)(split(f=factor(x_rows[ovlps[,1L]], levels=x_rows), x=anno_rows[ovlps[,2L]]))
}
