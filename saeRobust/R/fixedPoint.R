#' Fixed Point Algorithm Infrastructure
#'
#' @description A fixed-point function supplied by the user is iteratively
#'   evaluated. \code{addAverageDamp} can be used to add average damping to the
#'   function - this may have a positive effect on the speed of convergence.
#'
#' @param fun the function to be evaluated in the algorithm
#' @param x0 starting value
#' @param convCrit a function returning a logical scalar. Is called with two
#'   arguments; the first is the value from iteration n; the second is the value
#'   from iteration n-1
#' @param tolerance a numeric value > 0
#' @param value (numeric)
#'
#' @export
#' @rdname fixedPoint
#' @examples
#' \dontrun{
#' vignette("fixedPoint", "saeRobust")
#' }
fixedPoint <- function(fun, x0, convCrit) {
    assert_that(is.function(fun))
    assert_that(is.function(convCrit))
    x1 <- NULL
    repeat {
        x1 <- fun(x0)
        if (convCrit(x1, x0)) {
            break
        } else {
            x0 <- x1
        }
    }
    x1
}

#' @details \code{addAverageDamp} adds average damping to an arbitrary fixed point
#'   function.
#' @export
#' @rdname fixedPoint
addAverageDamp <- function(fun) {
    assert_that(is.function(fun))
    function(x) (x + fun(x)) / 2
}

#' @details \code{addConstraintMin} takes care that values are not below a
#'   minimum value.
#' @export
#' @rdname fixedPoint
addConstraintMin <- function(fun, value) {
    assert_that(is.function(fun))
    function(x) pmax(value, fun(x))
}

#' @details \code{addConstraintMax} takes care that values are not larger than
#'   maximum value.
#' @export
#' @rdname fixedPoint
addConstraintMax <- function(fun, value) {
    assert_that(is.function(fun))
    function(x) pmin(value, fun(x))
}

#' @details \code{convCritAbsolute} absolute difference as convergence criterion.
#' @export
#' @rdname fixedPoint
convCritAbsolute <- function(tolerance = 1e-6) {
    assert_that(tolerance > 0)
    function(xn1, xn0) log(1 + mean(abs(xn0 - xn1))) < tolerance
}

#' @details \code{convCritRelative} relative (to previous iteration) absolute
#'   difference as convergence criterion. If value is smaller than 1, absolute
#'   difference is used.
#' @export
#' @rdname fixedPoint
convCritRelative <- function(tolerance = 1e-6) {
    assert_that(tolerance > 0)
    function(xn1, xn0) log(1 + mean(abs(xn0 - xn1) / max(1, abs(xn0)))) < tolerance
}

#' @details \code{addMaxIter} can be used to modify convergence criterion functions.
#'
#' @param maxIter maximum number of iterations
#'
#' @rdname fixedPoint
#'
#' @export
addMaxIter <- function(fun, maxIter) {
    assert_that(is.function(fun))
    assert_that(maxIter > 0)
    count <- 1
    function(...) {
        count <<- count + 1
        if (count > maxIter) TRUE else fun(...)
    }
}

#' @details \code{addCounter} can be used to count the number of calls of a function.
#'
#' @export
#' @rdname fixedPoint
addCounter <- function(fun) {
    assert_that(is.function(fun))
    count <- 0
    function(...) {
        count <<- count + 1
        addAttr(fun(...), count, "count")
    }
}

#' @details \code{addHistory} can be used to save a history of results of a
#'   function. The history is stored as a matrix, so this works best if the
#'   return value of \code{fun} is numeric.
#'
#' @export
#' @rdname fixedPoint
addHistory <- function(fun) {
    force(fun)
    history <- NULL
    function(...) {
        res <- fun(...)
        history <<- rbind(history, res)
        addAttr(res, history, "history")
    }
}

#' @details \code{addStorage} will add a storage to a function. The storage is a
#'   list in which each result is stored. The function will coerce the return
#'   value into a numeric.
#'
#' @export
#' @rdname fixedPoint
addStorage <- function(fun) {
    force(fun)
    storage <- list()
    function(...) {
        res <- fun(...)
        storage <<- c(storage, list(res))
        addAttr(as.numeric(unlist(res)), storage, "storage")
    }
}

#' @details \code{newtonRaphson} finds zeroes of a function. The user can supply
#'   the function and its first derivative. Note that the Newton Raphson
#'   Algorithm is a special case of a fixed point algorithm thus it is
#'   implemented using \code{\link{fixedPoint}} and is only a convenience.
#'
#' @param funList (list) the functions to be evaluated in the algorithm. First
#'   element is typically the score function, second is the derivative of the
#'   score.
#' @param ... arguments passed to \code{\link{fixedPoint}}
#'
#' @export
#' @rdname fixedPoint
newtonRaphson <- function(funList, ...) {
    fixedPoint(newtonRaphsonFunction(funList), ...)
}

#' @export
#' @rdname fixedPoint
newtonRaphsonFunction <- function(funList) {
    force(funList)
    function(x) as.numeric(x - solve(funList$f1(x)) %*% funList$f(x))
}
