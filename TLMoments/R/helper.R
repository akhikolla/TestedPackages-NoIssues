# automatically select computation.method:
#   L-moments: recursive
#   TL-moments: recurrence
#   NI-TL-moments: direct
select_computation <- function(leftrim, rightrim) {
  stopifnot(is.numeric(leftrim), is.numeric(rightrim))

  if (isTRUE(all.equal(leftrim, 0L)) & isTRUE(all.equal(rightrim, 0L))) {
    return("recursive")
  } else if (are.integer.like(leftrim, rightrim)) {
    return("recurrence")
  } else {
    return("direct")
  }
}

# calculate pseudo observations
pseudo <- function(x, k) {
  as.matrix(pseudo_C(sort(stats::na.omit(x)), k)[rank(x, na.last = "keep"), ])
}

# Functions for general error catching
are.integer.like <- function(...) {
  all(vapply(list(...), function(a) isTRUE(all.equal(a, as.integer(a))), logical(1)))
}
are.numeric <- function(...) {
  all(vapply(list(...), function(a) is.numeric(a), logical(1)))
}

`%in.equal%` <- function(a, b) {
  r <- vapply(a, function(.a) vapply(b, function(.b) { isTRUE(all.equal(.a, .b)) }, logical(1)), logical(length(b)))
  if (is.null(dim(r))) {
    r
  } else {
    apply(r, 2, any)
  }
}

# infix for intervals
`%-+%` <- function(a, b) { cbind(a - b, a + b) }

# @description Calculates TL-Moments for every distribution in distr for the specified trim.
# Used by plot.TLMoments
# @return data.frame
# @examples
# getTLMomsByDistr(c("gev", "gpd", "gum", "norm"), c(2, 1))
getTLMomsByDistr <- function(distr, trim) {
  # "Lines" & "Points"
  .lines <- function(distr, shapes, trim) {
    sapply(shapes, function(shape) {
      r <- tryCatch(
        calcTLMom(4, trim[1], trim[2], qfunc = getQ.character(distr, loc = 0, scale = 1, shape = shape)),
        error = function(e) rep(NA, 4)
      )
      setNames(r, paste0("L", 1:4))
    })
  }
  .points  <- function(distr, trim) tryCatch(calcTLMom(4, trim[1], trim[2], qfunc = getQ.character(distr)), error = function(e) rep(NA, 4))

  if ("gum" %in% distr) {
    gum <- data.frame(
      distr = "gum",
      leftrim = trim[1],
      rightrim = trim[2],
      shape = NA,
      data.frame(t(setNames(.points("gum", trim), paste0("L", 1:4))))
    )
  } else gum <- NULL

  if ("exp" %in% distr) {
    exp <- data.frame(
      distr = "exp",
      leftrim = trim[1],
      rightrim = trim[2],
      shape = NA,
      data.frame(t(setNames(.points("exp", trim), paste0("L", 1:4))))
    )
  } else exp <- NULL

  if ("norm" %in% distr) {
    norm <- data.frame(
      distr = "norm",
      leftrim = trim[1],
      rightrim = trim[2],
      shape = NA,
      data.frame(t(setNames(.points("norm", trim), paste0("L", 1:4))))
    )
  } else norm <- NULL

  if ("gev" %in% distr) {
    shapes <- seq(-.99, .99, .01)
    gev <- data.frame(
      distr = "gev",
      leftrim = trim[1],
      rightrim = trim[2],
      cbind(
        shape = shapes,
        t(.lines("gev", shapes, trim))
      )
    )
  } else gev <- NULL

  if ("gpd" %in% distr) {
    shapes <- seq(-.99, .99, .01)
    gpd <- data.frame(
      distr = "gpd",
      leftrim = trim[1],
      rightrim = trim[2],
      cbind(
        shape = shapes,
        t(.lines("gpd", shapes, trim))
      )
    )
  } else gpd <- NULL

  if ("ln3" %in% distr) {
    shapes <- seq(-.99, .99, .01)
    ln3 <- data.frame(
      distr = "ln3",
      leftrim = trim[1],
      rightrim = trim[2],
      cbind(
        shape = shapes,
        t(.lines("ln3", shapes, trim))
      )
    )
  } else ln3 <- NULL

  dat <- rbind(gum, exp, norm, gev, gpd, ln3)
  dat$T3 <- dat$L3 / dat$L2
  dat$T4 <- dat$L4 / dat$L2
  dat <- dat[!(is.na(dat$T3) | is.na(dat$T4)), ]
  dat
}


# @description calculates the first maxr TL-moments of order s, t given a quantile function qfunc
# @examples
# calcTLMom(4, 0, 2, function(u) evd::qgumbel(u, loc = 0, scale = 1))
# lmomco::theoTLmoms(lmomco::vec2par(c(0, 1), type = "gum"), nmom = 4, leftrim = 0, rightrim = 2)$lambdas
#
# calcTLMom(4, 2, 1, function(u) evd::qgev(u, loc = 0, scale = 1, shape = .2))
# lmomco::theoTLmoms(lmomco::vec2par(c(0, 1, -.2), type = "gev"), nmom = 4, leftrim = 2, rightrim = 1)$lambdas
#
# microbenchmark::microbenchmark(
#   calcTLMom(4, 2, 1, function(u) evd::qgev(u, loc = 0, scale = 1, shape = .2)),
#   lmomco::theoTLmoms(lmomco::vec2par(c(0, 1, -.2), type = "gev"), nmom = 4, leftrim = 2, rightrim = 1)$lambdas
# )
calcTLMom <- function(maxr, s, t, qfunc, ...) {
  #if (!are.integer.like(maxr, s, t)) stop("s, t, and maxr have to be integer-like. ")
  if (s < 0 | t < 0 | maxr < 1) stop("s and t have to be >= 0, maxr has to be >= 1. ")
  if (!is.function(qfunc)) stop("qfunc has to be a function. ")

  vapply(1L:maxr, function(r) {
    sum(vapply(0:(r-1), function(j) {
      tryCatch(
        i <- stats::integrate(f <- function(u) u^(s+j) * (1-u)^(t+r-j-1) * qfunc(u, ...), lower = 0, upper = 1),
        error = function(e) stop(e)
      )
      if (i$message != "OK") stop("Error occurred while integrating. ")
      (-1)^(r-j-1) * factorial(r-1)*factorial(r+s+t) /r /factorial(j) /factorial(r+t-j-1) /factorial(r-j-1) /factorial(s+j) * i$value
    }, numeric(1)))
  }, numeric(1))
}

# @param args vector of characters giving the argument names to check
# @param distr character
checkParameterNames <- function(args, distr) {
  stopifnot(length(distr) == 1)
  stopifnot(is.character(args))

  if (grepl("::", x = distr)) { # if pkg::func
    f <- sub("^([a-zA-Z0-9]*)::([a-zA-Z0-9]*)$", "\\1::q\\2", x = distr)
    q <- eval(parse(text = paste0("match.fun(", f,")")))
  } else { # falls nur func
    q <- eval(parse(text = paste0("match.fun(q", distr, ")")))
  }
  .args <- names(formals(q))[sapply(formals(q), is.numeric)]

  args <- .args[pmatch(args, .args, nomatch = NA)]
  if (any(is.na(args))) stop("Names do not match distribution arguments or are ambiguous. ")

  args
}

# @description extracts the quantile function of a parameters-object or a character string (like evd::gev)
# @examples
# getQ(as.parameters(loc = 9, scale = 5, shape = .3, distr = "evd::gev"))
# getQ("evd::gev", loc = 10, scale = 4, shape = .2)
getQ <- function(x, ...) {
  if (!inherits(x, c("parameters", "character")))
    stop("x must be of class parameters or character vector")

  UseMethod("getQ")
}
getQ.character <- function(x, ...) {
  distr <- x
  args <- list(...)

  if (grepl("::", x = distr)) { # if pkg::func
    f <- sub("^([a-zA-Z0-9]*)::([a-zA-Z0-9]*)$", "\\1::q\\2", x = distr)
    q <- eval(parse(text = paste0("match.fun(", f,")")))
  } else { # falls nur func
    q <- eval(parse(text = paste0("match.fun(q", distr, ")")))
  }
  if (!is.function(q)) stop(paste0("Found no q-function for ", distr))

  if (any(!(names(args) %in% names(formals(q))))) stop("Wrong arguments given.")
  formals(q)[names(args)] <- args
  q
}
getQ.parameters <- function(x) {
  if (!inherits(x, "numeric")) stop("By now only for parameters, numeric!")

  distr <- attr(x, "distribution")
  args <- as.list(x)

  do.call(getQ.character, c(x = distr, args))
}

correctNames <- function(x, forbidden_pattern, preceding) {

  if (inherits(x, "data.frame")) {
    forbids <- grepl(forbidden_pattern, names(x))
    if (any(forbids)) {
      forbidden_names <- names(x)[forbids]
      new_names <- paste0(preceding, forbidden_names)
      names(x)[forbids] <- new_names
      warning("Renamed variables due to invalid names. ")
    }

  } else if (inherits(x, "formula")) {
    forbids <- grepl(forbidden_pattern, all.vars(x))
    if (any(forbids)) {
      forbidden_names <- all.vars(x)[forbids]
      strformula <- deparse(x)
      for (i in seq_along(forbidden_names)) {
        strformula <- sub(forbidden_names[i], paste0(preceding, forbidden_names[i]), strformula)
      }
      x <- as.formula(strformula)
      warning("Renamed variables due to invalid names. ")
    }

  }

  x
}

calcRatios <- function(lambdas) {
  if (length(lambdas) > 2L) {
    out <- c(NA, lambdas[2]/lambdas[1], lambdas[3:length(lambdas)]/lambdas[2])
  } else if (length(lambdas) == 2L) {
    out <- c(NA, lambdas[2]/lambdas[1])
  } else {
    out <- NA
  }
  setNames(out, paste0("T", seq_along(out)))
}

calcLambdas <- function(ratios, L1) {
  if (is.na(L1)) {
    L1 <- 1
    warning("L1 is assumed to be 1. ")
  }

  if (length(ratios) > 1L) {
    out <- c(L1, ratios[1]*L1, ratios[2:length(ratios)]*(ratios[1]*L1))
  } else if (length(ratios) == 1L) {
    out <- c(L1, ratios[1]*L1)
  } else {
    out <- L1
  }
  setNames(out, paste0("L", seq_along(out)))
}

getFormulaSides <- function(formula, names = NULL) {

  lhs <- all.vars(update(formula, .~0))
  all <- all.vars(formula)
  rhs <- all[!(all %in% lhs)]

  if (!is.null(names)) {
    if (length(lhs) == 1 && lhs == ".") {
      lhs <- names[!(names %in% rhs)]
    }
    if (length(rhs) == 1 && rhs == ".") {
      rhs <- names[!(names %in% lhs)]
    }
    if ("." %in% all) all <- names

    #if (!all(lhs %in% names)) stop("Formula error")
    #if (!all(rhs %in% names)) stop("Formula error")
  }

  list(lhs = lhs, rhs = rhs, all = all,
       new.formula = stats::as.formula(paste0("cbind(", paste0(lhs, collapse = ","), ") ~ ", paste0(rhs, collapse = "+"))))
}
# getFormulaSides(z ~ x + y)
# getFormulaSides(cbind(z1, z2) ~ x + y)
# getFormulaSides(. ~ x + y)
# getFormulaSides(cbind(z1, z2) ~ .)
# getFormulaSides(. ~ x + y, names = c("z1", "z2", "x", "y"))
# getFormulaSides(cbind(z1, z2) ~ ., names = c("z1", "z2", "x", "y"))


blockdiag <- function(x, j, back = NULL) {
  d <- dim(as.matrix(x))
  if (!is.null(back) & (dim(back)[1] != d[1]*j || dim(back)[2] != d[2]*j)) {
    warning("Wrong dimensions of background matrix. Set to Zero-Matrix. ")
    back <- NULL
  }
  if (is.null(back)) {
    X <- matrix(0, nrow = d[1] * j, ncol = d[2] * j)
  } else {
    X <- back
  }
  for (i in 0:(j-1)) {
    X[(i*d[1]+1):((i+1)*d[1]), (i*d[2]+1):((i+1)*d[2])] <- x
  }
  X
}
# blockdiag(matrix(1:4, 2), 3)
# blockdiag(matrix(1:8, 2), 3)
# blockdiag(matrix(1:4, 2), 3, matrix(NA, nr = 6, nc = 6))
# blockdiag(matrix(1:4, 2), 3, matrix(NA, nr = 6, nc = 7))
blockdiag_list <- function(x, back = NULL) {
  dims <- lapply(x, dim)
  dim <- c(sum(sapply(dims, getElement, 1)), sum(sapply(dims, getElement, 2)))

  if (!is.null(back) & (dim(back)[1] != dim[1] || dim(back)[2] != dim[2])) {
    warning("Wrong dimensions of background matrix. Set to Zero-Matrix. ")
    back <- NULL
  }
  if (is.null(back)) {
    X <- matrix(0, nrow = dim[1], ncol = dim[2])
  } else {
    X <- back
  }
  pos <- list(
    c(0, cumsum(sapply(dims, getElement, 1))),
    c(0, cumsum(sapply(dims, getElement, 2)))
  )
  for (i in 1:length(x)) {
    X[(pos[[1]][i]+1):pos[[1]][i+1], (pos[[2]][i]+1):pos[[2]][i+1]] <- x[[i]]
  }
  X
}
# blockdiag_list(list(matrix(1:4, 2), matrix(1:9, 3)))
# blockdiag_list(list(matrix(1:4, 2), matrix(1:8, 2)), back = matrix(NA, nr = 4, nc = 6))


removeAttributes <- function(x) {
  attr(x, "source") <- NULL
  attr(x, "class") <- NULL
  attr(x, "order") <- NULL
  attr(x, "distribution") <- NULL
  attr(x, "leftrim") <- NULL
  attr(x, "rightrim") <- NULL
  attr(x, "computation.method") <- NULL
  x
}


buildNames <- function(prefix, order, stations = NULL) {
  if (is.null(stations) || length(stations) == 1) {
    paste0(prefix, order)
  } else {
    paste0(rep(paste0(prefix, order), length(stations)), "_", rep(stations, each = length(order)))
  }
}
