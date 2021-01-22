#' #### OurICsurvPH is the same as fast.PH.ICsurv.EM in ICsurv with a couple of distinctions:
#' #### 1)
#' # All observations are assumed to be interval-censored, hence d1,d2,d3 are not needed.
#' # This is done for consistency with the other methods where zeros are replaced with a very small number and Inf 
#' # are replaced with 200 (which did not effected the other methods when I took 2000 or 20000)
#' #### 
#' #' @title FUNCTION_TITLE
#' #' @description FUNCTION_DESCRIPTION
#' #' @param d1 PARAM_DESCRIPTION
#' #' @param d2 PARAM_DESCRIPTION
#' #' @param d3 PARAM_DESCRIPTION
#' #' @param Li PARAM_DESCRIPTION
#' #' @param Ri PARAM_DESCRIPTION
#' #' @param Xp PARAM_DESCRIPTION
#' #' @param n.int PARAM_DESCRIPTION
#' #' @param order PARAM_DESCRIPTION
#' #' @param g0 PARAM_DESCRIPTION
#' #' @param b0 PARAM_DESCRIPTION
#' #' @param tol PARAM_DESCRIPTION
#' #' @param t.seq PARAM_DESCRIPTION
#' #' @param equal PARAM_DESCRIPTION, Default: TRUE
#' #' @return OUTPUT_DESCRIPTION
#' #' @details DETAILS
#' #' @examples 
#' #' \dontrun{
#' #' if(interactive()){
#' #'  #EXAMPLE1
#' #'  }
#' #' }
#' #' @seealso 
#' #'  \code{\link[ICsurv]{Ispline}}
#' # @rdname OurICsurvPH
#' #' @importFrom ICsurv Ispline
#' OurICsurvPH <- function (d1, d2, d3, Li, Ri, Xp, n.int, order, g0, b0, tol, 
#'           t.seq, equal = TRUE) 
#' {
#'   P <- length(b0)
#'   L <- length(g0)
#'   N <- length(Li)
#'   ti <- c(Li, Ri)
#'   if (equal == TRUE) {
#'     ti.max <- max(ti) + 1e-05
#'     ti.min <- min(ti) - 1e-05
#'     knots <- seq(ti.min, ti.max, length.out = (n.int + 2))
#'   }
#'   if (equal == FALSE) {
#'     id <- seq(0, 1, length.out = (n.int + 2))
#'     id <- id[-c(1, (n.int + 2))]
#'     ti.max <- max(ti) + 1e-05
#'     ti.min <- min(ti) - 1e-05
#'     knots <- c(ti.min, stats::quantile(ti, id), ti.max)
#'   }
#'   bRi <- t(ICsurv::Ispline(x = Ri, order = order, knots = knots))
#'   bLi <- t(ICsurv::Ispline(x = Li, order = order, knots = knots))
#'   bt <- t(ICsurv::Ispline(x = t.seq, order = order, knots = knots))
#'   Q1 <- function(b0, b1, g0, Xp, bRi, bLi, d1, d2, d3, L) {
#'     g1 <- rep(-99, L)
#'     xb0 <- Xp %*% b0
#'     xb1 <- Xp %*% b1
#'     dz <- 1 - exp(-(bRi %*% g0) * exp(xb0))
#'     dw <- 1 - exp(-(bRi %*% g0 - bLi %*% g0) * exp(xb0))
#'     dw[d2 == 0] = 1
#'     EZil <- t(t(d1 * bRi) * g0) * as.vector(exp(xb0)/dz)
#'     EWil <- t(t(d2 * (bRi - bLi)) * g0) * as.vector(exp(xb0)/dw)
#'     num <- EZil + (d2 + d3) * EWil
#'     den <- ((d1 + d2) * bRi + d3 * bLi) * as.vector(exp(xb1))
#'     g1 <- apply(num, 2, sum)/apply(den, 2, sum)
#'     return(g1)
#'   }
#'   Q2 <- function(b1, b0, g0, Xp, bRi, bLi, d1, d2, d3, L) {
#'     xb0 <- Xp %*% b0
#'     xb1 <- Xp %*% b1
#'     dz <- 1 - exp(-(bRi %*% g0) * exp(xb0))
#'     dw <- 1 - exp(-(bRi %*% g0 - bLi %*% g0) * exp(xb0))
#'     dw[d2 == 0] = 1
#'     EZi <- d1 * (bRi %*% g0) * exp(xb0)/dz
#'     EWi <- d2 * (bRi %*% g0 - bLi %*% g0) * exp(Xp %*% b0)/dw
#'     p1 <- sum((EZi + EWi) * (Xp %*% b1))
#'     EZil <- t(t(d1 * bRi) * g0) * as.vector(exp(xb0)/dz)
#'     EWil <- t(t(d2 * (bRi - bLi)) * g0) * as.vector(exp(xb0)/dw)
#'     num <- EZil + (d2 + d3) * EWil
#'     den <- ((d1 + d2) * bRi + d3 * bLi) * as.vector(exp(xb1))
#'     g1 <- apply(num, 2, sum)/apply(den, 2, sum)
#'     p2 <- sum(t(EZil) * log(g1) + t(EWil) * log(g1))
#'     p3 <- sum(((d1 + d2) * (bRi %*% g1) + d3 * (bLi %*% g1)) * 
#'                 exp(xb1))
#'     res <- -(p1 + p2 - p3)
#'     return(res)
#'   }
#'   b1 <- stats::optim(b0, Q2, method = "Nelder-Mead", b0 = b0, g0 = g0, 
#'               Xp = Xp, bRi = bRi, bLi = bLi, d1 = d1, d2 = d2, d3 = d3, 
#'               L = L)$par
#'   g1 <- Q1(b0, b1, g0, Xp, bRi, bLi, d1, d2, d3, L)
#'   while (max(abs(c(b0, g0) - c(b1, g1))) > tol) {
#'     b0 <- b1
#'     g0 <- g1
#'     b1 <- stats::optim(b0, Q2, method = "Nelder-Mead", b0 = b0, 
#'                 g0 = g0, Xp = Xp, bRi = bRi, bLi = bLi, d1 = d1, 
#'                 d2 = d2, d3 = d3, L = L)$par
#'     g1 <- Q1(b0, b1, g0, Xp, bRi, bLi, d1, d2, d3, L)
#'   }
#'   GRi <- bRi %*% g1
#'   GLi <- bLi %*% g1
#'   xb <- Xp %*% b1
#'   ll <- sum(log(exp(-GLi * exp(xb)) - exp(-GRi * exp(xb))))
#'    AIC <- 2 * (length(b1) + length(g1)) - 2 * ll
#'   BIC <- (length(b1) + length(g1)) * log(N) - 2 * ll
#'   v <- ICsurv::fast.PH.Louis.ICsurv(b1, g1, bLi, bRi, d1, d2, d3, Xp)
#'   Hessian = v
#'  # flag = is.non.singular.matrix(v)
#'   flag <- class(try(solve(v),silent=T))=="matrix"
#'   if (flag) {
#'     var.b = solve(v)[1:P, 1:P]
#'   }
#'   else {
#'     A = v[1:P, 1:P]
#'     B = v[1:P, (P + 1):(P + L)]
#'     C = v[(P + 1):(P + L), 1:P]
#'     D = v[(P + 1):(P + L), (P + 1):(P + L)]
#'     var.b = ginv(A - B %*% ginv(D) %*% C)
#'   }
#'   hz <- bt %*% g1
#'   return(list(b = b1, g = g1, hz = hz, Hessian = Hessian, var.b = var.b, 
#'               flag = flag, ll = ll, AIC = AIC, BIC = BIC))
#' }
