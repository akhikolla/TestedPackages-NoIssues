#' #### OurICsurvPHgrad calculates the gradient of the survival likleihood under the I-spline model of Wang et al. (2016) 
#' # The code heavily uses code from fast.PH.ICsurv.EM function in ICsurv. See also  OurICsurvPH
#' #' @title FUNCTION_TITLE
#' #' @description FUNCTION_DESCRIPTION
#' #' @param Li PARAM_DESCRIPTION
#' #' @param Ri PARAM_DESCRIPTION
#' #' @param Xp PARAM_DESCRIPTION
#' #' @param n.int PARAM_DESCRIPTION
#' #' @param order PARAM_DESCRIPTION
#' #' @param g1 PARAM_DESCRIPTION
#' #' @param b1 PARAM_DESCRIPTION
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
#' # @rdname OurICsurvPHgrad
#' OurICsurvPHgrad <- function (Li, Ri, Xp, n.int, order, g1, b1, tol, t.seq, equal = TRUE) 
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
#'     knots <- c(ti.min, quantile(ti, id), ti.max)
#'   }
#'   bRi <- t(Ispline(x = Ri, order = order, knots = knots))
#'   bLi <- t(Ispline(x = Li, order = order, knots = knots))
#'   bt <- t(Ispline(x = t.seq, order = order, knots = knots))
#'  
#'   GRi <- bRi %*% g1
#'   GLi <- bLi %*% g1
#'   xb <- Xp %*% b1
#'   LL <- exp(-GLi * exp(xb)) - exp(-GRi * exp(xb))
#'     ll <- log(LL)
#' 
#'   u.b <- as.vector((1/LL)*(exp(-GLi * exp(xb)) * ( -GLi * exp(xb)) - exp(-GRi * exp(xb)) * (-GRi * exp(xb)))) * Xp
#'   u.g <- as.vector((1/LL)*(exp(-GLi * exp(xb)) * ( -GLi * exp(xb)) ))* bLi - 
#'                      as.vector((1/LL)*(exp(-GRi * exp(xb)) * (-GRi * exp(xb))))*bRi
#'   u <- cbind(u.b, u.g)
#'     return(u)
#' }
#' 
