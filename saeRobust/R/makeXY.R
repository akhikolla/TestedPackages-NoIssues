#' makeXY
#'
#' @description Extract respone vector and design matrix from data with given
#'   formula.
#'
#' @param .formula (formula)
#' @param .data (data.frame) data from which design matrix and response to
#'   extract from
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' makeXY(y ~ I(x^2), data.frame(x = rnorm(10), y = rnorm(10)))
makeXY <- function(.formula, .data){

    assert_that(inherits(.formula, "formula"), is.data.frame(.data))
    .mf <- model.frame(.formula, .data)
    x <- Matrix(model.matrix(attr(.mf, "terms"), data = .mf))
    y <- model.response(.mf)

    retList() %>% stripSelf

}
