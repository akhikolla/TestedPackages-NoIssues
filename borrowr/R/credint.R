#  borrowr: estimate population average treatment effects with borrowing between data sources.
#  Copyright (C) 2019  Jeffrey A. Verdoliva Boatman
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.


#'Posterior Credible Interval for Population Average Treatment Effect (PATE)
#'
#'Computes an equal-tailed credible interval for an object of class 'pate'
#'@param object An object of class pate fit by the \code{pate} function.
#'@param level the credible level required
#'
#'@examples
#'data(adapt)
#'
#'est <- pate(y ~ treatment*x + treatment*I(x ^ 2), data = adapt,
#'  estimator = "bayesian_lm", src_var = "source", primary_source = "Primary",
#'  trt_var = "treatment")
#'
#'credint(est)
#'
credint <- function(object, level = 0.95){
  if(class(object) != 'pate')
    stop("'object' must be of class 'pate'")
  cf <- coef(object)
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  quants <- quantile(object$pate_post, a)
  ci <- array(quants, dim = c(1L, 2L), dimnames = list("PATE",
    names(quants)))
  ci
}

