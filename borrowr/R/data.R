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


#' Data set used in the package vignette
#'
#' A simulated data set used in the package vignette. Data generating
#' mechanism is adapted from Hill (2011). Includes 3 data sources to illustrate
#' borrowing information in estimate the population average treatment effect.
#'
#' @format A data frame with 180 rows and 5 variables
#' \describe{
#'   \item{y}{the outcome}
#'   \item{x}{a variable associated with y and treatment}
#'   \item{source}{character variable with 3 levels: "Primary", "Supp1", "Supp2"}
#'   \item{treatment}{treatment variable. 1 = treated, 0 = untreated}
#'   \item{compliant}{compliance indicator. 1 = compliant to treatment,
#'   0 = noncompliant to treatment}
#' }
#' @examples
#' data(adapt)
#' head(adapt)
#' @references
#' Hill, Jennifer L. (2011) Bayesian Nonparametric Modeling for Causal Inference.
#' Journal of Computational and Graphical Statistics, 20:1, 217-240.
"adapt"
