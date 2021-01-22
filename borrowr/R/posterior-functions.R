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


#Sample from the weighted posterior.
#
#Given a matrix with samples from the MEM posteriors and a correspodning vector
#of log marginal likelihoods, returns a sample from the weighted posterior.
#@param posts An \emph{n} by \emph{p} matrix with \emph{n} samples from
#             the posterior of each of the \emph{p} MEMs.
#@param log_marg_like A vector of length \emph{p} with the log marginal marginal likelihoods
#for each of the \emph{p} MEMs.
#@return A numeric vector of length \emph{n} containing a
#sample from the weighted posterior.
sample_posterior <- function(posts, probs) {

  n_post <- nrow(posts)
  n_mems <- ncol(posts)

  row_select <- seq_len(n_post)
  col_select <- sample(seq_len(n_mems), n_post, replace = TRUE, prob = probs)

  out <- posts[cbind(row_select, col_select)]
  out

}

posterior_summary <- function(x) {
  c(
    "Mean Treatment Effect" = mean(x),
    "Std. Dev."   = sd(x),
    "Pr(PATE > 0)" = mean(x > 0)
  )
}

#For a numeric vector \code{x}, the function returns a single Monte Carlo
#draw from the Bayesian bootstrap distribution of the posterior mean.
#@param x A numeric vector.
#@return A numeric value: a single Monte Carlo draw from the Bayesian bootstrap
#distribution of the posterior mean
#@references
#Rubin, D.B. (1981) The Bayesian Bootstrap. The Annals of Statistics; 9(1): 130-134.
#
bayes_boot_mean <- function(x) {
  n <- length(x)
  u <- sort(runif(n - 1))
  g <- diff(c(0, u, 1))
  sum(g * x)
}

