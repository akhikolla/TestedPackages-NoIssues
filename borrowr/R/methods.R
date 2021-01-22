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


print.pate <- function(x, ...) {
  nc <- x$non_compliance
  nc <- ifelse(nc, "Yes", "No")
  cat("Population Average Treatment Effect (PATE)\n")
  cat("\nAdjustment for confounding due to noncompliance: ",nc,"\n")
  cat("\nPATE Posterior Summary Statistics (Treated minus Untreated):\n\n")
  # posterior statistics
  post <- x$pate_post
  psum <- posterior_summary(post)
  print(psum)
  out <- list(posterior_summary = psum)
  cat("\n")
  invisible(out)
}

summary.pate <- function(object, ...) {
  x <- object
  out <- list()
  out$call <- x$call
  out$posterior_summary <- posterior_summary(x$pate_post)
  out$mems <- x$MEMs + 0
  out$exch_prob <- x$exch_prob
  out$post_probs <- x$post_probs
  out$non_compliance <- x$non_compliance
  class(out) <- "summary.pate"
  out
}

print.summary.pate <- function(x, ...) {
  cat("\nPopulation Average Treatment Effect (PATE)\n")
  cat("\nCall:\n\n")
  print(x$call)
  nc <- ifelse(x$non_compliance, "Yes", "No")
  cat("\nAdjustment for confounding due to noncompliance: ",nc,"\n")
  cat("\nPATE Posterior Summary Statistics (Treated minus Untreated):\n\n")
  # posterior statistics
  post <- x$pate_post
  psum <- x$posterior_summary
  print(psum)

  cat("\nPrior Probability of Exchangeability:\n\n")
  print(x$exch_prob)
  cat("\nExchangeability Matrix (1 == Exchangeable with primary source):\n\n")
  print(x$mems)
  cat("\nMEM Posterior Probability:\n")
  print(x$post_probs)
  cat("\n")
  invisible(x)
}

