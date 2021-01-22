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


make_prior_trees <- function(nprior, ntree, cut_lens, var_probs = NULL, depth = 8, alpha = 0.95, beta = 2, sigma_mu) {

  p <- length(cut_lens)
  # if(missing(var_probs))
  #   var_probs <- rep(1 / p, p)
  if(is.null(var_probs))
    var_probs <- rep(1 / p, p)
  sigma_mu_list <- as.list(rep(sigma_mu, each = ntree))

  # debugonce(make_one_prior_treeC)
  # tree_list <- vapply(1:(nprior * ntree), make_one_prior_treeC, character(1),
  #   cut_lens  = cut_lens,
  #   var_probs = var_probs,
  #   sigma_mu  = 0.5 / (2 * sqrt(ntree)),
  #   #sigma_mu = 1 / 4,
  #   tiers     = depth,
  #   beta      = beta
  #   )

  # original implementation:
  # tree_list <- mapply(make_one_prior_tree,
  #   # seed      = 1:(nprior * ntree),
  #   cut_lens  = list(cut_lens),
  #   var_probs = list(var_probs),
  #   sigma_mu  = sigma_mu_list,
  #   tiers     = depth,
  #   beta      = beta
  # )

  # rewrite as for-loop
  nt <- nprior * ntree
  tree_list <- character(nt)
  # tree_list <- as.list(tree_list)
  for(ii in seq_len(nt)) {
    tree_list[ii] <- make_one_prior_tree(cut_lens = cut_lens,
      var_probs = var_probs, sigma_mu = sigma_mu_list[[ii]],
      tiers = depth, beta = beta)
  }

  stream <- paste(tree_list, collapse = "\n")
  line_1 <- paste0(paste(nprior, ntree, p, collapse = " "), "\n")
  stream <- paste0(line_1, stream)

  stream

}

make_one_prior_tree <- function(cut_lens, var_probs, sigma_mu, tiers = 8, alpha = 0.95, beta = 2) {

  # set.seed(seed)

  # res <- priortree(cut_lens = cut_lens,
  #   cum_sum_var_probs = cumsum(var_probs),
  #   sigma_mu          = sigma_mu,
  #   tiers             = tiers,
  #   alpha             = alpha,
  #   beta              = beta)


  # nodes  <- res$nodes
  # killed <- res$killed
  # var    <- res$var
  # cut    <- res$cut
  # leaves <- res$leaves
  #
  # nodes  <- nodes[!killed]
  # var    <- var[!killed]
  # cut    <- cut[!killed]
  # leaves <- leaves[!killed]

  # nn <- length(nodes)
  # stream    <- character(2)
  # stream[1] <- paste0(nn, "\n")
  # stream[2] <- paste(nodes, var, cut, leaves, collapse = "\n")

  res <- priortree(cut_lens = cut_lens, alpha = alpha, beta = beta)

  nn <- length(res$node)
  stream    <- character(2)
  stream[1] <- paste0(nn, "\n")
  stream[2] <- paste(res$node, res$var, res$cut, res$leave, collapse = "\n")

  out <- paste(stream, collapse = "")
  out
}
