##' Compute an optimal segmentation (change in Gaussian mean model,
##' square loss), which is consistent with the given labels, and with
##' a penalty for each changepoint outside of labeled regions.
##'
##' Provides a high-level interface to LOPART_interface R function and
##' LOPART C code.
##' @title Labeled Optimal PARTitioning
##' @return list with named elements, all of which are data
##'   tables. loss has one row with loss/cost values. cost is the
##'   output from LOPART_interface. changes has one row for each
##'   predicted changepoint (e.g. change=1.5 means a change between
##'   data points 1 and 2). segments has one row for each segment.
##' @author Toby Dylan Hocking
##' @param x numeric vector of data to fit a Gaussian mean model.
##' @param labels data frame with at least three columns: start, end,
##'   changes. start/end should be indices of x, from 1 to
##'   length(x). changes should be either 0 or 1. The prediced
##'   changepoints are guaranteed to be consistent with these labels.
##' @param penalty_unlabeled non-negative penalty constant (larger for
##'   fewer changes, smaller for more changes). penalty=0 means a
##'   change in every unlabeled region, penalty=Inf means no changes
##'   in unlabeled regions.
##' @param n_updates how many dynamic programming updates to compute?
##'   Must be at least 1 and at most length(x).
##' @param penalty_labeled non-negative penalty constant to use for
##'   changes in positive labels. 
##' @example inst/examples/LOPART.R
LOPART <- function
(x, labels, penalty_unlabeled,
  n_updates=length(x),
  penalty_labeled=penalty_unlabeled
){
  last_change <- cost_optimal <- . <- start <- end <- changes <- NULL
  ## above to avoid CRAN NOTE.
  Inf.penalty <- penalty_unlabeled == Inf
  if(Inf.penalty){
    ## special case for efficient computation of model with infinite
    ## penalty.
    pos.dt <- data.table(labels)[changes==1][order(start)]
    neg.dt <- pos.dt[, data.table(
      start=c(1L, end),
      end=c(start, length(x)),
      changes=0L)]
    labels <- rbind(
      pos.dt[, .(start, end, changes)], neg.dt
    )[start < end][order(start)]
  }
  out_df <- LOPART_interface(
    x,
    labels$start-1L,
    labels$end-1L,
    labels$changes,
    n_updates,
    if(Inf.penalty)0 else penalty_unlabeled,
    if(Inf.penalty)0 else penalty_labeled)
  out_dt <- data.table(out_df)
  change.vec <- out_dt[0 <= last_change, last_change+1L]
  changes_total <- length(change.vec)
  changes_labeled <- sum(labels$changes==1)
  changes_unlabeled <- changes_total-changes_labeled
  multiply <- function(changes, pen){
    if(changes==0)0 else changes*pen
  }
  complexity_labeled <- multiply(changes_labeled, penalty_labeled)
  complexity_unlabeled <- multiply(changes_unlabeled, penalty_unlabeled)
  complexity_total <- complexity_labeled+complexity_unlabeled
  if(Inf.penalty){
    total_loss <- out_dt[.N, cost_optimal]
    penalized_cost <- total_loss+complexity_total
  }else{
    penalized_cost <- out_dt[.N, cost_optimal]
    total_loss <- penalized_cost-complexity_total
  }
  list(
    loss=data.table(
      changes_total,
      changes_labeled,
      changes_unlabeled,
      penalty_labeled,
      penalty_unlabeled,
      penalized_cost,
      total_loss),
    cost=out_dt,
    changes=data.table(
      change=change.vec+0.5),
    segments=out_dt[last_change != -2, .(
      start=c(1L, change.vec+1L),
      end=c(change.vec, nrow(out_dt)),
      mean)])
}
