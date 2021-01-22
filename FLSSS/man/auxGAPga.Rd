\name{auxGAPga}
\alias{auxGAPga}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multithreaded generalized assignment problem solver via genetic algorithm
}
\description{
A genetic algorithm with local heuristics for GAP.
}
\usage{
auxGAPga(
  cost,
  profitOrLoss,
  budget,
  trials,
  populationSize,
  generations,
  randomSeed = NULL,
  maxCore = 7,
  optim = "max"
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cost}{
A numeric matrix. Dimensionality = N(agents) \code{x} N(tasks).
}
  \item{profitOrLoss}{
A numeric matrix of the same dimensionality of \code{cost}. Profit for maximum GAP. Loss for minimum GAP.
}
  \item{budget}{
A numeric vector. Size = N(agents).
}
  \item{trials}{
An integer. Number of trials, aka the number of population sets.
}
  \item{populationSize}{
An integer. Size of each population.
}
  \item{generations}{
An integer. As reproduction iterates, if there have been \code{generations} many children produced and accepted in population but no update on the current optimum, function quits.
}
  \item{randomSeed}{
An integer or \code{NULL}. \code{randomSeed} seeds the random number generator in R, generates \code{trials} many integers to seed the \code{mt19937_64} (Mersenne Twister) engine for each trial.
}
  \item{maxCore}{
Maximal threads to invoke. No greater than the number of logical CPUs on machine. The algorithm multithreads over trials.
}
  \item{optim}{
A string. \code{optim = "max"} (\code{"min"}) solves the maximum (minimum) GAP.
}
}
\details{
This algorithm is based on a foundational paper by Chu and Beasley (1997) and is carefully engineered towards speed. Besides the standard cross-over and mutation operations, the algorithm applies two local heuristics for educating the new borns. The first is to randomly pick a task from each overloaded agent and reassign the task to the next budget-sufficient agent --- if there is any. The second is to raise the total profit by reassigning another agent for each task --- if the reassignment would not result in overload. The algorithm outperforms most peer metaheuristics such as variants of simulated annealing and tabu search (Osman), and is highly effective for large and hard instances.
}
\value{
A list of 4:
\item{totalProfitOrLoss }{Total profit or loss generated from the assignment. Negative infinity if no solution found.}

\item{agentCost }{A numeric vector of total costs for each agent. Empty if no solution found.}

\item{assignment }{An integer vector. \code{assignment[i]} indexes the agent assigned to the \code{i}th task. Empty if no solution found.}

\item{populationInfo }{A list of 3:}

\item{allGenes }{An N(task) \code{x} (\code{populationSize x trials}) integer matrix recording genes in all population sets upon completion. Each column represents a gene, namely a tentative assignment.}

\item{allBudgetExceedance }{A numeric vector of the size of \code{populationSize x trials}. \code{allBudgetExceedance[i]} equals the total budget exceedance of tentative assignment \code{allGenes[, i]}.}

\item{allProfitOrLoss }{A numeric vector of the size of \code{allBudgetExceedance}. \code{allProfitOrLoss[i]} equals the total profit or loss of tentative assignment \code{allGenes[, i]}.}
}

\note{
The C++ implementation is fully independent and borrows no code from any commercial or open source.
}


\examples{
# =============================================================================
# A trivial instance
# =============================================================================
profit = c(17,21,22,18,24,15,20,18,19,18,16,22,24,24,16,23,16,21,16,17,16,19,
           25,18,21,17,15,25,17,24,16,20,16,25,24,16,17,19,19,18,20,16,17,21,
           24,19,19,22,22,20,16,19,17,21,19,25,23,25,25,25,18,19,15,15,21,25,
           16,16,23,15,22,17,19,22,24)
profit = t(matrix(profit, ncol = 5))
cost = c(8,15,14,23,8,16,8,25,9,17,25,15,10,8,24,15,7,23,22,11,11,12,10,17,16,
         7,16,10,18,22,21,20,6,22,24,10,24,9,21,14,11,14,11,19,16,20,11,8,14,
         9,5,6,19,19,7,6,6,13,9,18,8,13,13,13,10,20,25,16,16,17,10,10,5,12,23)
cost = t(matrix(cost, ncol = 5))
budget = c(36, 34, 38, 27, 33)


Nagent = 5L; Ntask = 15L
rst = FLSSS::auxGAPga(
  cost, profit, budget, trials = 2, populationSize = 100, generations = 10000,
  randomSeed = 42, maxCore = 2, optim = "max")


# =============================================================================
# A relatively hard instance.
# =============================================================================
# Download gapInstances.Rdata from
# https://github.com/WhateverLiu/gapInstances. Load it in R.
\donttest{
cost = gapC[[3]]$cost
loss = gapC[[3]]$loss
budget = gapC[[3]]$budget
# Intel CPU i7-4770 3.4GHz, g++ '-Ofast', 64-bit Windows 7.
system.time({rst = FLSSS::auxGAPga(cost, loss, budget, trials = 7,
             randomSeed = 42, populationSize = 100, generations = 500000,
             optim = "min", maxCore = 7)})
rst$totalProfitOrLoss # 1416
#  user  system elapsed
# 69.24    0.17   11.61
# The known optimum equals 1402 as the total loss.
}
}




















