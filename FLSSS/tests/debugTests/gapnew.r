

set.seed(42)
weight = runif(100, min = 1e3, max = 1e6)
value = weight ^ 0.5 * 100 # Higher correlation between item weights and values
# typically implies a harder knapsack problem.
caps = runif(10, min(weight), sum(weight))
rst = FLSSS::auxKnapsack01bb(weight, value, caps, maxCore = 2, tlimit = 2000, ub = "HS")
str(rst)



























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
# 1 1 3 2 0 4 0 1 0 3 3 3 0 4 2


Nagent = 5L; Ntask = 15L
sink("sink.txt")
rst = FLSSS:::auxGAPga(cost, profit, budget, trials = 1, randomSeed = NULL, populationSize = 100, generations = 10000, optim = "max", maxCore = 1)
sink()
str(rst)




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
  cost, profit, budget, trials = 2, randomSeed = 42, populationSize = 100,
  generations = 10000, optim = "max", maxCore = 2)







# Examination
tmp = apply(rst$populationInfo$allGenes, 2, function(x)
{
  s = 0
  for(i in 1L : length(x)) s = s + profit[x[i], i]
  s
})
range(abs(tmp - rst$populationInfo$allProfitOrLoss))


tmp2 = apply(rst$populationInfo$allGenes, 2, function(x)
{
  s = 0
  b = numeric(Nagent)
  for(i in 1L : length(x)) b[x[i]] = b[x[i]] + cost[x[i], i]
  sum(pmax(0, b - budget))
})
range(abs(tmp2 - rst$populationInfo$allBudgetExceedance))


FLSSS::auxGAPbb(cost, profit, budget, optim = "min")









load("C:/Users/i56087/Downloads/gapInstances.Rdata")
cost = gapC[[3]]$cost
loss = gapC[[3]]$loss
budget = gapC[[3]]$budget
system.time({rst = FLSSS:::auxGAPga(cost, loss, budget, trials = 10, randomSeed = 42, populationSize = 100, generations = 500000, optim = "min", maxCore = 7)})


# Examination
tmp = apply(rst$populationInfo$allGenes, 2, function(x)
{
  s = 0
  for(i in 1L : length(x)) s = s + loss[x[i], i]
  s
})


tmp = apply(rst$populationInfo$allGenes, 2, function(x)
{
  s = 0
  for(i in 1L : length(x)) s = s + cost[x[i], i]
  s -
})























