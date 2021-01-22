

weight = as.integer(runif(10) * 100)
value = as.integer(runif(10) * 1000)
cap = as.integer(mean(weight))
tmp = order(-value / weight)
weight = weight[tmp]
value = value[tmp]


value <- c(15, 100, 90, 60, 40, 15, 10,  1)
weight =  c( 2,  20, 20, 30, 40, 30, 60, 10)
cap = 102

sink("sink.txt")
tmp = FLSSS:::bkpOrdered(weight, value, cap) + 1L
sink()


adagio::knapsack(weight, value, cap)



weight = c(1, 382745,	799601,	909247,	729069,	467902,	44328,	34610,	698150,	823460,	903959,	853665,	551830,	610856,	670702,	488960,	951111,	323046,	446298,	931161,	31385,	496951,	264724,	224916,	169684)
value = c(1000, 825594,	1677009,	1676628,	1523970,	943972,	97426,	69666,	1296457,	1679693,	1902996,	1844992,	1049289,	1252836,	1319836,	953277,	2067538,	675367,	853655,	1826027,	65731,	901489,	577243,	466257,	369261)
value = value[order(weight)]
cap = c(1e10, 6404180 * seq(1, by = -0.01, len = 20), 1e6)


sink("sink.txt")
# tmp = FLSSS:::extra01knapsackBaB(weight, value, cap, itemNcaps = rep(11L, length(cap)))
# tmp = FLSSS:::extra01knapsackBaB(weight, value, cap)
tmp = FLSSS:::extra01knapsackDP(weight, value, cap, maxCore = 7)
sink()




weight = as.integer(c(100000, 382745,	799601,	909247,	729069,	467902,	44328,	34610,	698150,	823460,	903959,	853665,	551830,	610856,	670702,	488960,	951111,	323046,	446298,	931161,	31385,	496951,	264724,	224916,	169684) / 1000)
value = c(1000000, 825594,	1677009,	1676628,	1523970,	943972,	97426,	69666,	1296457,	1679693,	1902996,	1844992,	1049289,	1252836,	1319836,	953277,	2067538,	675367,	853655,	1826027,	65731,	901489,	577243,	466257,	369261)
caps = as.integer(c(sum(weight) + 1L, 6404180 / 1000 * seq(1, by = -0.01, len = 20), 100))
sink("sink.txt")
tmp = FLSSS:::auxKnapsack01dp(weight, value, caps, maxCore = 7, tlimit = 0.01)
sink()




tmp2 = FLSSS:::auxKnapsack01bb(weight, value, caps, maxCore = 7)
tmp2$selection = lapply(tmp2$selection, function(x) sort(x))
unlist(lapply(tmp$selection, function(x) sum(value[x])))



#
# https://drive.google.com/file/d/13ioVxTUn9ptdnbkA16FWNgPt9v0ms69V/view?usp=sharing
# https://drive.google.com/file/d/18mjaLg9bb3JZegrXVDQsdFrEOja-_h_F/view?usp=sharing
# https://drive.google.com/file/d/1rKdRx227ctu-EMKpBC5hV6J8QTbzxLeH/view?usp=sharing
# https://drive.google.com/file/d/1lvSwHfRF497cIevrRE-sG7kbJK3ipPdo/view?usp=sharing
# https://drive.google.com/file/d/1FkdfHaVt1m7piR68MTAvCOcG3bKJdpV4/view?usp=sharing


load(url("https://drive.google.com/file/d/1XTzYWUXSvTL_i1c6TaoqJvAkYsHvVdha/view?usp=sharing"))












weight = as.integer(c(382745,	799601,	909247,	729069,	467902,	44328,	34610,	698150,	823460,	903959,	853665,	551830,	610856,	670702,	488960,	951111,	323046,	446298,	931161,	31385,	496951,	264724,	224916,	169684) / 1000)
value = c(825594,	1677009,	1676628,	1523970,	943972,	97426,	69666,	1296457,	1679693,	1902996,	1844992,	1049289,	1252836,	1319836,	953277,	2067538,	675367,	853655,	1826027,	65731,	901489,	577243,	466257,	369261)
value = sample(value, length(value))
cap = as.integer(6404180 / 1000)
tmp = FLSSS:::auxBiknapsackDP(weight, value, cap)
tmp2 = FLSSS:::auxBiknapsackBaB(weight, value, cap)
tmp$maxValue
abs(tmp$maxValue / tmp2$maxVal - 1)
sum(abs(sort(tmp2$selection) - tmp$selection))


# The equivalent minimization problem
tmp3 = FLSSS:::extra01knapsackBaB(-weight, -value, cap - sum(weight))





sink("sink.txt")
tmp = FLSSS:::extra01knapsackBaB(weight, value, c(6404180, 6404180 * 0.5, 6404180 * 100000, 100))
sink()




{
totalP <- function(x, profit)
{
  sum(apply(cbind(1L : length(x), x), 1, function(x) profit[x[1], x[2]]))
}
totalW <- function(x, weight)
{
  S = numeric(ncol(weight))
  for(i in 1L : length(x))
  {
    S[x[i]] = S[x[i]] + weight[i, x[i]]
  }
  S
}
}
weight = matrix(c(21, 13 , 9 , 5 , 7 ,15 , 5,24, 20  ,8 ,18, 25 , 6 , 6 , 9, 6, 16, 16, 18 ,24, 11 ,11, 16, 18), ncol = 3)
profit = matrix(c(27, 12 ,12, 16 ,24 ,31 ,41, 13, 14,  5 ,37,  9 ,36, 25 , 1 ,34, 34 ,34, 20,  9, 19, 19 , 3 ,34), ncol = 3)
budget = c(26, 25, 34)
sink("sink.txt")
tmp = FLSSS:::extraGAP(weight, profit, budget, maxCore = 7)
sink()




weight = t(matrix(c(21, 13 , 9 , 5 , 7 ,15 , 5,24, 20  ,8 ,18, 25 , 6 , 6 , 9, 6, 16, 16, 18 ,24, 11 ,11, 16, 18), ncol = 3))
profit = t(matrix(c(27, 12 ,12, 16 ,24 ,31 ,41, 13, 14,  5 ,37,  9 ,36, 25 , 1 ,34, 34 ,34, 20,  9, 19, 19 , 3 ,34), ncol = 3))
budget = c(26, 25, 34)
sink("sink.txt")
tmp = FLSSS:::auxGAPbab(weight, profit, budget, maxCore = 1)
sink()


sink("sink.txt")
tmp2 = FLSSS:::auxGAPbabDp(weight, profit, budget, maxCore = 1)
sink()























sol = c(3, 3, 1, 1, 2, 2, 1, 2)
totalP(sol, profit)
totalW(sol, weight)
budget




value = matrix(c(14, 38, 1, 26, 14, 49, 1, 20, 46, 49, 1, 45, 16, 10, 1), ncol = 3)
value = 50 - value
weight = matrix(c(12, 19, 1000, 11, 18, 6, 1000, 11, 15, 18, 1000, 10, 14, 22, 1000), ncol = 3)
budget = c(28, 28, 28)
sink("sink.txt")
tmp = FLSSS:::extraGAP(weight, value, budget, maxCore = 1)
sink()


tmp = FLSSS:::extra01knapsackBaB(c(12, 19, 18), c(35, 7, 35), sum(c(12, 19, 18)) - 21)














agents = 4L
tasks = 10L
costs = t(as.data.frame(lapply(1L : agents, function(x) pmax(1L, as.integer(runif(tasks) * 100)))))
budgets = as.integer(apply(costs, 1, function(x) runif(1, min(x), sum(x) / 2)))
profits = t(as.data.frame(lapply(1L : agents, function(x)
  abs(rnorm(tasks) + runif(1, 0, 4)) * 10000)))
dimnames(costs) = NULL
dimnames(profits) = NULL
names(budgets) = NULL
rst = FLSSS::GAP(maxCore = 7L, agentsCosts = costs, agentsProfits = profits,
                 agentsBudgets = budgets, heuristic = FALSE, tlimit = 60,
                 threadLoad = 8L, verbose = TRUE)
# Function also saves the assignment costs and profits


sink("sink.txt")
tmp = FLSSS:::auxGAPbbDp(costs, profits, budgets, maxCore = 7, greedyBranching = TRUE)
tmp2 = FLSSS:::auxGAPbb(costs, profits, budgets, maxCore = 7, greedyBranching = TRUE)
sink()
rst$assignedAgents$agent; tmp$assignment
rst$assignmentProfit; tmp$totalProfit




w = c(352.719, 608.055, 3.68912)
p = c(20580.8, 31248.7, 22221.2)
FLSSS:::extra01knapsackBaB(w, p, sum(w) - 100.413)




agents = 4L
tasks = 8L
while(T)
{
  costs = t(as.data.frame(lapply(1L : agents, function(x) pmax(1L, as.integer(runif(tasks) * 100)))))
  budgets = as.integer(apply(costs, 1, function(x) runif(1, min(x), sum(x) / 1.1)))
  profits = t(as.data.frame(lapply(1L : agents, function(x)
    abs(rnorm(tasks) + runif(1, 0, 4)) * 10000)))
  dimnames(costs) = NULL
  dimnames(profits) = NULL
  names(budgets) = NULL
  save.image()
  sol = FLSSS:::auxGAPbb(costs, profits, budgets, maxCore = 7, tlimit = 60, ub = "MT", greedyBranching = TRUE, optim = "max")
}















