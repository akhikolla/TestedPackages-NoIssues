## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(fig.width = 6, fig.height = 5, fig.align = "center")
options(rmarkdown.html_vignette.check_title = FALSE)
load("vignettedata/vignettedata.Rdata")

## ----CRAN, eval = FALSE-------------------------------------------------------
#  install.packages("cutpointr")

## -----------------------------------------------------------------------------
library(cutpointr)
data(suicide)
head(suicide)
cp <- cutpointr(suicide, dsi, suicide, 
                method = maximize_metric, metric = sum_sens_spec)

## -----------------------------------------------------------------------------
summary(cp)

## -----------------------------------------------------------------------------
plot(cp)

## -----------------------------------------------------------------------------
opt_cut <- cutpointr(suicide, dsi, suicide, direction = ">=", pos_class = "yes",
                     neg_class = "no", method = maximize_metric, metric = youden)

## -----------------------------------------------------------------------------
plot_metric(opt_cut)

## -----------------------------------------------------------------------------
predict(opt_cut, newdata = data.frame(dsi = 0:5))

## ----separate subgroups and bootstrapping, eval = FALSE-----------------------
#  set.seed(12)
#  opt_cut_b <- cutpointr(suicide, dsi, suicide, boot_runs = 1000)

## -----------------------------------------------------------------------------
opt_cut_b

## -----------------------------------------------------------------------------
opt_cut$boot

## -----------------------------------------------------------------------------
summary(opt_cut)
plot(opt_cut)

## ---- eval = FALSE------------------------------------------------------------
#  library(doParallel)
#  cl <- makeCluster(2) # 2 cores
#  registerDoParallel(cl)
#  registerDoRNG(12) # Reproducible parallel loops using doRNG
#  opt_cut <- cutpointr(suicide, dsi, suicide, gender, pos_class = "yes",
#                       direction = ">=", boot_runs = 1000, allowParallel = TRUE)
#  stopCluster(cl)
#  opt_cut

## ---- cache=TRUE--------------------------------------------------------------
set.seed(100)
cutpointr(suicide, dsi, suicide, gender, 
          method = maximize_boot_metric,
          boot_cut = 200, summary_func = mean,
          metric = accuracy, silent = TRUE)

## -----------------------------------------------------------------------------
opt_cut <- cutpointr(suicide, dsi, suicide, gender, method = minimize_metric,
                     metric = misclassification_cost, cost_fp = 1, cost_fn = 10)

## -----------------------------------------------------------------------------
plot_metric(opt_cut)

## ---- message = FALSE---------------------------------------------------------
opt_cut <- cutpointr(suicide, dsi, suicide, gender, 
                     method = minimize_loess_metric,
                     criterion = "aicc", family = "symmetric", 
                     degree = 2, user.span = 0.7,
                     metric = misclassification_cost, cost_fp = 1, cost_fn = 10)

## -----------------------------------------------------------------------------
plot_metric(opt_cut)

## -----------------------------------------------------------------------------
library(ggplot2)
exdat <- iris
exdat <- exdat[exdat$Species != "setosa", ]
opt_cut <- cutpointr(exdat, Petal.Length, Species,
                     method = minimize_gam_metric,
                     formula = m ~ s(x.sorted, bs = "cr"),
                     metric = abs_d_sens_spec)
plot_metric(opt_cut)

## -----------------------------------------------------------------------------
opt_cut <- cutpointr(suicide, dsi, suicide, gender, 
                     method = minimize_spline_metric, spar = 0.4,
                     metric = misclassification_cost, cost_fp = 1, cost_fn = 10)
plot_metric(opt_cut)

## -----------------------------------------------------------------------------
cutpointr(suicide, dsi, suicide, gender, method = oc_youden_normal)

## -----------------------------------------------------------------------------
cutpointr(suicide, dsi, suicide, gender, method = oc_youden_kernel)

## ---- fig.width=4, fig.height=3-----------------------------------------------
roc_curve <- roc(data = suicide, x = dsi, class = suicide,
    pos_class = "yes", neg_class = "no", direction = ">=")
auc(roc_curve)
head(roc_curve)
plot_roc(roc_curve)

## -----------------------------------------------------------------------------
dat <- data.frame(outcome = c("neg", "neg", "neg", "pos", "pos", "pos", "pos"),
                  pred    = c(1, 2, 3, 8, 11, 11, 12))

## -----------------------------------------------------------------------------
opt_cut <- cutpointr(dat, x = pred, class = outcome, use_midpoints = TRUE)
plot_x(opt_cut)

## ---- echo = FALSE------------------------------------------------------------
plotdat_nomidpoints <- structure(list(sim_nr = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 
4L, 4L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 
6L, 6L, 6L, 6L, 6L, 6L, 6L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 8L, 
8L, 8L, 8L, 8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 10L, 
10L, 10L, 10L, 10L, 10L, 10L, 10L, 11L, 11L, 11L, 11L, 11L, 11L, 
11L, 11L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 13L, 13L, 13L, 
13L, 13L, 13L, 13L, 13L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 
15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 16L, 16L, 16L, 16L, 16L, 
16L, 16L, 16L, 17L, 17L, 17L, 17L, 17L, 17L, 17L, 17L, 18L, 18L, 
18L, 18L, 18L, 18L, 18L, 18L, 19L, 19L, 19L, 19L, 19L, 19L, 19L, 
19L, 20L, 20L, 20L, 20L, 20L, 20L, 20L, 20L, 21L, 21L, 21L, 21L, 
21L, 21L, 21L, 21L, 22L, 22L, 22L, 22L, 22L, 22L, 22L, 22L, 23L, 
23L, 23L, 23L, 23L, 23L, 23L, 23L, 24L, 24L, 24L, 24L, 24L, 24L, 
24L, 24L, 25L, 25L, 25L, 25L, 25L, 25L, 25L, 25L, 26L, 26L, 26L, 
26L, 26L, 26L, 26L, 26L, 27L, 27L, 27L, 27L, 27L, 27L, 27L, 27L, 
28L, 28L, 28L, 28L, 28L, 28L, 28L, 28L, 29L, 29L, 29L, 29L, 29L, 
29L, 29L, 29L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 31L, 31L, 
31L, 31L, 31L, 31L, 31L, 31L, 32L, 32L, 32L, 32L, 32L, 32L, 32L, 
32L, 33L, 33L, 33L, 33L, 33L, 33L, 33L, 33L, 34L, 34L, 34L, 34L, 
34L, 34L, 34L, 34L, 35L, 35L, 35L, 35L, 35L, 35L, 35L, 35L, 36L, 
36L, 36L, 36L, 36L, 36L, 36L, 36L), method = structure(c(1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 
2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L), .Label = c("emp", 
"normal", "loess", "boot", "spline", "spline_20", "kernel", "gam"
), class = "factor"), n = c(30, 30, 30, 30, 30, 30, 30, 30, 50, 
50, 50, 50, 50, 50, 50, 50, 75, 75, 75, 75, 75, 75, 75, 75, 100, 
100, 100, 100, 100, 100, 100, 100, 150, 150, 150, 150, 150, 150, 
150, 150, 250, 250, 250, 250, 250, 250, 250, 250, 500, 500, 500, 
500, 500, 500, 500, 500, 750, 750, 750, 750, 750, 750, 750, 750, 
1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 30, 30, 30, 30, 
30, 30, 30, 30, 50, 50, 50, 50, 50, 50, 50, 50, 75, 75, 75, 75, 
75, 75, 75, 75, 100, 100, 100, 100, 100, 100, 100, 100, 150, 
150, 150, 150, 150, 150, 150, 150, 250, 250, 250, 250, 250, 250, 
250, 250, 500, 500, 500, 500, 500, 500, 500, 500, 750, 750, 750, 
750, 750, 750, 750, 750, 1000, 1000, 1000, 1000, 1000, 1000, 
1000, 1000, 30, 30, 30, 30, 30, 30, 30, 30, 50, 50, 50, 50, 50, 
50, 50, 50, 75, 75, 75, 75, 75, 75, 75, 75, 100, 100, 100, 100, 
100, 100, 100, 100, 150, 150, 150, 150, 150, 150, 150, 150, 250, 
250, 250, 250, 250, 250, 250, 250, 500, 500, 500, 500, 500, 500, 
500, 500, 750, 750, 750, 750, 750, 750, 750, 750, 1000, 1000, 
1000, 1000, 1000, 1000, 1000, 1000, 30, 30, 30, 30, 30, 30, 30, 
30, 50, 50, 50, 50, 50, 50, 50, 50, 75, 75, 75, 75, 75, 75, 75, 
75, 100, 100, 100, 100, 100, 100, 100, 100, 150, 150, 150, 150, 
150, 150, 150, 150, 250, 250, 250, 250, 250, 250, 250, 250, 500, 
500, 500, 500, 500, 500, 500, 500, 750, 750, 750, 750, 750, 750, 
750, 750, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000), mean_err = c(0.532157164015659, 
0.0344907054484091, 1.09430750651166, 0.847845162156675, 1.72337372126503, 
0.893756658507988, 0.0430309247027736, 0.785821459035346, 0.368063404388512, 
0.0256197760404459, 0.54480529648463, 0.54385597929651, 0.657325657699579, 
0.578611116865437, 0.0400491342691897, 0.515688005217413, 0.256713912589642, 
0.0444582875885996, 0.326975493112402, 0.371128780921122, 0.473515115741104, 
0.389519558405289, 0.105044360789378, 0.301924717299333, 0.207750921776918, 
-0.00318128936770314, 0.215170156089776, 0.27218780048926, 0.260519564021842, 
0.236792923882582, 0.0209319074923902, 0.232957055204834, 0.0726605917614469, 
-0.00282823355849125, 0.0753216783313991, 0.147121931849656, 
0.0986417724955371, 0.10048009778446, -0.0117861260923649, 0.0650845904350442, 
0.0985144485083747, 0.00601003227249322, 0.107439979908118, 0.120777421732797, 
0.098470489820427, 0.0940946984227826, 0.0340166854141625, 0.107851118082414, 
0.0249685210781582, -0.00275219600614378, 0.0258069390207201, 
0.0303381972274654, -0.000994602151869198, 0.00196854833683764, 
-0.0172319489562159, 0.0230957871932473, 0.00787486424680835, 
-0.018438041997315, -0.000808033567394628, -0.00151904153864496, 
-0.0258118523805697, -0.020984156892953, -0.0411584927473141, 
-0.0462075435919094, 0.0481149217843661, -0.0115241085997692, 
0.0278045708419731, 0.0358588316426625, 0.0424473909450939, 0.0379773233328197, 
0.0298772985879321, 0.0494939492036379, 0.561286668337778, 0.0210874502384648, 
0.607711822769155, 0.944733906256477, 1.32069801051061, 0.623604782428556, 
0.0138075109769806, 0.640859854412358, 0.284873604303057, -0.0170357985365701, 
0.273426417633118, 0.524432737895336, 0.355003110979807, 0.312837607951434, 
-0.0316296929553873, 0.270109834098986, 0.174110335581819, 0.0253101719615279, 
0.199956222742702, 0.375485416120771, 0.278956745944806, 0.244245525945888, 
0.0325126314233263, 0.253352659868514, 0.154840760004461, 0.00231589709639472, 
0.154412480165179, 0.264847742842386, 0.196572744185608, 0.182934783774992, 
0.0207139021497755, 0.17351041376412, 0.14190910348156, -0.000766834484010096, 
0.159975205214477, 0.191222128926019, 0.119768252112669, 0.12033372914036, 
-0.00429047209392067, 0.120982527821078, 0.0756304869484406, 
-0.00890884219048113, 0.0727782693168392, 0.118690444738942, 
0.0814898789647033, 0.0799724348001957, 0.0182926240912726, 0.0887155007804252, 
0.00799604720502299, 0.000599148388616836, -0.00567769035990384, 
0.0358412514670032, 0.0308474979074875, 0.0341668723768997, -0.000180318451026095, 
0.0180733341290925, 0.00456876236626807, 0.00150574966485876, 
0.011152095953916, 0.0176039119729626, 0.00608274255434991, 0.0146257828313115, 
-0.0108877417404102, 0.00341198000323035, -0.00198459880370283, 
0.0026551895445694, 0.00199040664534129, 0.0150165794544221, 
0.00646287144368147, 0.00999205240904708, -0.00850278571195971, 
0.000833666619266177, 0.714730067273087, -0.00916546079360956, 
0.662799490366986, 1.18552468844156, 1.25901933062308, 0.672701515532179, 
0.0311066140197676, 0.699068058809396, 0.451908043813962, 0.0131716226592205, 
0.429551369887738, 0.697928133235757, 0.445367768988423, 0.408463448982185, 
0.0318707241721211, 0.406284953982951, 0.257736307754364, -0.00525924423719458, 
0.236977000055322, 0.429144726596141, 0.291381752107184, 0.267557606428613, 
0.0103657879176852, 0.254728590646094, 0.187783398487578, -0.00216064381479362, 
0.209025103860707, 0.318293592390017, 0.216751610346408, 0.195630579126633, 
-0.00355723971644246, 0.174111826408428, 0.151010324964235, 0.0152409223899092, 
0.159002511320467, 0.214643583389694, 0.136211731513269, 0.138948149207635, 
0.00736196817594524, 0.115637867729083, 0.0491348055596302, -0.00133957946235943, 
0.0507437758212659, 0.103956325245849, 0.0641182216839426, 0.0721933081297794, 
-0.0124376134651938, 0.0632317888879588, 0.0322195712438111, 
0.00170122889182022, 0.0287526766624194, 0.0589662164030242, 
0.0348535721709848, 0.039527944642463, -0.00617539706415593, 
0.0274246010641889, 0.0325877909680824, 0.00530528253248245, 
0.0221776555499961, 0.0389702052631117, 0.0221602091288215, 0.0254478639695596, 
-0.0016189234058987, 0.0197144417326668, -0.00632485262604172, 
-0.00364979854195596, -0.00276076468388984, 0.0126267527301874, 
0.0123592498266038, 0.0154921632247644, -0.00591512196680815, 
0.0098685016547149, 1.19276916750486, -0.0296831640401583, 0.99406393593888, 
1.95758669445116, 1.45842790978446, 0.916899913902239, -0.0240222410217233, 
1.00771193034927, 0.748151091428865, 0.001671855025917, 0.665180535306263, 
1.1777049634557, 0.578603609273264, 0.546625362141714, 0.0292152981607387, 
0.615230814912951, 0.417886753756131, 0.00324593885807739, 0.406076310942717, 
0.732191741449251, 0.352684769616612, 0.326901376027897, -0.000759357576337989, 
0.350075431324921, 0.310927617707656, -0.0107255472998434, 0.28102101085112, 
0.514683023356017, 0.24913510139508, 0.235155452507568, -0.0220885572014814, 
0.243370611433649, 0.209652330609093, -0.00502865663759991, 0.2172246261346, 
0.356540958804122, 0.172121720418057, 0.17487914828986, 0.00365942442127361, 
0.176594681455494, 0.126956927057327, -0.00270525933073803, 0.120116234221594, 
0.210827536708082, 0.101520409193932, 0.101379097920023, 0.00238043252144371, 
0.113027315928011, 0.0598624378953727, -0.00538838415690431, 
0.0568400730102315, 0.0978115258288965, 0.0454207684906316, 0.0473140143579152, 
-0.00165813015281622, 0.0521772135812508, 0.0530224090961669, 
-0.000416993405198653, 0.0353236911458531, 0.0605493601241619, 
0.0316204159297213, 0.0344789374555544, -0.00446984887315054, 
0.0328807595695966, 0.0396438546423947, -0.00331466369719113, 
0.0379029847219126, 0.0572435100638761, 0.0253269328104989, 0.0235663211070417, 
0.00220241478536399, 0.0307132312422208), youden = c(0.2, 0.2, 
0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.6, 
0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 
0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 
0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 
0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 
0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 
0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 
0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 
0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 
0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 
0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 
0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8
)), row.names = c(NA, -288L), group_sizes = c(8L, 8L, 8L, 8L, 
8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 
8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L
), biggest_group_size = 8L, class = c("grouped_df", "tbl_df", 
"tbl", "data.frame"), groups = structure(list(sim_nr = 1:36, 
    .rows = list(1:8, 9:16, 17:24, 25:32, 33:40, 41:48, 49:56, 
        57:64, 65:72, 73:80, 81:88, 89:96, 97:104, 105:112, 113:120, 
        121:128, 129:136, 137:144, 145:152, 153:160, 161:168, 
        169:176, 177:184, 185:192, 193:200, 201:208, 209:216, 
        217:224, 225:232, 233:240, 241:248, 249:256, 257:264, 
        265:272, 273:280, 281:288)), row.names = c(NA, -36L), class = c("tbl_df", 
"tbl", "data.frame"), .drop = TRUE))

## ---- echo = FALSE------------------------------------------------------------
library(dplyr)
ggplot(plotdat_nomidpoints %>% filter(!(method %in% c("spline_20"))), 
       aes(x = n, y = mean_err, color = method, shape = method)) + 
    geom_line() + geom_point() +
    facet_wrap(~ youden, scales = "fixed") +
    scale_shape_manual(values = 1:nlevels(plotdat_nomidpoints$method)) +
    scale_x_log10(breaks = c(30, 50, 75, 100, 150, 250, 500, 750, 1000)) +
    ggtitle("Bias of all methods when use_midpoints = FALSE",
            "normally distributed data, 10000 repetitions of simulation")

## -----------------------------------------------------------------------------
opt_cut <- cutpointr(suicide, dsi, suicide, metric = sum_sens_spec, 
                     tol_metric = 0.05, break_ties = c)
library(tidyr)
opt_cut %>% 
    select(optimal_cutpoint, sum_sens_spec) %>% 
    unnest(cols = c(optimal_cutpoint, sum_sens_spec))

## ---- eval = FALSE------------------------------------------------------------
#  set.seed(100)
#  opt_cut_manual <- cutpointr(suicide, dsi, suicide, method = oc_manual,
#                         cutpoint = mean(suicide$dsi), boot_runs = 1000)
#  set.seed(100)
#  opt_cut_mean <- cutpointr(suicide, dsi, suicide, method = oc_mean, boot_runs = 1000)

## ---- eval = FALSE------------------------------------------------------------
#  myvar <- "dsi"
#  cutpointr(suicide, !!myvar, suicide)

## -----------------------------------------------------------------------------
mcp <- multi_cutpointr(suicide, class = suicide, pos_class = "yes", 
                use_midpoints = TRUE, silent = TRUE) 
summary(mcp)

## ---- eval = FALSE, message = FALSE-------------------------------------------
#  set.seed(123)
#  opt_cut_b_g <- cutpointr(suicide, dsi, suicide, gender, boot_runs = 1000)

## ---- message = FALSE---------------------------------------------------------
# Using dplyr and tidyr
library(tidyr)
opt_cut_b_g %>% 
  group_by(subgroup) %>% 
  select(subgroup, boot) %>%
  unnest(cols = boot) %>% 
  summarise(sd_oc_boot = sd(optimal_cutpoint),
            m_oc_boot  = mean(optimal_cutpoint),
            m_acc_oob  = mean(acc_oob))

## -----------------------------------------------------------------------------
cutpointr(suicide, dsi, suicide, gender, metric = youden, silent = TRUE) %>% 
    add_metric(list(ppv, npv)) %>% 
    select(subgroup, optimal_cutpoint, youden, ppv, npv)

## -----------------------------------------------------------------------------
roc(data = suicide, x = dsi, class = suicide, pos_class = "yes",
    neg_class = "no", direction = ">=") %>% 
  add_metric(list(cohens_kappa, F1_score)) %>% 
  select(x.sorted, tp, fp, tn, fn, cohens_kappa, F1_score) %>% 
  head()

## ---- eval = FALSE------------------------------------------------------------
#  mean_cut <- function(data, x, ...) {
#      oc <- mean(data[[x]])
#      return(data.frame(optimal_cutpoint = oc))
#  }

## -----------------------------------------------------------------------------
misclassification_cost

## ---- fig.width=4, fig.height=3-----------------------------------------------
plot_cut_boot(opt_cut_b_g)
plot_metric(opt_cut_b_g, conf_lvl = 0.9)
plot_metric_boot(opt_cut_b_g)
plot_precision_recall(opt_cut_b_g)
plot_sensitivity_specificity(opt_cut_b_g)
plot_roc(opt_cut_b_g)

## ---- fig.width=4, fig.height=3-----------------------------------------------
p <- plot_x(opt_cut_b_g)
p + ggtitle("Distribution of dsi") + theme_minimal() + xlab("Depression score")

## ---- fig.width=4, fig.height=3, cache=FALSE----------------------------------
plot_cutpointr(opt_cut_b, xvar = cutpoints, yvar = sum_sens_spec, conf_lvl = 0.9)
plot_cutpointr(opt_cut_b, xvar = fpr, yvar = tpr, aspect_ratio = 1, conf_lvl = 0)
plot_cutpointr(opt_cut_b, xvar = cutpoint, yvar = tp, conf_lvl = 0.9) + geom_point()

## ---- fig.width=4, fig.height=3-----------------------------------------------
opt_cut_b_g %>% 
    select(data, subgroup) %>% 
    unnest(cols = data) %>% 
    ggplot(aes(x = suicide, y = dsi)) + 
    geom_boxplot(alpha = 0.3) + facet_grid(~subgroup)

## ---- eval = FALSE------------------------------------------------------------
#  # Return cutpoint that maximizes the sum of sensitivity and specificiy
#  # ROCR package
#  rocr_sensspec <- function(x, class) {
#      pred <- ROCR::prediction(x, class)
#      perf <- ROCR::performance(pred, "sens", "spec")
#      sens <- slot(perf, "y.values")[[1]]
#      spec <- slot(perf, "x.values")[[1]]
#      cut <- slot(perf, "alpha.values")[[1]]
#      cut[which.max(sens + spec)]
#  }
#  
#  # pROC package
#  proc_sensspec <- function(x, class) {
#      r <- pROC::roc(class, x, algorithm = 2, levels = c(0, 1), direction = "<")
#      pROC::coords(r, "best", ret="threshold", transpose = FALSE)[1]
#  }

## ---- eval = FALSE, echo = FALSE----------------------------------------------
#  library(OptimalCutpoints)
#  library(ThresholdROC)
#  n <- 100
#  set.seed(123)
#  dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
#  x_pos <- dat$x[dat$y == 1]
#  x_neg <- dat$x[dat$y == 0]
#  bench_100 <- microbenchmark::microbenchmark(
#      cutpointr(dat, x, y, pos_class = 1, neg_class = 0,
#                direction = ">=", metric = youden, break_ties = mean),
#      rocr_sensspec(dat$x, dat$y),
#      proc_sensspec(dat$x, dat$y),
#      optimal.cutpoints(X = "x", status = "y", tag.healthy = 0, methods = "Youden",
#                        data = dat),
#      thres2(k1 = x_neg, k2 = x_pos, rho = 0.5,
#             method = "empirical", ci = FALSE),
#      times = 100, unit = "ms"
#  )
#  
#  n <- 1000
#  set.seed(123)
#  dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
#  x_pos <- dat$x[dat$y == 1]
#  x_neg <- dat$x[dat$y == 0]
#  bench_1000 <- microbenchmark::microbenchmark(
#      cutpointr(dat, x, y, pos_class = 1, neg_class = 0,
#                direction = ">=", metric = youden, break_ties = mean),
#      rocr_sensspec(dat$x, dat$y),
#      proc_sensspec(dat$x, dat$y),
#      optimal.cutpoints(X = "x", status = "y", tag.healthy = 0, methods = "Youden",
#                        data = dat),
#      thres2(k1 = x_neg, k2 = x_pos, rho = 0.5,
#             method = "empirical", ci = FALSE),
#      times = 100, unit = "ms"
#  )
#  
#  n <- 10000
#  set.seed(123)
#  dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
#  x_pos <- dat$x[dat$y == 1]
#  x_neg <- dat$x[dat$y == 0]
#  bench_10000 <- microbenchmark::microbenchmark(
#      cutpointr(dat, x, y, pos_class = 1, neg_class = 0,
#                direction = ">=", metric = youden, break_ties = mean, silent = TRUE),
#      rocr_sensspec(dat$x, dat$y),
#      optimal.cutpoints(X = "x", status = "y", tag.healthy = 0, methods = "Youden",
#                        data = dat),
#      proc_sensspec(dat$x, dat$y),
#      thres2(k1 = x_neg, k2 = x_pos, rho = 0.5,
#             method = "empirical", ci = FALSE),
#      times = 100
#  )
#  
#  n <- 1e5
#  set.seed(123)
#  dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
#  bench_1e5 <- microbenchmark::microbenchmark(
#      cutpointr(dat, x, y, pos_class = 1, neg_class = 0,
#                direction = ">=", metric = youden, break_ties = mean),
#      rocr_sensspec(dat$x, dat$y),
#      proc_sensspec(dat$x, dat$y),
#      times = 100, unit = "ms"
#  )
#  
#  n <- 1e6
#  set.seed(123)
#  dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
#  bench_1e6 <- microbenchmark::microbenchmark(
#      cutpointr(dat, x, y, pos_class = 1, neg_class = 0,
#                direction = ">=", metric = youden, break_ties = mean),
#      rocr_sensspec(dat$x, dat$y),
#      proc_sensspec(dat$x, dat$y),
#      times = 30, unit = "ms"
#  )
#  
#  n <- 1e7
#  set.seed(123)
#  dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
#  bench_1e7 <- microbenchmark::microbenchmark(
#      cutpointr(dat, x, y, pos_class = 1, neg_class = 0,
#                direction = ">=", metric = youden, break_ties = mean),
#      rocr_sensspec(dat$x, dat$y),
#      proc_sensspec(dat$x, dat$y),
#      times = 30, unit = "ms"
#  )
#  
#  results <- rbind(
#      data.frame(time = summary(bench_100)$median,
#                 Solution = summary(bench_100)$expr,
#                 n = 100),
#      data.frame(time = summary(bench_1000)$median,
#                 Solution = summary(bench_1000)$expr,
#                 n = 1000),
#      data.frame(time = summary(bench_10000)$median,
#                 Solution = summary(bench_10000)$expr,
#                 n = 10000),
#      data.frame(time = summary(bench_1e5)$median,
#                 Solution = summary(bench_1e5)$expr,
#                 n = 1e5),
#      data.frame(time = summary(bench_1e6)$median,
#                 Solution = summary(bench_1e6)$expr,
#                 n = 1e6),
#      data.frame(time = summary(bench_1e7)$median,
#                 Solution = summary(bench_1e7)$expr,
#                 n = 1e7)
#  )
#  results$Solution <- as.character(results$Solution)
#  results$Solution[grep(pattern = "cutpointr", x = results$Solution)] <- "cutpointr"
#  results$Solution[grep(pattern = "rocr", x = results$Solution)] <- "ROCR"
#  results$Solution[grep(pattern = "optimal", x = results$Solution)] <- "OptimalCutpoints"
#  results$Solution[grep(pattern = "proc", x = results$Solution)] <- "pROC"
#  results$Solution[grep(pattern = "thres", x = results$Solution)] <- "ThresholdROC"
#  
#  results$task <- "Cutpoint Estimation"

## ---- echo = FALSE------------------------------------------------------------
# These are the original results on our system
# dput(results)
results <- structure(list(time = c(4.5018015, 1.812802, 0.662101, 2.2887015, 
1.194301, 4.839401, 2.1764015, 0.981001, 45.0568005, 36.2398515, 
8.5662515, 5.667101, 2538.612001, 4.031701, 2503.8012505, 45.384501, 
43.118751, 37.150151, 465.003201, 607.023851, 583.0950005, 5467.332801, 
7850.2587, 7339.356101), Solution = c("cutpointr", "ROCR", "pROC", 
"OptimalCutpoints", "ThresholdROC", "cutpointr", "ROCR", "pROC", 
"OptimalCutpoints", "ThresholdROC", "cutpointr", "ROCR", "OptimalCutpoints", 
"pROC", "ThresholdROC", "cutpointr", "ROCR", "pROC", "cutpointr", 
"ROCR", "pROC", "cutpointr", "ROCR", "pROC"), n = c(100, 100, 
100, 100, 100, 1000, 1000, 1000, 1000, 1000, 10000, 10000, 10000, 
10000, 10000, 1e+05, 1e+05, 1e+05, 1e+06, 1e+06, 1e+06, 1e+07, 
1e+07, 1e+07), task = c("Cutpoint Estimation", "Cutpoint Estimation", 
"Cutpoint Estimation", "Cutpoint Estimation", "Cutpoint Estimation", 
"Cutpoint Estimation", "Cutpoint Estimation", "Cutpoint Estimation", 
"Cutpoint Estimation", "Cutpoint Estimation", "Cutpoint Estimation", 
"Cutpoint Estimation", "Cutpoint Estimation", "Cutpoint Estimation", 
"Cutpoint Estimation", "Cutpoint Estimation", "Cutpoint Estimation", 
"Cutpoint Estimation", "Cutpoint Estimation", "Cutpoint Estimation", 
"Cutpoint Estimation", "Cutpoint Estimation", "Cutpoint Estimation", 
"Cutpoint Estimation")), row.names = c(NA, -24L), class = "data.frame")

## ---- eval = FALSE------------------------------------------------------------
#  # ROCR package
#  rocr_roc <- function(x, class) {
#      pred <- ROCR::prediction(x, class)
#      perf <- ROCR::performance(pred, "sens", "spec")
#      return(NULL)
#  }
#  
#  # pROC package
#  proc_roc <- function(x, class) {
#      r <- pROC::roc(class, x, algorithm = 2, levels = c(0, 1), direction = "<")
#      return(NULL)
#  }

## ---- eval = FALSE, echo = FALSE----------------------------------------------
#  n <- 100
#  set.seed(123)
#  dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
#  bench_100 <- microbenchmark::microbenchmark(
#      cutpointr::roc(dat, "x", "y", pos_class = 1,
#                     neg_class = 0, direction = ">="),
#      rocr_roc(dat$x, dat$y),
#      proc_roc(dat$x, dat$y),
#      times = 100, unit = "ms"
#  )
#  n <- 1000
#  set.seed(123)
#  dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
#  bench_1000 <- microbenchmark::microbenchmark(
#      cutpointr::roc(dat, "x", "y", pos_class = 1, neg_class = 0,
#                     direction = ">="),
#      rocr_roc(dat$x, dat$y),
#      proc_roc(dat$x, dat$y),
#      times = 100, unit = "ms"
#  )
#  n <- 10000
#  set.seed(123)
#  dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
#  bench_10000 <- microbenchmark::microbenchmark(
#      cutpointr::roc(dat, "x", "y", pos_class = 1, neg_class = 0,
#                     direction = ">="),
#      rocr_roc(dat$x, dat$y),
#      proc_roc(dat$x, dat$y),
#      times = 100, unit = "ms"
#  )
#  n <- 1e5
#  set.seed(123)
#  dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
#  bench_1e5 <- microbenchmark::microbenchmark(
#      cutpointr::roc(dat, "x", "y", pos_class = 1, neg_class = 0,
#                     direction = ">="),
#      rocr_roc(dat$x, dat$y),
#      proc_roc(dat$x, dat$y),
#      times = 100, unit = "ms"
#  )
#  n <- 1e6
#  set.seed(123)
#  dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
#  bench_1e6 <- microbenchmark::microbenchmark(
#      cutpointr::roc(dat, "x", "y", pos_class = 1, neg_class = 0,
#                     direction = ">="),
#      rocr_roc(dat$x, dat$y),
#      proc_roc(dat$x, dat$y),
#      times = 30, unit = "ms"
#  )
#  n <- 1e7
#  set.seed(123)
#  dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
#  bench_1e7 <- microbenchmark::microbenchmark(
#      cutpointr::roc(dat, "x", "y", pos_class = 1, neg_class = 0,
#                     direction = ">="),
#      rocr_roc(dat$x, dat$y),
#      proc_roc(dat$x, dat$y),
#      times = 30, unit = "ms"
#  )
#  
#  results_roc <- rbind(
#      data.frame(time = summary(bench_100)$median,
#                 Solution = summary(bench_100)$expr,
#                 n = 100),
#      data.frame(time = summary(bench_1000)$median,
#                 Solution = summary(bench_1000)$expr,
#                 n = 1000),
#      data.frame(time = summary(bench_10000)$median,
#                 Solution = summary(bench_10000)$expr,
#                 n = 10000),
#      data.frame(time = summary(bench_1e5)$median,
#                 Solution = summary(bench_1e5)$expr,
#                 n = 1e5),
#      data.frame(time = summary(bench_1e6)$median,
#                 Solution = summary(bench_1e6)$expr,
#                 n = 1e6),
#      data.frame(time = summary(bench_1e7)$median,
#                 Solution = summary(bench_1e7)$expr,
#                 n = 1e7)
#  )
#  results_roc$Solution <- as.character(results_roc$Solution)
#  results_roc$Solution[grep(pattern = "cutpointr", x = results_roc$Solution)] <- "cutpointr"
#  results_roc$Solution[grep(pattern = "rocr", x = results_roc$Solution)] <- "ROCR"
#  results_roc$Solution[grep(pattern = "proc", x = results_roc$Solution)] <- "pROC"
#  results_roc$task <- "ROC curve calculation"

## ---- echo = FALSE------------------------------------------------------------
# Our results
results_roc <- structure(list(time = c(0.7973505, 1.732651, 0.447701, 0.859301, 
2.0358515, 0.694802, 1.878151, 5.662151, 3.6580505, 11.099251, 
42.8208515, 35.3293005, 159.8100505, 612.471901, 610.4337005, 
2032.693551, 7806.3854515, 7081.897251), Solution = c("cutpointr", 
"ROCR", "pROC", "cutpointr", "ROCR", "pROC", "cutpointr", "ROCR", 
"pROC", "cutpointr", "ROCR", "pROC", "cutpointr", "ROCR", "pROC", 
"cutpointr", "ROCR", "pROC"), n = c(100, 100, 100, 1000, 1000, 
1000, 10000, 10000, 10000, 1e+05, 1e+05, 1e+05, 1e+06, 1e+06, 
1e+06, 1e+07, 1e+07, 1e+07), task = c("ROC curve calculation", 
"ROC curve calculation", "ROC curve calculation", "ROC curve calculation", 
"ROC curve calculation", "ROC curve calculation", "ROC curve calculation", 
"ROC curve calculation", "ROC curve calculation", "ROC curve calculation", 
"ROC curve calculation", "ROC curve calculation", "ROC curve calculation", 
"ROC curve calculation", "ROC curve calculation", "ROC curve calculation", 
"ROC curve calculation", "ROC curve calculation")), row.names = c(NA, 
-18L), class = "data.frame")

## ---- echo = FALSE------------------------------------------------------------
results_all <- dplyr::bind_rows(results, results_roc)

ggplot(results_all, aes(x = n, y = time, col = Solution, shape = Solution)) +
  geom_point(size = 3) + geom_line() +
  scale_y_log10(breaks = c(0.5, 1, 2, 3, 5, 10, 25, 100, 250, 1000, 5000, 1e4, 15000)) +
  scale_x_log10(breaks = c(100, 1000, 1e4, 1e5, 1e6, 1e7)) +
  ylab("Median Time (milliseconds, log scale)") + xlab("Sample Size (log scale)") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.8, "cm"), 
        panel.spacing = unit(1, "lines")) +
  facet_grid(~task)

## ---- echo = FALSE------------------------------------------------------------
res_table <- tidyr::spread(results_all, Solution, time) %>% 
  arrange(task)
knitr::kable(res_table)

