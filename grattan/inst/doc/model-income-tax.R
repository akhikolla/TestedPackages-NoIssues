## ----no-run-vignette, include=FALSE-------------------------------------------
library(knitr)
if (!identical(Sys.getenv("R_GRATTAN_BUILD_MAIN_VIGNETTE"), "true")) {
  do_eval <- function() return(FALSE)
  opts_hooks$set(inline = function(x) invisible(NULL))
  knit_hooks$set(inline = function(x) deparse(substitute(x)),
                 evaluate.inline = function(code, envir) {
                   deparse(substitute(code))
                 })
  
  opts_knit$set(eval = FALSE, error = TRUE, warning = FALSE)
  opts_chunk$set(eval = FALSE, error = TRUE, warning = FALSE)
}

## ----no-includes, include=FALSE-----------------------------------------------
#  do_eval <- function() {
#    library(hutils)
#    NEITHER(nzchar(Sys.getenv("CIRCLECI")),
#            nzchar(Sys.getenv("TRAVIS")))
#  }

## ----loadPackages-------------------------------------------------------------
#  library(knitr)
#  library(data.table)
#  library(magrittr)
#  library(hutils)
#  library(grattan)
#  require_taxstats1516()
#  
#  # Use the actual sample file if you've got it
#  s1516 <- as.data.table(sample_file_1516_synth)
#  s1516[, WEIGHT := 50L]

## ----detach-taxstats, include=FALSE-------------------------------------------
#  # memory
#  detach("package:taxstats1516", unload = TRUE)

## ----dollar-------------------------------------------------------------------
#  #' @return Number formatted as dollar e.g. 30e3 => $30,000
#  dollar <- function (x, digits = 0) {
#    nsmall <- digits
#    commaz <- format(abs(x), nsmall = nsmall, trim = TRUE, big.mark = ",",
#                     scientific = FALSE, digits = 1L)
#    if_else(x < 0,
#            paste0("\U2212","$", commaz),
#            paste0("$", commaz))
#  }

## ----baseline-fy--------------------------------------------------------------
#  s1516 %>%
#    model_income_tax(baseline_fy = "2015-16") %>%
#    select_grep("tax$", "Taxable_Income") %>%  # just look at relevant cols
#    head %>%
#    kable

## ----baseline-fy-int----------------------------------------------------------
#  is_all_equal <- function(x, y) {
#    if (is.integer(x) && is.integer(y)) {
#      all(x == y)
#    } else {
#      isTRUE(all.equal(x, y))
#    }
#  }
#  
#  s1516 %>%
#    model_income_tax(baseline_fy = "2015-16",
#                     return. = "sample_file.int") %>%
#    select_grep("tax$", "Taxable_Income") %T>%
#    .[, stopifnot(is_all_equal(baseline_tax, new_tax))] %>%
#    head %>%
#    kable

## ----s1516_no_changes---------------------------------------------------------
#  s1516_no_changes <-
#    # Temp budget repair levy not refundable against SBTO
#    s1516 %>%
#    model_income_tax(baseline_fy = "2015-16",
#                     ordinary_tax_thresholds = c(0, 18200, 37000, 80000, 180000),
#                                                                 # temp budget
#                                                                 # repair levy
#                     ordinary_tax_rates = c(0, 0.19, 0.325, 0.37, 0.45 + 0.02),
#                     return. = "sample_file.int")

## ----medicare-levy-rate-increase-a, eval=do_eval()----------------------------
#  m1516a <-
#    s1516 %>%
#    model_income_tax("2015-16",
#                     # Increase to 3%
#                     medicare_levy_rate = 0.03)

## ----medicare-levy-rate-increase-b, eval=do_eval()----------------------------
#  m1516a <-
#    s1516 %>%
#    model_income_tax("2015-16",
#                     # Increase to 3%
#                     medicare_levy_rate = 0.03,
#                     medicare_levy_upper_threshold = 30479,
#                     medicare_levy_upper_sapto_threshold = 48197)

## ----medicare-levy-rate-increase-c, eval=do_eval()----------------------------
#  m1516a <-
#    s1516 %>%
#    model_income_tax("2015-16",
#                     # Increase to 3%
#                     medicare_levy_rate = 0.03,
#                     # but keep the upper threshold the same
#                     medicare_levy_upper_threshold = 26670,
#                     medicare_levy_upper_sapto_threshold = 48197)

## ----medicare-levy-rate-increase-d, eval=do_eval()----------------------------
#  m1516a <-
#    s1516 %>%
#    model_income_tax("2015-16",
#                     # Increase to 3%
#                     medicare_levy_rate = 0.03,
#                     # but keep the upper threshold the same
#                     medicare_levy_lower_threshold = 21335,
#                     medicare_levy_upper_threshold = 26670,
#                     medicare_levy_upper_sapto_threshold = 48197)

## ----lito-a-------------------------------------------------------------------
#  L1516a <-
#    s1516 %>%
#    model_income_tax("2015-16",
#                     lito_max_offset = 1000)
#  revenue_foregone(L1516a)

## ----cleanup-before-project, include=FALSE------------------------------------
#  # reduce memory usage (e.g. travis, CCI)
#  rm(L1516a, m1516a, s1516_no_changes)

## ----project-a----------------------------------------------------------------
#  s1819 <- project(s1516, h = 3L)

## ----project-wage-lf-series---------------------------------------------------
#  s1819_lf2pc_wage2pc <-
#    s1516 %>%
#    project(h = 3L,
#            lf.series = 0.02,
#            wage.series = 0.02)

## ----compare-2pc-to-default---------------------------------------------------
#  tax_Grattan_1819 <-
#    s1819 %$%
#    income_tax(Taxable_Income, "2018-19", .dots.ATO = copy(s1819)) %>%
#    sum %>%
#    # Weight (equi-weighted so do now)
#    multiply_by(s1819[["WEIGHT"]][1L])
#  tax_2pc_1819 <-
#    s1819_lf2pc_wage2pc %$%
#    income_tax(Taxable_Income, "2018-19", .dots.ATO = copy(s1819)) %>%
#    sum %>%
#    # Weight (equi-weighted so do now)
#    multiply_by(s1819[["WEIGHT"]][1L])

## ----s2021--------------------------------------------------------------------
#  s2021 <- project(s1516, h = 5L)

## ----s1819_wage80pc-----------------------------------------------------------
#  s2021_wage80pc <-
#    s1516 %>%
#    copy %>%
#    .[, Sw_amt := wage_inflator(Sw_amt,
#                                from_fy = "2015-16",
#                                to_fy = "2020-21",
#                                forecast.level = 80,
#                                forecast.series = "upper")] %>%
#    .[] %>%
#    project(h = 5L,
#            excl_vars = "Sw_amt",
#            .copyDT = FALSE) %>%  # just for memory frugality
#    .[]

## ----compare-Sw_amt-----------------------------------------------------------
#  s2021[, mean(Sw_amt)] %>% dollar
#  s2021_wage80pc[, mean(Sw_amt)] %>% dollar
#  s2021[, mean(Taxable_Income)] %>% dollar
#  s2021_wage80pc[, mean(Taxable_Income)] %>% dollar

## ----cgt_25pc_fwd_estimates, eval=do_eval()-----------------------------------
#  cgt_25pc_fwd_estimates <-
#    lapply(yr2fy(2019:2022), function(fy) {
#      s1516 %>%
#        project_to(to_fy = fy) %>%
#        model_income_tax("2018-19",
#                         cgt_discount_rate = 0.25) %>%
#        .[, fy_year := fy]
#    }) %>%
#    rbindlist

## ----cgt_25pc_fwd_estimates-deciles, eval=do_eval()---------------------------
#  cgt_25pc_fwd_estimates %>%
#    mutate_ntile("Taxable_Income", n = 5L, keyby = "fy_year") %>%
#    .[, delta := new_tax - baseline_tax] %>%
#    .[, .(totDelta = sum(delta),
#          avgDelta = mean(delta)),
#      keyby = .(fy_year, Taxable_IncomeQuintile)] %>%
#    # cosmetic
#    .[, lapply(.SD, round), keyby = key(.)] %>%
#    kable

## ----lito_multi_201516--------------------------------------------------------
#  s1516 %>%
#    model_income_tax("2015-16",
#                     lito_multi = list(x = c(-Inf, 37e3, 200e3/3, Inf),
#                                       y = c(445, 445, 0, 0)),
#                     return. = "sample_file.int") %>%
#    .[new_tax != baseline_tax]
#  

## ----sapto_abolished1819------------------------------------------------------
#  sapto_abolished1819 <-
#    s1819 %>%
#    model_income_tax("2018-19",
#                     sapto_eligible = FALSE)

## ----sapto_abolished_abv27k_1819----------------------------------------------
#  sapto_abolished_abv27k_1819 <-
#    s1819 %>%
#    model_income_tax("2018-19",
#                     sapto_lower_threshold = 27000)

## ----sapto-age-of-entitlement-------------------------------------------------
#  s1718_AgeOfEntitlement <-
#      project(s1516,
#              h = 2L) %>%
#      model_income_tax("2017-18",
#                       sapto_lower_threshold = 27e3,
#                       sapto_lower_threshold_married = 42e3,
#                       sapto_max_offset = 1160,
#                       sapto_max_offset_married = 390,
#                       medicare_levy_lower_sapto_threshold = 27000,
#                       medicare_levy_upper_sapto_threshold = 33750,
#                       medicare_levy_upper_family_threshold = 46361,
#                       medicare_levy_lower_family_sapto_threshold = 42000,
#                       medicare_levy_upper_family_sapto_threshold = 52500)
#  revenue_foregone(s1718_AgeOfEntitlement)

