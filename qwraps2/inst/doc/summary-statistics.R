## ----label = "setup", include = FALSE-----------------------------------------
knitr::opts_chunk$set(collapse = TRUE)

## -----------------------------------------------------------------------------
set.seed(42)
  library(qwraps2)
# define the markup language we are working in.
# options(qwraps2_markup = "latex") is also supported.
options(qwraps2_markup = "markdown")

## -----------------------------------------------------------------------------
library(qwraps2)

## -----------------------------------------------------------------------------
data(mtcars2)
str(mtcars2)

## -----------------------------------------------------------------------------
mean_sd(mtcars2$mpg)
mean_sd(mtcars2$mpg, denote_sd = "paren")

## -----------------------------------------------------------------------------
mci <- mean_ci(mtcars2$mpg)
str(mci)
mci
print(mci, show_level = TRUE)

## -----------------------------------------------------------------------------
median_iqr(mtcars2$mpg)

## -----------------------------------------------------------------------------
n_perc(mtcars2$cyl == 4)
n_perc0(mtcars2$cyl == 4)

n_perc(mtcars2$cyl_factor == 4)  # this returns 0 (0.00%)
n_perc(mtcars2$cyl_factor == "4 cylinders")
n_perc(mtcars2$cyl_factor == levels(mtcars2$cyl_factor)[2])

# The count and percentage of 4 or 6 cylinders vehicles in the data set is
n_perc(mtcars2$cyl %in% c(4, 6))

## -----------------------------------------------------------------------------
x <- runif(6, min = 4, max = 70)

# geometric mean
mu_g <- prod(x) ** (1 / length(x))
mu_g
exp(mean(log(x)))
1.2 ** mean(log(x, base = 1.2))

# geometric standard deviation
exp(sd(log(x)))  ## This is wrong

# these equations are correct
sigma_g <- exp(sqrt(sum(log(x / mu_g) ** 2) / length(x)))
sigma_g

exp(sqrt((length(x) - 1) / length(x)) * sd(log(x)))

## -----------------------------------------------------------------------------
gmean(x)
all.equal(gmean(x), mu_g)

gvar(x)
all.equal(gvar(x), sigma_g^2)  # This is supposed to be FALSE
all.equal(gvar(x), exp(log(sigma_g)^2))

gsd(x)
all.equal(gsd(x), sigma_g)

## -----------------------------------------------------------------------------
gmean_sd(x)

## -----------------------------------------------------------------------------
our_summary1 <-
  list("Miles Per Gallon" =
       list("min"       = ~ min(mpg),
            "max"       = ~ max(mpg),
            "mean (sd)" = ~ qwraps2::mean_sd(mpg)),
       "Displacement" =
       list("min"       = ~ min(disp),
            "median"    = ~ median(disp),
            "max"       = ~ max(disp),
            "mean (sd)" = ~ qwraps2::mean_sd(disp)),
       "Weight (1000 lbs)" =
       list("min"       = ~ min(wt),
            "max"       = ~ max(wt),
            "mean (sd)" = ~ qwraps2::mean_sd(wt)),
       "Forward Gears" =
       list("Three" = ~ qwraps2::n_perc0(gear == 3),
            "Four"  = ~ qwraps2::n_perc0(gear == 4),
            "Five"  = ~ qwraps2::n_perc0(gear == 5))
       )

## ----results = "asis"---------------------------------------------------------
### Overall
whole <- summary_table(mtcars2, our_summary1)
whole

## ----results = "asis"---------------------------------------------------------
### By number of Cylinders
by_cyl <- summary_table(dplyr::group_by(mtcars2, cyl_factor), our_summary1)
by_cyl

## ----results = "asis"---------------------------------------------------------
summary_table(mtcars2, summaries = our_summary1, by = c("cyl_factor"))

## -----------------------------------------------------------------------------
by_cyl_am <- summary_table(mtcars2, summaries = our_summary1, by = c("cyl_factor", "am"))
by_cyl_am

## -----------------------------------------------------------------------------
all.equal(summary_table(dplyr::group_by(mtcars2, cyl_factor, am), summaries = our_summary1),
          by_cyl_am)

## ----results = "asis"---------------------------------------------------------
summary_table(dplyr::group_by(mtcars2, carb), summaries = our_summary1, by = c("cyl_factor", "am"))

## ----results = "asis"---------------------------------------------------------
both <- cbind(whole, by_cyl)
both

## ----results = "asis"---------------------------------------------------------
print(both,
      rtitle = "Summary Statistics",
      cnames = c("Col 0", "Col 1", "Col 2", "Col 3"))

## -----------------------------------------------------------------------------
qsummary(mtcars2[, c("mpg", "cyl_factor", "wt")])

## ----label="summary_table_mtcars2_default", results = "asis"------------------
summary_table(mtcars2[, c("mpg", "cyl_factor", "wt")])

## -----------------------------------------------------------------------------
new_summary <-
  qsummary(mtcars2[, c("mpg", "cyl_factor", "wt")],
           numeric_summaries = list("Minimum" = "~ min(%s)",
                                    "Maximum" = "~ max(%s)"),
           n_perc_args = list(digits = 1, show_symbol = TRUE, show_denom = "always"))
str(new_summary)

## ----results = "asis"---------------------------------------------------------
summary_table(mtcars2, new_summary)

## ----results = "asis"---------------------------------------------------------
summary_table(mtcars2, new_summary, by = c("cyl_factor"))

## -----------------------------------------------------------------------------
str(both)

## -----------------------------------------------------------------------------
# difference in means
mpvals <-
  sapply(
         list(lm(mpg ~ cyl_factor,  data = mtcars2),
              lm(disp ~ cyl_factor, data = mtcars2),
              lm(wt ~ cyl_factor,   data = mtcars2)),
         extract_fpvalue)

# Fisher test
fpval <- frmtp(fisher.test(table(mtcars2$gear, mtcars2$cyl_factor))$p.value)

## -----------------------------------------------------------------------------
both <- cbind(both, "P-value" = "")
both[grepl("mean \\(sd\\)", rownames(both)), "P-value"] <- mpvals
a <- capture.output(print(both))
a[grepl("Forward Gears", a)] <-
  sub("&nbsp;&nbsp;\\ \\|$", paste(fpval, "|"), a[grepl("Forward Gears", a)])

## ----results = "asis"---------------------------------------------------------
cat(a, sep = "\n")

## ----results = "asis"---------------------------------------------------------
gear_summary <-
  list("Forward Gears" =
       list("Three" = ~ qwraps2::n_perc0(gear == 3),
            "Four"  = ~ qwraps2::n_perc0(gear == 4),
            "Five"  = ~ qwraps2::n_perc0(gear == 5)),
       "Transmission" =
       list("Automatic" = ~ qwraps2::n_perc0(am == 0),
            "Manual"    = ~ qwraps2::n_perc0(am == 1))
       )

gear_summary <-
setNames(gear_summary,
         c(
         paste("Forward Gears: ", frmtp(fisher.test(xtabs( ~ gear + cyl_factor, data = mtcars2))$p.value)),
         paste("Transmission: ",  frmtp(fisher.test(xtabs( ~ am + cyl_factor, data = mtcars2))$p.value)))
         )

summary_table(mtcars2, gear_summary, by = "cyl_factor")

## -----------------------------------------------------------------------------
new_data_frame <-
  data.frame(age = c(18, 20, 24, 17, 43),
             edu = c(1, 3, 1, 5, 2),
             rt  = c(0.01, 0.04, 0.02, 0.10, 0.06))

# Set a label for the variables
attr(new_data_frame$age, "label") <- "Age in years"
attr(new_data_frame$rt,  "label") <- "Reaction time"

# mistakenly set the attribute to name instead of label
attr(new_data_frame$edu, "name") <- "Education"

## -----------------------------------------------------------------------------
qsummary(new_data_frame)

## ----results = "asis"---------------------------------------------------------
summary_table(new_data_frame)

## -----------------------------------------------------------------------------
print(sessionInfo(), local = FALSE)

