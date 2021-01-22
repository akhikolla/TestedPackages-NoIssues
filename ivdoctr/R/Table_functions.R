#' @import data.table
NULL
#' Rounds x to two decimal places
#' @param x Number to be rounded
#' @return Number rounded to 2 decimal places
myformat <- function(x){
  x <- ifelse(is.na(x), NA, format(round(x, 2), n_digits = 2, nsmall = 2))
}

#' Creates LaTeX code for parameter estimates
#' @param est Number
#' @return LaTeX string for the number
format_est <- function(est) {
  paste0('$', myformat(est), '$')
}

#' Creates LaTeX code for the standard error
#' @param se Standard error
#' @return LaTeX string for the standard error
format_se <- function(se) {
  paste0('$(', myformat(se), ')$')
}

#' Creates LaTeX code for the HPDI
#' @param bounds 2-element vector of the upper and lower HPDI bounds
#' @return LaTeX string of the HPDI
format_HPDI <- function(bounds) {
  paste0('$[', myformat(bounds[1]), ', ', myformat(bounds[2]), ']$')
}

#' Makes LaTeX code to make a row of a table and shift by some amount of columns
#'   if necessary
#' @param char_vec Vector of characters to be collapsed into a LaTeX table
#' @param shift Number of columns to shift over
#' @return LaTeX string of the whole row of the table
make_tex_row <- function(char_vec, shift = 0) {
  out <- paste0(char_vec, collapse = ' & ')
  if (identical(shift, 0)) {
    out <- paste(out, '\\\\')
  } else {
    out <- paste(paste0(rep('&', shift), collapse = ''), out, '\\\\')
  }
  return(out)
}

#' Takes the OLS and IV estimates and converts it to a row of the LaTeX table
#' @param stats List with OLS and IV estimates and the bounds on kappa and r_uz
#' @param example_name Character string detailing the example
#' @return LaTeX code passed to makeTable()
make_full_row <- function(stats, example_name) {
  example_name <- paste(example_name, paste0('($n=', stats$n, '$)'))
  if (stats$binary == 1) {
    est <- with(stats, sapply(c(b_OLS, b_IV, a0, a1, psi_lower), format_est))
    se <- with(stats, c(format_se(se_OLS), format_se(se_IV), "", "", ""))
  } else {
    est <- with(stats, sapply(c(b_OLS, b_IV, k_lower), format_est))
    se <- with(stats, c(format_se(se_OLS), format_se(se_IV), ""))
  }

  if (stats$p_empty > 0) {
    panel_II_III <- c(format_est(stats$p_empty), "$-$", "$[-, -]$", "$[-, -]$",
                      "$-$", "$-$")
    se <- c(se, "", "", "", "", "$[-, -]$", "$[-, -]$")
  } else {
    panel_II_III <- with(stats, c(format_est(p_empty), format_est(p_valid),
                                  format_HPDI(r_uz_full_interval),
                                  format_HPDI(beta_full_interval),
                                  format_est(r_uz_median), format_est(beta_median)))
    se <- with(stats, c(se, "", "", "", "", format_HPDI(r_uz_partial_interval),
                        format_HPDI(beta_partial_interval)))
  }

  est <- make_tex_row(c(example_name, est, panel_II_III))
  scenario <- with(stats, paste0("\\quad $(\\kappa, \\rho_{u\\xi^*}) \\in (",
                                 ifelse(min(k_restriction) < 0.01, 0, min(k_restriction)),
                                 ",", max(k_restriction), "] \\times [",
                                 min(r_TstarU_restriction), ",",
                                 max(r_TstarU_restriction), "]$"))
  se <- make_tex_row(c(scenario, se))
  paste(est, se)
}

#' Generates parameter estimates given user restrictions and data
#'
#' @param y_name Character string with the column name of the dependent variable
#' @param T_name Character string with the column name of the endogenous regressor(s)
#' @param z_name Character string with the column name of the instrument(s)
#' @param data Data frame
#' @param example_name Character string naming estimation
#' @param controls Vector of character strings specifying the exogenous variables
#' @param robust Indicator for heteroskedasticity-robust standard errors
#' @param r_TstarU_restriction 2-element vector of min and max of r_TstarU.
#' @param k_restriction 2-element vector of min and max of kappa.
#' @param n_draws Number of draws when generating frequentist-friendly draws of the covariance matrix
#' @param n_RF_draws Number of reduced-form draws
#' @param n_IS_draws Number of draws on the identified set
#' @param resample Indicator of whether or not to resample using magnification factor
#' @return List with elements:
#'\itemize{
#'  \item ols: lm object of OLS estimation,
#'  \item iv: ivreg object of the IV estimation
#'  \item n: Number of observations
#'  \item b_OLS: OLS point estimate
#'  \item se_OLS: OLS standard errors
#'  \item b_IV: IV point estimate
#'  \item se_IV: IV standard errors
#'  \item k_lower: lower bound of kappa
#'  \item p_empty: fraction of parameter draws that yield an empty identified set
#'  \item p_valid: fraction of parameter draws compatible with a valid instrument
#'  \item r_uz_full_interval: 90\% posterior credible interval for fully identified set of rho
#'  \item beta_full_interval: 90\% posterior credible interval for fully identified set of beta
#'  \item r_uz_median: posterior median for partially identified rho
#'  \item r_uz_partial_interval: 90\% posterior credible interval for partially identified set of rho under a conditionally uniform reference prior
#'  \item beta_median: posterior median for partially identified beta
#'  \item beta_partial_interval: 90\% posterior credible interval for partially identified set of beta under a conditionally uniform reference prior
#'  \item a0: If treatment is binary, mis-classification probability of no-treatment case. NULL otherwise
#'  \item a1: If treatment is binary, mis-classification probability of treatment case. NULL otherwise
#'  \item psi_lower: lower bound for psi
#'  \item binary: logical indicating if treatment is binary
#'  \item k_restriction: User-specified bounds on kappa
#'  \item r_TstarU_restriction: User-specified bounds on r_TstarU
#'  }
#'@examples
#'library(ivdoctr)
#'endog <- c(0, 0.9)
#'meas <- c(0.6, 1)
#'
#'colonial_example1 <- ivdoctr(y_name = "logpgp95", T_name = "avexpr",
#'                             z_name = "logem4", data = colonial,
#'                             controls = NULL, robust = FALSE,
#'                             r_TstarU_restriction = endog,
#'                             k_restriction = meas,
#'                             example_name = "Colonial Origins")
#'
#' @export
ivdoctr <- function(y_name, T_name, z_name, data, example_name,
                    controls = NULL, robust = FALSE,
                    r_TstarU_restriction = c(-1, 1), k_restriction = c(0.0001, 1),
                    n_draws = 5000, n_RF_draws = 1000, n_IS_draws = 1000,
                    resample = FALSE) {
  binary <- ifelse(uniqueN(data[[T_name]]) == 2, 1, 0)

  summary_stats <- get_estimates(y_name, T_name, z_name, data, controls, robust)
  obs <- get_observables(y_name, T_name, z_name, data, controls)
  bounds_unrest <- get_bounds_unrest(obs)

  if (binary) {
    p <- mean(data[[T_name]])
    alpha_bounds <- get_alpha_bounds(obs, p)
    psi_lower <- get_psi_lower(obs$s2_T, p, bounds_unrest$k$Lower)
  } else {
    alpha_bounds <- NULL
    psi_lower <- NULL
  }

  # Running through simulations
  bounds <- draw_bounds(y_name, T_name, z_name, data, controls,
                        r_TstarU_restriction, k_restriction, n_draws)
  freq <- summarize_bounds(bounds)
  posterior <- draw_posterior(y_name, T_name, z_name, data, controls,
                              r_TstarU_restriction, k_restriction,
                              n_RF_draws, n_IS_draws, resample)
  if (binary) {
    bayes <- summarize_posterior_binary(posterior, p)
  } else {
    bayes <- summarize_posterior(posterior)
  }

  # Compute covering beta interval
  if (binary) {
    beta_bounds <- get_beta_bounds_binary_post(posterior, n_RF_draws)
    beta_center <- posterior$beta_center
  } else {
    beta_center <- bounds$beta_center
    beta_bounds <- cbind(bounds$restricted$beta_lower, bounds$restricted$beta_upper)
  }
  beta_interval <- getInterval(beta_bounds, beta_center)

  # Compute covering r_uz interval
  r_uz_center <- bounds$r_uz_center
  r_uz_bounds <- cbind(bounds$restricted$r_uz_lower, bounds$restricted$r_uz_upper)
  r_uz_interval <- getInterval(r_uz_bounds, r_uz_center)
  r_uz_interval <- c(max(-1, r_uz_interval[1]), min(1, r_uz_interval[2]))


  full_stats <- list(ols = summary_stats$ols,
                     iv = summary_stats$iv,
                     n = summary_stats$n,
                     b_OLS = summary_stats$b_OLS,
                     se_OLS = summary_stats$se_OLS,
                     b_IV = summary_stats$b_IV,
                     se_IV = summary_stats$se_IV,
                     k_lower = bounds_unrest$k$Lower,
                     a0 = alpha_bounds$a0,
                     a1 = alpha_bounds$a1,
                     psi_lower = psi_lower,
                     binary = binary,
                     p_empty = freq$p_empty,
                     p_valid = freq$p_valid,
                     r_uz_full_interval = r_uz_interval,
                     beta_full_interval = beta_interval,
                     r_uz_median = bayes$HPDI$median[1],
                     r_uz_partial_interval = c(bayes$HPDI$lower[1], bayes$HPDI$upper[1]),
                     beta_median = bayes$HPDI$median[2],
                     beta_partial_interval = c(bayes$HPDI$lower[2], bayes$HPDI$upper[2]),
                     k_restriction = k_restriction,
                     r_TstarU_restriction = r_TstarU_restriction,
                     example_name = example_name)
  return(full_stats)
}

# Generates header LaTeX code for continuous table
table_header_cts <- function() {
  "\\begin{tabular}{lccccccccc}
  \\hline
  \\hline
  &\\multicolumn{3}{c}{(I) Summary Statistics}
  &\\multicolumn{4}{c}{(II) Inference for $\\Theta$}
  &\\multicolumn{2}{c}{(III) Inference for $\\theta$} \\\\
  \\cmidrule(lr){2-4}\\cmidrule(lr){5-8}\\cmidrule(lr){9-10}
  & OLS & IV & $L$ & $\\mathbb{P}(\\varnothing)$ & $\\mathbb{P}(\\mbox{Valid})$ & $\\rho_{u \\zeta}$ & $\\beta$ & $\\rho_{u \\zeta}$ & $\\beta$ \\\\"
}

# Generates header LaTeX code for continuous table
table_header_bin <- function() {
  "\\begin{tabular}{lccccccccccc}
  \\hline
  \\hline
  &\\multicolumn{5}{c}{(I) Summary Statistics}
  &\\multicolumn{4}{c}{(II) Inference for $\\Theta$}
  &\\multicolumn{2}{c}{(III) Inference for $\\theta$} \\\\
  \\cmidrule(lr){2-6}\\cmidrule(lr){7-10}\\cmidrule(lr){11-12}
  & OLS & IV & $\\bar{\\alpha_0}$ & $\\bar{\\alpha_1}$ & \\underbar{$\\psi$} & $\\mathbb{P}(\\varnothing)$ & $\\mathbb{P}(\\mbox{Valid})$ & $\\rho_{u \\zeta}$ & $\\beta$ & $\\rho_{u \\zeta}$ & $\\beta$ \\\\
  \\\\"
}

# Generates footer LaTeX code for continuous table
table_footer_fn <- function() {
  "\\hline
  \\end{tabular}"
}

#' Generates table of parameter estimates given user restrictions and data
#'
#' @param ... Arguments of TeX code for individual examples to be combined into a single table
#' @param output File name to write
#' @return LaTeX code that generates output table with regression results
#'
#'@examples
#'library(ivdoctr)
#'endog <- c(0, 0.9)
#'meas <- c(0.6, 1)
#'
#'colonial_example1 <- ivdoctr(y_name = "logpgp95", T_name = "avexpr",
#'                             z_name = "logem4", data = colonial,
#'                             controls = NULL, robust = FALSE,
#'                             r_TstarU_restriction = endog,
#'                             k_restriction = meas,
#'                             example_name = "Colonial Origins")
#'makeTable(colonial_example1, output = file.path(tempdir(), "colonial.tex"))
#'
#' @export
makeTable <- function(..., output) {

  table_examples <- c()
  list_examples <- list(...)

  for (i in 1:length(list_examples)) {
    table_examples[i] <- make_full_row(list_examples[[i]],
                                       list_examples[[i]]$example_name)
  }

  if (list_examples[[1]]$binary) {
    cat(table_header_bin(), '\\\\', table_examples, table_footer_fn(), sep = '\n')
    cat(table_header_bin(), '\\\\', table_examples, table_footer_fn(), sep = '\n', file = output)
  } else {
    cat(table_header_cts(), '\\\\', table_examples, table_footer_fn(), sep = '\n')
    cat(table_header_cts(), '\\\\', table_examples, table_footer_fn(), sep = '\n', file = output)
  }
}
