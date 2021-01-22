\name{Spearman}
\alias{Spearman}

\title{Spearman test for dependence}

\description{
  Test statistic is computed by cor.test\{stats\}. See \link[stats]{cor.test}.
  Note that test statistic returned is the pivot t that approximately follows Student's t distribution.
  Spearman test cannot handle tie. Since bootstrap resamples with replacement which generates ties, bootstrap confidnece interval does not apply. Setting BS.CI > 0 throughs warning message.
}
