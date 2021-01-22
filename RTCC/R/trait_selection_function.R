#' Trait selection
#'
#' This function determines whether the selected traits exhibit or not a clustering/overdispersion
#' signal on the tested samples. For each trait, compares the observed Mean Pairwise Distance (MPD)
#' of each sample against a distribution of synthetic commmunities MPDs obtained by a randomization
#' test. Each synthetic community is build maintaining the original sample richness and randomly
#' selecting organisms form the global pool.
#'
#' @param table1 A data frame containing organisms names on the first column and its trait values on the
#' consecutive ones. It also has to contain two columns with the maximum and the minimum values of the tested
#' environmental variable where the organisms have been observed.
#'
#' @param table2 A presence-absence observations table with the organisms names on the first column and the
#' sample names as consecutive colnames.
#'
#' @param table3 A dataframe containing sample names on the first column and environmental parameters on the
#' consecutive ones.
#'
#' @param traits_columns Table 1 column numbers where different trait values appear.
#'
#'
#' @param repetitions Number of simulated synthetic communities distributions.
#'
#' @return The function returns a dataframe with trait names as colnames and the p-value distribution of the different
#' traits.
#'
#' @examples
#'
#' \donttest{
#' data(group_information)
#' data(table_presence_absence)
#' data(metadata)
#' rtcc1(group_information, table_presence_absence, metadata, 2:11, 100)
#' }
#'
#'
#' @export


rtcc1 <- function(table1, table2, table3, traits_columns, repetitions){

  dataset <- table1
  dataset <- as.data.frame(dataset)
  table_presence_absence <- table2
  metadata <- table3

  colnames(dataset)[1] <- "tax"

  colnames(metadata)[1] <- "sample_ID"

  colnames(table_presence_absence)[1] <- "tax"

  pool <- as.vector(dataset$tax)
  l_pool <- length(pool)
  richnes <- vegan::specnumber(t(table_presence_absence[,-1]))

  local_communities <- metadata$sample_ID

  y <- 1
  h_iteration <- 1
  results <- as.data.frame(matrix(nrow = length(local_communities), ncol = length(traits_columns), NA))
  colnames(results) <- colnames(dataset)[traits_columns]

  for(i in traits_columns){

    # Calculate MPDs for observed communities
    real_MPDs <- real_mpd(local_communities, table_presence_absence, dataset, i)

    trait <- as.numeric(dataset[,i])
    p <- 1
    species_presence <- rowSums(table_presence_absence[,colnames(table_presence_absence) %in% local_communities])
    p_values <- rep(NA, length(local_communities))
    species_presence[species_presence != 0] <- 1
    prob_vector <- species_presence
    for(m in local_communities){
      MPDs <- rep(NA, (repetitions))
      rich_m <- richnes[m]
      ncomps <- rich_m * (rich_m-1)
      MPDs <- sample_and_mpd(repetitions, h_iteration, trait, l_pool, rich_m, prob_vector, MPDs, ncomps)
      p_values[p] <- sum(MPDs < real_MPDs[p,2])/length(MPDs)
      p <- p + 1
    }
    results[,y] <- p_values
    y <- y + 1
  }
  return(results)
}
