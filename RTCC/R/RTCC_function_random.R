#' Clustering signal significance.
#'
#' For a given trait and environmental variable, this function creates a null model of the
#' clustering/overdispersion pattern in order to test if the observed pattern statistically
#' differs from the expected by random. For this, it sequentially remove random samples from
#' the metacommunity and computes at each step the remaining metacommunity h-index. This index
#' is based on the percentage of samples on a metacoomunity presenting significant trait
#' clustering/overdispersion. After h iterations, computes a 95% confidence interval of the
#' obtained h-index for each point of the environmental gradient.
#'
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
#' @param species_abundances A vector containing the relative abundance of the organisms on the whole data set on
#' the same order as appear on Table 1.
#'
#' @param trait_col_number Table 1 column number of the tested trait.
#'
#' @param min_env_col Table 1 column number indicating the minimum value of the environmental variable were each
#' organism has been observed.
#'
#' @param max_env_col Table 1 column number indicating the maximum value of the environmental variable were each
#' organism has been observed.
#'
#' @param env_var_col Table 2 column number indicating the tested environmental variable.
#'
#' @param h_iteration Number of h-index calculations for computing a confidence interval.

#' @param repetitions Number of simulated synthetic communities distributions.
#'
#' @param model Model selection. All models build synthetic communities based on the organisms richness of the
#' observed communities.
#'
#' - Model 1: organism are selected randomly from the global pool.
#' - Model 2: organism are selected randomly with a probability based on its relative abundance on the global pool.
#' - Model 3: organism are selected randomly, but only those whose environmental range includes the value of
#' the simulated community are elegible.
#' - Model 4: organism are selected randomly, but only those whose environmental range includes the value of
#' the simulated community are elegible and the selection probability is based on its relative abundance on the global pool.
#'
#' @return The function returns a dataframe with the maximum value of environmental variable corresponding to the same
#' number of samples on the ordered remova, h-index calculation for each environmental value, and the percentiles 0.025,
#' 0.5 and 0.975 of the obtained distribution for each point (mean value and 95% confidence interval).
#'
#' @examples
#'
#' \donttest{
#' data(group_information)
#' data(table_presence_absence)
#' data(metadata)
#' rtcc3(group_information, table_presence_absence, metadata, group_information$sums,
#' 9, 12, 13, 2, 50, 20, model = 1)
#' }
#'
#' @export


rtcc3 <- function(table1, table2, table3, species_abundances, trait_col_number, min_env_col,
                  max_env_col, env_var_col, h_iteration, repetitions, model){

  dataset <- table1
  dataset <- as.data.frame(dataset)
  table_presence_absence <- table2
  metadata <- as.data.frame(table3)

  colnames(dataset)[min_env_col] <- "min"
  colnames(dataset)[max_env_col] <- "max"
  colnames(dataset)[1] <- "tax"

  colnames(metadata)[1] <- "sample_ID"
  colnames(metadata)[env_var_col] <- "env_variable"

  colnames(table_presence_absence)[1] <- "tax"
  table_presence_absence <- table_presence_absence[order(table_presence_absence$tax),]

  pool <- as.vector(dataset$tax)
  l_pool <- length(pool)
  richnes <- vegan::specnumber(t(table_presence_absence[,-1]))
  trait <- dataset[,trait_col_number]

  metadata <- metadata[order(metadata$env_variable, decreasing = T),]
  local_communities <- metadata$sample_ID

  # Calculate MPDs for observed communities
  i <- trait_col_number
  real_MPDs <- real_mpd(local_communities, table_presence_absence, dataset, i)

  max_sal <- rep(NA, length(local_communities))
  h_indexes <- rep(NA, length(local_communities))
  h_sd <- rep(NA, length(local_communities))
  observed_MPDs <- real_MPDs

  if (model == 1) { # Model I: random
    res <- matrix(nrow = h_iteration*length(local_communities), ncol = 2, NA)
    w <- 1
    time <- Sys.time()
    for(h in 1:h_iteration){
      local_communities <- metadata$sample_ID
      random_shuffle <- sample(length(local_communities))
      local_communities <- local_communities[random_shuffle]
      observed_MPDs <- observed_MPDs[random_shuffle,]
      for(i in 1:(length(local_communities)-5)){
        h_index <- rep(NA, h_iteration)
        p <- 1
        sum_quantiles <- rep(NA, length(local_communities))
        species_presence <- rowSums(table_presence_absence[,colnames(table_presence_absence) %in% local_communities])
        species_presence[species_presence != 0] <- 1
        for(m in local_communities){
          MPDs <- rep(NA, repetitions)
          rich_m <- richnes[m]
          ncomps <- rich_m * (rich_m-1)
          prob_vector <- species_presence
          MPDs <- sample_and_mpd(repetitions, h_iteration, trait, l_pool, rich_m, prob_vector, MPDs, ncomps)
          sum_quantiles[p] <- unname(stats::quantile(MPDs, probs = 0.05))
          p <- p + 1
        }
        real_MPDs_vector <- real_MPDs[real_MPDs$V1 %in% local_communities,]$V2
        comparision_vector <- real_MPDs_vector <= sum_quantiles
        res[w,1] <- metadata$env_variable[i]
        res[w,2] <- sum(comparision_vector)/length(local_communities)
        local_communities <- local_communities[-1]
        observed_MPDs <- observed_MPDs[-1,]
        w <- w + 1
      }
    }
    res.1 <- as.data.frame(res)
    a <- stats::aggregate(res.1$V2, by = list(res.1$V1), function(x) stats::quantile(x, probs = .025))
    b <- stats::aggregate(res.1$V2, by = list(res.1$V1), function(x) stats::quantile(x, probs = .5))
    c <- stats::aggregate(res.1$V2, by = list(res.1$V1), function(x) stats::quantile(x, probs = .975))
    d <- merge(a, b, by = 1)
    d <- merge(d, c, by = 1)
    d <- d[stats::complete.cases(d),]
    colnames(d) <- c("env_variable_value", "percentile_0.025", "percentile_0.5", "percentile_0.975")
    return(d)

  } else if ( model == 2) { # Model II : abundance

    res <- matrix(nrow = h_iteration*length(local_communities), ncol = 2, NA)
    w <- 1
    time <- Sys.time()
    for(h in 1:h_iteration){
      local_communities <- metadata$sample_ID
      random_shuffle <- sample(length(local_communities))
      local_communities <- local_communities[random_shuffle]
      observed_MPDs <- observed_MPDs[random_shuffle,]
      for(i in 1:(length(local_communities)-5)){
        h_index <- rep(NA, h_iteration)
        p <- 1
        sum_quantiles <- rep(NA, length(local_communities))
        species_presence <- rowSums(table_presence_absence[,colnames(table_presence_absence) %in% local_communities])
        species_presence[species_presence != 0] <- 1
        for(m in local_communities){
          MPDs <- rep(NA, repetitions)
          rich_m <- richnes[m]
          ncomps <- rich_m * (rich_m-1)
          prob_vector <- species_presence * species_abundances
          MPDs <- sample_and_mpd(repetitions, h_iteration, trait, l_pool, rich_m, prob_vector, MPDs, ncomps)
          sum_quantiles[p] <- unname(stats::quantile(MPDs, probs = 0.05))
          p <- p + 1
        }
        real_MPDs_vector <- real_MPDs[real_MPDs$V1 %in% local_communities,]$V2
        comparision_vector <- real_MPDs_vector <= sum_quantiles
        res[w,1] <- metadata$env_variable[i]
        res[w,2] <- sum(comparision_vector)/length(local_communities)
        local_communities <- local_communities[-1]
        observed_MPDs <- observed_MPDs[-1,]
        w <- w + 1
      }
    }
    res.1 <- as.data.frame(res)
    a <- stats::aggregate(res.1$V2, by = list(res.1$V1), function(x) stats::quantile(x, probs = .025))
    b <- stats::aggregate(res.1$V2, by = list(res.1$V1), function(x) stats::quantile(x, probs = .5))
    c <- stats::aggregate(res.1$V2, by = list(res.1$V1), function(x) stats::quantile(x, probs = .975))
    d <- merge(a, b, by = 1)
    d <- merge(d, c, by = 1)
    d <- d[stats::complete.cases(d),]
    colnames(d) <- c("env_variable_value", "percentile_0.025", "percentile_0.5", "percentile_0.975")
    return(d)

  } else if ( model == 3) { # Model III : salinity range

    res <- matrix(nrow = h_iteration*length(local_communities), ncol = 2, NA)
    w <- 1
    time <- Sys.time()
    for(h in 1:h_iteration){
      local_communities <- metadata$sample_ID
      random_shuffle <- sample(length(local_communities))
      local_communities <- local_communities[random_shuffle]
      observed_MPDs <- observed_MPDs[random_shuffle,]
      for(i in 1:(length(local_communities)-5)){
        h_index <- rep(NA, h_iteration)
        p <- 1
        sum_quantiles <- rep(NA, length(local_communities))
        species_presence <- rowSums(table_presence_absence[,colnames(table_presence_absence) %in% local_communities])
        species_presence[species_presence != 0] <- 1
        for(m in local_communities){
          MPDs <- rep(NA, repetitions)
          rich_m <- richnes[m]
          ncomps <- rich_m * (rich_m-1)
          sample_salinity_vector <- rep(as.numeric(metadata[metadata$sample_ID == m, env_var_col]), l_pool)
          presence <- as.numeric(dataset$min <= sample_salinity_vector & dataset$max >= sample_salinity_vector)
          prob_vector <- species_presence * presence
          MPDs <- sample_and_mpd(repetitions, h_iteration, trait, l_pool, rich_m, prob_vector, MPDs, ncomps)
          sum_quantiles[p] <- unname(stats::quantile(MPDs, probs = 0.05))
          p <- p + 1
        }
        real_MPDs_vector <- real_MPDs[real_MPDs$V1 %in% local_communities,]$V2
        comparision_vector <- real_MPDs_vector <= sum_quantiles
        res[w,1] <- metadata$env_variable[i]
        res[w,2] <- sum(comparision_vector)/length(local_communities)
        local_communities <- local_communities[-1]
        observed_MPDs <- observed_MPDs[-1,]
        w <- w + 1
      }
    }
    res.1 <- as.data.frame(res)
    a <- stats::aggregate(res.1$V2, by = list(res.1$V1), function(x) stats::quantile(x, probs = .025))
    b <- stats::aggregate(res.1$V2, by = list(res.1$V1), function(x) stats::quantile(x, probs = .5))
    c <- stats::aggregate(res.1$V2, by = list(res.1$V1), function(x) stats::quantile(x, probs = .975))
    d <- merge(a, b, by = 1)
    d <- merge(d, c, by = 1)
    d <- d[stats::complete.cases(d),]
    colnames(d) <- c("env_variable_value", "percentile_0.025", "percentile_0.5", "percentile_0.975")
    d <- d[stats::complete.cases(d),]
    return(d)

  } else if ( model == 4) { # Model IV : abundance and salinity range

    res <- matrix(nrow = h_iteration*length(local_communities), ncol = 2, NA)
    w <- 1
    time <- Sys.time()
    for(h in 1:h_iteration){
      local_communities <- metadata$sample_ID
      random_shuffle <- sample(length(local_communities))
      local_communities <- local_communities[random_shuffle]
      observed_MPDs <- observed_MPDs[random_shuffle,]
      for(i in 1:(length(local_communities)-5)){
        h_index <- rep(NA, h_iteration)
        p <- 1
        sum_quantiles <- rep(NA, length(local_communities))
        species_presence <- rowSums(table_presence_absence[,colnames(table_presence_absence) %in% local_communities])
        species_presence[species_presence != 0] <- 1
        for(m in local_communities){
          MPDs <- rep(NA, repetitions)
          rich_m <- richnes[m]
          ncomps <- rich_m * (rich_m-1)
          sample_salinity_vector <- rep(as.numeric(metadata[metadata$sample_ID == m, env_var_col]), l_pool)
          presence <- as.numeric(dataset$min <= sample_salinity_vector & dataset$max >= sample_salinity_vector)
          prob_vector <- species_presence * presence * species_abundances
          MPDs <- sample_and_mpd(repetitions, h_iteration, trait, l_pool, rich_m, prob_vector, MPDs, ncomps)
          sum_quantiles[p] <- unname(stats::quantile(MPDs, probs = 0.05))
          p <- p + 1
        }
        real_MPDs_vector <- real_MPDs[real_MPDs$V1 %in% local_communities,]$V2
        comparision_vector <- real_MPDs_vector <= sum_quantiles
        res[w,1] <- metadata$env_variable[i]
        res[w,2] <- sum(comparision_vector)/length(local_communities)
        local_communities <- local_communities[-1]
        observed_MPDs <- observed_MPDs[-1,]
        w <- w + 1
      }
    }
    res.1 <- as.data.frame(res)
    a <- stats::aggregate(res.1$V2, by = list(res.1$V1), function(x) stats::quantile(x, probs = .025))
    b <- stats::aggregate(res.1$V2, by = list(res.1$V1), function(x) stats::quantile(x, probs = .5))
    c <- stats::aggregate(res.1$V2, by = list(res.1$V1), function(x) stats::quantile(x, probs = .975))
    d <- merge(a, b, by = 1)
    d <- merge(d, c, by = 1)
    d <- d[stats::complete.cases(d),]
    colnames(d) <- c("env_variable_value", "percentile_0.025", "percentile_0.5", "percentile_0.975")
    return(d)

  }
}
