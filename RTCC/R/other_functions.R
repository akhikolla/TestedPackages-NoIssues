real_mpd <- function(local_communities, table_presence_absence, dataset, i){
  real_MPDs <- as.data.frame(matrix(nrow=length(colnames(table_presence_absence))-1, ncol = 2, NA))
  a <- 1
  for(p in local_communities){
    temp <- as.data.frame(table_presence_absence[,colnames(table_presence_absence) %in% c("tax", as.character(p))])
    temp <- temp[temp[,2] != 0,]
    real_MPDs[a,1] <- p

    real_MPDs[a,2] <- (2/(nrow(dataset[dataset$tax %in% temp$tax,])*(nrow(dataset[dataset$tax %in% temp$tax,])-1)))*dist1(dataset[dataset$tax %in% temp$tax,i])
    a <- a+1
  }
  return(real_MPDs)
}


sample_and_mpd <- function(repetitions, h_iteration, trait, l_pool, rich_m, prob_vector, MPDs, ncomps){
  for(n in 1:(repetitions*h_iteration)){
    sintetic_community <- trait[sample.int(l_pool, rich_m, replace = F, prob = prob_vector)]
    MPDs[n] <- 2/ncomps * dist1(sintetic_community)
  }
  return(MPDs)
}


