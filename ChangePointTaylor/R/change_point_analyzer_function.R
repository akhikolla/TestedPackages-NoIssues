#' @useDynLib ChangePointTaylor
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr "%>%"
#' @importFrom stats "complete.cases"
#' @importFrom stats "quantile"
#' @importFrom rlang ":="
#' @importFrom rlang "!!"
#' @importFrom rlang ".data"
NULL

if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

bootstrap_S_diff <- function(x, n_smpls = 1000){

  X_bar <- mean.default(x)
  length_x <- length(x)

  # 1. Generate a bootstrap sample of n units, denoted X01, X02, ..., Xn, by randomly reordering the original n values.
  x_btstrps <- purrr::map(1:n_smpls, function(i) sample_cpp(x,length_x))

  #2 Based on the bootstrap sample, calculate the bootstrap CUSUM, denoted S00, S01, ..., Sn.
  S_btstrps <- purrr::map(x_btstrps, cusum, X_bar)
  #Calculate the maximum, minimum and difference of the bootstrap CUSUM, denoted S0max, S0min, and S0diff
  S_diff_btstrp <- purrr::map_dbl(S_btstrps, S_diff)

  return(S_diff_btstrp)
}


confidence_change_occurred <- function(x, n_bootstraps = 1000){
  #1. First calculate the average
  X_bar <- mean.default(x)

  S_x <- cusum(x, X_bar)
  # plot(S_x)
  S_diff_x <- S_diff(S_x)
  # print(S_diff_x)

  S_diff_btstrp <- bootstrap_S_diff(x, n_bootstraps)
  return(sum(S_diff_btstrp < S_diff_x)/n_bootstraps)
}





# Once a change has been detected, an estimate of when the change occurred can be made. One such estimator is the CUSUM estimator

# Sm is the point furthest from zero in the CUSUM chart. The point m estimates last point before the change occurred. The point m+1 estimates the first point after the change.
Sm_ix <- function(x){
  S <- cusum(x, mean.default(x))[2:(length(x) + 1)] #remove S_0 at ix 1
  which(abs(S) == max(abs(S)))
}




get_change <- function(x, min_conf = 0.5, recursive = T, mthd, level = 1, recursive_call = F, original_indices = NA, n_bootraps = 1000){
  if(is.na(original_indices[1]) & !recursive_call){
    # message("no indices passed to function. Assuming all original values were provided...")
    original_x_indices <- 1:length(x)
  }else{
    original_x_indices <- original_indices
  }
  minimum_confidence <- min_conf
  change_conf <- confidence_change_occurred(x,n_bootraps)
  # print(x)
  # print(change_conf)
  if(change_conf < min_conf | length(x) <5){#minimum segment length is 5
    return(data.frame('change_ix' = NA, 'change_conf' = NA, 'level' = NA))
  }
  if(mthd == "MSE"){
    change_ix <- min_MSE_ix(x)

  }else if(mthd == "CUSUM"){
    change_ix <- Sm_ix(x)

  }else{
    stop("method must be 'MSE' or 'CUSUM'")
  }

  original_x_change_ix <- original_x_indices[change_ix]

  changes <- data.frame('change_ix' = original_x_change_ix, 'change_conf' = change_conf, 'level' = level)

  if(recursive == F){

    return(changes)
  }else{
    x_indices_before_change <- original_x_indices[1:(change_ix)]
    x_indices_after_change <- original_x_indices[(change_ix+1):length(x)]

    x_before_change <- x[1:(change_ix)]
    x_after_change <- x[(change_ix+1):length(x)]
    changes_before <- get_change(x_before_change, original_indices = x_indices_before_change, min_conf = minimum_confidence, mthd = mthd, recursive_call = T, level = level + 1,n_bootraps = n_bootraps)
    changes_after <- get_change(x_after_change, original_indices = x_indices_after_change, min_conf = minimum_confidence, mthd = mthd, recursive_call = T, level = level + 1,n_bootraps = n_bootraps)

    return(
      dplyr::bind_rows(changes, changes_before,changes_after) %>%
        dplyr::filter(complete.cases(.))
    )
  }
}


#make function that will speed up bootrapping CI
get_change_ix_only <- function(x,mthd, original_indices){
  if(mthd == "MSE"){
    change_ix <- min_MSE_ix(x)

  }else if(mthd == "CUSUM"){
    change_ix <- Sm_ix(x)

  }else{
    stop("method must be 'MSE' or 'CUSUM'")
  }
  return(original_indices[change_ix])
}



get_change_sections <- function(x, change_df, min_conf, mthd, recursive = F,filter_by,filter_values){
  sorted_change_df <- change_df %>%
    dplyr::arrange(.data$change_ix)

  change_ixs <- sorted_change_df$change_ix
  sorted_initial_change_ixs <- change_ixs
  change_sections <- list()
  change_section_start_stops <- c(0,sorted_initial_change_ixs, length(x))
  change_section_orig_ixs <- list()
  # print(change_section_start_stops)
  for(i in 1:length(change_ixs)){
    change_section_seq <- (change_section_start_stops[i]+1):(change_section_start_stops[i+2])
    change_sections[[i]] <- x[change_section_seq]
    change_section_orig_ixs[[i]] <- change_section_seq
  }

  if(filter_by == "levels"){
    reestimate_lgl <- sorted_change_df$level %in% c(filter_values)

  }else if(filter_by == "index"){
    reestimate_lgl <- sorted_change_df$change_ix %in% c(filter_values)

  }

  return(list(x = change_sections[reestimate_lgl], original_indices = change_section_orig_ixs[reestimate_lgl]
              , recursive = recursive, min_conf = min_conf,mthd = mthd
              , change_ix = sorted_initial_change_ixs[reestimate_lgl], level = sorted_change_df$level[reestimate_lgl]))
}


reestimate_changes_by_level <- function(x,df_changes_for_reestimation, min_conf, mthd, levels_to_reestimate){
  changes_in_non_reestimation_levels <- df_changes_for_reestimation %>%
    dplyr::filter(!(.data$level %in% c(levels_to_reestimate)))

  changes_in_reestimation_levels <- df_changes_for_reestimation %>%
    dplyr::filter(.data$level %in% c(levels_to_reestimate))

  changes_for_reestimation <- get_change_sections(x,df_changes_for_reestimation, min_conf = min_conf, mthd = mthd, filter_values = levels_to_reestimate, filter_by = "levels")

  # print(changes_for_reestimation[c('x','original_indices','recursive','min_conf','mthd')])

  reestimated_changes_df <- purrr::pmap_df(changes_for_reestimation[c('x','original_indices','recursive','min_conf','mthd')],get_change) %>%
    dplyr::mutate(level = changes_for_reestimation$level) %>%
    dplyr::bind_rows(changes_in_non_reestimation_levels)

  return(reestimated_changes_df)
}

reestimate_change_level_seq <- function(x,df_changes_for_reestimation, min_conf, mthd, level_sequence){
  for(level_to_reestimate in level_sequence){
    df_changes_for_reestimation <- reestimate_changes_by_level(x,df_changes_for_reestimation, min_conf, mthd, level_to_reestimate) %>%
      dplyr::filter(complete.cases(.)) #sometimes we lose one due to the min 5 length segment.
  }
  df_changes_for_reestimation %>%
    dplyr::group_by(.data$change_ix) %>%
    dplyr::summarize(change_conf = max(.data$change_conf)
                     ,level = min(.data$level)) %>%
    return()
}


drop_any_under_threshold_in_given_level <- function(change_df, drop_level,min_tbl_conf){
  change_df %>%
    dplyr::filter(!(.data$level == drop_level & .data$change_conf < min_tbl_conf)) %>%
    return()
}

drop_lowest_change_conf_in_given_level <- function(change_df, drop_level){
  min_change_conf_in_given_level <- min(change_df$change_conf[change_df$level == drop_level])
  change_df %>%
    dplyr::filter(!(.data$level == drop_level & .data$change_conf == min_change_conf_in_given_level)) %>%
    return()
}


get_all_changes <- function(x, mthd, labels = NA, n_bootraps = 1000, min_candidate_conf = 0.5, min_tbl_conf = 0.9){
  changes_for_reestimation_df <- get_change(x, min_conf = min_candidate_conf, mthd = mthd, n_bootraps = n_bootraps)
  labels_lookup_df <- data.frame('label' = labels, change_ix = 1:length(labels))

  # print(changes_for_reestimation_df %>%
  #         left_join(labels_lookup_df, by = c("change_ix"))%>%
  #         mutate(data = "Change Candidates"))

  n_change_rows <- nrow(changes_for_reestimation_df)

  if(n_change_rows == 1 & is.na(changes_for_reestimation_df$change_conf[1])){
    # print("0 Changes Identified After Re-Estimation")
    return(data.frame())
  }
  # print(paste0(n_change_rows, " Candidate Change(s) Identified"))
  #if there is only one change initially, check to see if its above our tbl threshold.
  if(n_change_rows == 1){
    if(changes_for_reestimation_df$change_conf >min_tbl_conf) {
      reestimated_changes_df <- changes_for_reestimation_df %>%
        dplyr::mutate(change_ix = .data$change_ix + 1)
      changes_df_w_labels <- labels_lookup_df %>%
        dplyr::right_join(reestimated_changes_df, by = c("change_ix"))
      # print("1 Change Identified After Re-Estimation")
      return(changes_df_w_labels)
    }else{
      # print("0 Changes Identified After Re-Estimation")
      return(data.frame())
    }
  }
  # print(changes_for_reestimation_df %>%
  #         left_join(labels_lookup_df, by = c("change_ix"))%>%
  #         mutate(data = "Change Candidates"))
  changes_for_reestimation_df <- drop_any_under_threshold_in_given_level(changes_for_reestimation_df, max(changes_for_reestimation_df$level),min_tbl_conf)

  all_changepoints_above_tbl_threshold <- FALSE
  while(!all_changepoints_above_tbl_threshold){

    # print(changes_for_reestimation_df %>%
    #         left_join(labels_lookup_df, by = c("change_ix"))%>%
    #         mutate(data = "While Loop Start"))

    if(nrow(changes_for_reestimation_df) == 0){
      # print("nrow = 0")
      changes_for_reestimation_df <- get_change(x, min_conf = min_tbl_conf, mthd = mthd, recursive = F, n_bootraps = n_bootraps) %>%
        # print() %>%
        dplyr::filter(complete.cases(.)) %>%
        dplyr::group_by(.data$change_ix) %>%
        dplyr::summarize(change_conf = max(.data$change_conf))

    }else{
      ### reestimate all change candidates top to bottom ###
      level_sequence_top_to_bottom <- max(changes_for_reestimation_df$level):min(changes_for_reestimation_df$level)
      changes_for_reestimation_df <- reestimate_change_level_seq(x,changes_for_reestimation_df, min_conf = 0, mthd = mthd, level_sequence_top_to_bottom)
    }


    all_changepoints_above_tbl_threshold <- sum(changes_for_reestimation_df$change_conf > min_tbl_conf) == nrow(changes_for_reestimation_df)

    #if all changes aren't above the change threshold, remove the change with the lowest confidence and re-estimate.
    if(!all_changepoints_above_tbl_threshold){

      # changes_for_reestimation_df <- drop_lowest_change_conf(changes_for_reestimation_df)

      #level based removal
      lowest_level_with_change_under_tbl_threshold <- max(changes_for_reestimation_df$level[changes_for_reestimation_df$change_conf < min_tbl_conf])
      # changes_for_reestimation_df <- drop_any_under_threshold_in_given_level(changes_for_reestimation_df,lowest_level_with_change_under_tbl_threshold)
      changes_for_reestimation_df <- drop_lowest_change_conf_in_given_level(changes_for_reestimation_df,lowest_level_with_change_under_tbl_threshold)

    }else{
      reestimated_changes_df <- changes_for_reestimation_df %>%
        dplyr::mutate(change_ix = .data$change_ix + 1)
    }
  }

  message(paste0(nrow(reestimated_changes_df), " Change(s) Identified"))
  if(any(!is.na(labels))){
    changes_df_w_labels <- labels_lookup_df %>%
      dplyr::right_join(reestimated_changes_df, by = c("change_ix")) %>%
      dplyr::arrange(.data$change_ix)
    return(changes_df_w_labels)

  }else{
    message("NA supplied to 'label' argument")
    reestimated_changes_df %>%
      dplyr::mutate(label = NA) %>%
      dplyr::arrange(.data$change_ix) %>%
      return()
  }
}








percent <- function(x, digits = 0, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}









# From Dr. Taylor:
# -Again, just the two intervals each side of the change up to the closest changes are used.
#  A bootstrap is a random reordering of the data on each side of the change point.
#  The change-point is than reestimated. The 2.5% percentile and 97.5% percentile of the estimates are used to construct the 95% confidence interval.
get_bootstraps_CIs <- function(chng_df,x,mthd, labels = NA, CI_level = 0.95,n_boots = 1000){
  # print(chng_df)
  # print("estimating confidence intervals....")
  if(is.na(labels[1])){
    labels <- as.character(1:length(x))
  }
  change_ixs <- chng_df$change_ix
  # print(change_ixs)
  sorted_initial_change_ixs <- sort(change_ixs)
  change_section_start_stops <- c(1,sorted_initial_change_ixs, length(x)+1)
  change_section_orig_ixs <- list()
  for(i in 1:length(change_ixs)){
    change_section_orig_ixs[[i]] <- (change_section_start_stops[i]):(change_section_start_stops[i+2]-1)
  }

  mp_bootstrap_CI <- function(change_sec_orig_ixs,chg_ix, orig_x,CI, mthd) {
    tryCatch({
      # print(change_sec_orig_ixs)
      # print(chg_ix)
      #get ix of the change for this particular section based on the original value
      sec_change_ix <- which(change_sec_orig_ixs == chg_ix)
      # print(sec_change_ix)
      # print(mthd)
      pre_change_ixs <- change_sec_orig_ixs[1:sec_change_ix-1]
      post_change_ixs <- change_sec_orig_ixs[sec_change_ix:length(change_sec_orig_ixs)]

      bootstraped_change_ixs <- purrr::map_int(1:n_boots, function(i) get_change_ix_only(orig_x[c(sample(pre_change_ixs), sample(post_change_ixs))],mthd = mthd, original_indices = change_sec_orig_ixs))

      # bootstraped_change_ixs <- replace_na(bootstraped_change_ixs, chg_ix)

      CI_labels <- labels[quantile(bootstraped_change_ixs + 1, c((1-CI)/2, 1-((1-CI)/2)),na.rm = T)]

      #if change butts up againsts another change, the conf interval can return an NA. replace with original label value
      CI_labels <- tidyr::replace_na(CI_labels, labels[chg_ix])
      # print(CI_labels)
      # print(order(CI_labels))
      # print(order(CI_labels)[labels])
      # print(sort(order(CI_labels)[labels]))
      CI_labels <- CI_labels[order(match(CI_labels,labels))] #sort based on original
      # CI_labels <- CI_labels[sort(order(CI_labels)[labels])]
      CI_labels_str <- paste0("(",paste0( CI_labels, collapse = " - "),")")
      return(CI_labels_str)},
      error = function(e) conditionMessage(e))
  }

  CIs <- purrr::map2_chr(change_section_orig_ixs,sorted_initial_change_ixs,  mp_bootstrap_CI, orig_x = x, CI = CI_level, mthd = mthd)
  # print(CIs)
  chng_df %>%
    dplyr::mutate(!!as.name(paste0("CI (",percent(CI_level),")")) := CIs)
}



get_change_levels <- function(chng_df,x){
  change_ixs <- unique(chng_df$change_ix)
  sorted_initial_change_ixs <- sort(change_ixs)
  change_section_start_stops <- c(1,sorted_initial_change_ixs, length(x)+1)
  change_section_orig_ixs <- list()
  change_sections <- list()
  for(i in 1:(length(change_ixs)+1)){
    change_section_seq <- (change_section_start_stops[i]):(change_section_start_stops[i+1]-1)
    change_sections[[i]] <- x[change_section_seq]
  }
  # print(change_sections)
  section_means <- purrr::map_dbl(change_sections, mean.default)
  from <- NA
  to <- NA

  for(i in 1:length(change_ixs)){
    from[i] <- section_means[i]
    to[i] <- section_means[i+1]
  }

  chng_df %>%
    dplyr::arrange(.data$change_ix) %>%
    dplyr::mutate(From = from) %>%
    dplyr::mutate(To = to) %>%
    return()
}


#' change_point_analyzer
#' 
#' 
#' a simple implementation of the change in mean detection \href{https://variation.com/wp-content/uploads/change-point-analyzer/change-point-analysis-a-powerful-new-tool-for-detecting-changes.pdf}{methods} developed by Wayne Taylor and utilized in his \href{https://variation.com/product/change-point-analyzer/}{Change Point Analyzer} software. The package recursively uses the 'MSE' change point calculation to identify candidate change points. Taylor's backwards elimination process is then employed to come up with a final set of change points.
#'
#' @param x a numeric vector
#' @param labels a vector the same length as \code{x}. Will generate labels for the change points in the output dataframe. 
#' @param n_bootstraps an integer value. Determines the number of bootstraps when calculating the change confidence level. 
#' @param min_candidate_conf a value between 0 and 1. The minimum change confidence level to become a candidate change point before re-estimation and backwards elimination. 
#' @param min_tbl_conf a value between 0 and 1. The minimum change confidence level below which a candidate change point will be eliminated after re-estimation and backwards elimination. 
#' @param CI a value between 0 and 1. The value of the confidence interval.
#'
#' @return a dataframe containing the change points, their confidence levels, and other relevant information
#' @export
#'
#' @references \href{https://variation.com/wp-content/uploads/change-point-analyzer/change-point-analysis-a-powerful-new-tool-for-detecting-changes.pdf}{Taylor, W. A. (2000). Change-point analysis: a powerful new tool for detecting changes.}
#' 
#' @examples
#' x <- US_Trade_Deficit$deficit_billions
#' label_vals <- US_Trade_Deficit$date
#' 
#' change_point_analyzer(x)
#' 
#' change_point_analyzer(x, label = label_vals)
#' 
#' change_point_analyzer(x, label = label_vals, n_bootstraps = 10000)
#' 
#' change_point_analyzer(x, label = label_vals, min_candidate_conf = 0.66,  min_tbl_conf = 0.95)
change_point_analyzer <- function(x, labels = NA, n_bootstraps = 1000, min_candidate_conf = 0.5, min_tbl_conf = 0.9, CI = 0.95){
  if(!is.numeric(x) | length(x) < 5){
    stop("Invalid x argument. 'x' must be a numeric vector with length(x) >= 5")
  }
  
  if(!is.na(labels[1]) & length(x) != length(labels)){
    stop("Invalid labels argument. length(x) != length(labels)")
  }
  
  if(!is.numeric(n_bootstraps) | !dplyr::between(n_bootstraps, 100,1000000)){
    stop("Invalid n_bootraps argument. n_bootraps must be a numeric value between 100 and 1,000,000")
  }
  
  if(!is.numeric(min_candidate_conf) | !dplyr::between(min_candidate_conf, 0.3,1)){
    stop("Invalid min_candidate_conf argument. min_candidate_conf must be a numeric value between 0.3 and 1")
  }
  
  if(!is.numeric(min_tbl_conf) | !dplyr::between(min_tbl_conf, 0.5,1)){
    stop("Invalid min_tbl_conf argument. min_tbl_conf must be a numeric value between 0.5 and 1")
  }
  
  if(!is.numeric(min_tbl_conf) | !dplyr::between(min_tbl_conf, 0.9,0.999)){
    stop("Invalid CI argument. CI must be a numeric value between 0.9 and 0.999")
  }
  
  method <-  "MSE"
  
  tryCatch({
    all_changes_df <- get_all_changes(x, mthd = method, labels, min_candidate_conf = min_candidate_conf,min_tbl_conf = min_tbl_conf)
    #if there aren't any changes return that.
    if(nrow(all_changes_df)==0){
      data.frame(change_ix = NA, label = NA, CI_init = NA
                 , change_conf = NA, From = NA, To= NA  ) %>%
        dplyr::rename(!!as.name(paste0("CI (",percent(CI),")")) := .data$CI_init) %>%
        return()
    }else{
      if(is.na(CI)){
        all_changes_df_with_CI <- all_changes_df %>%
          dplyr::mutate(CI = NA)
      }else{
        all_changes_df_with_CI <- all_changes_df %>%
          get_bootstraps_CIs(x, mthd = method, labels,CI , n_boots = n_bootstraps)
      }
      all_changes_df_with_CI %>%
        get_change_levels(x) %>%
        dplyr::select(.data$change_ix, .data$label, dplyr::matches("CI"), .data$change_conf, .data$From, .data$To) %>%
        return()
    }
  },
  error = function(e) {
    warning(e)
    data.frame(change_ix = NA, label = NA, CI_init = "error"
               , change_conf = NA, From = NA, To= NA  ) %>%
      dplyr::rename(!!as.name(paste0("CI (",percent(CI),")")) := .data$CI_init) %>%
      return()
  }
  )

}





