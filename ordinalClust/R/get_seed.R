get_seed <- function() {
  return(sample.int(.Machine$integer.max, 1))
}