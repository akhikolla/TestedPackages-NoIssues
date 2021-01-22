## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 5)

library(dplyr)
library(ipaddress)

## -----------------------------------------------------------------------------
integer_to_ip(c(0, 3232235521, 4294967295))

## -----------------------------------------------------------------------------
tibble(
  start = ip_address(c("192.168.0.0", "2001:db8::")),
  end = ip_address(c("192.168.0.15", "2001:db8::ffff:ffff:ffff"))
) %>%
  mutate(network = common_network(start, end))

## -----------------------------------------------------------------------------
my_networks <- tibble(
  network = ip_network(c("192.168.0.0/16", "2001:db8::/32")),
  label = c("Private", "Documentation")
)

my_networks

## -----------------------------------------------------------------------------
my_addresses <- tibble(
  address = ip_address(c("192.168.100.1", "1.2.3.4", "2001:db8::8a2e:370:7334", "::1"))
)

## -----------------------------------------------------------------------------
my_addresses %>%
  mutate(in_network = is_within_any(address, my_networks$network))

## -----------------------------------------------------------------------------
my_addresses %>%
  fuzzyjoin::fuzzy_left_join(my_networks, c("address" = "network"), is_within)

## ----accept-reject------------------------------------------------------------
sample_public <- function(size) {
  result <- sample_ipv4(size)
  
  all_public <- FALSE
  while (!all_public) {
    public <- is_global(result)
    n_replace <- sum(!public)
    
    if (n_replace == 0) {
      all_public <- TRUE
    } else {
      result[!public] <- sample_ipv4(n_replace)
    }
  }
  
  result
}

## -----------------------------------------------------------------------------
tibble(address = sample_public(10)) %>%
  mutate(is_public = is_global(address))

