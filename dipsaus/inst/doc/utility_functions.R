## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(dipsaus)

## ---- eval=FALSE--------------------------------------------------------------
#  cat2('Debug passed!', level = 'DEBUG')
#  #> ✔ Debug passed!
#  
#  cat2('You are all set.', level = 'INFO')
#  #> ♥ You are all set.
#  
#  cat2('Wait a second...', level = 'WARNING')
#  #> ⚠ Wait a second...
#  
#  cat2('Ooops', level = 'ERROR')
#  #> ✖ Ooops
#  
#  cat2('Bi--Doop---', level = 'FATAL')
#  #> ✖ Bi--Doop---
#  #> Error:
#  #> ...

## -----------------------------------------------------------------------------
parse_svec("7-10,14-15")

## -----------------------------------------------------------------------------
deparse_svec(c(2,5,3,1,7))

## -----------------------------------------------------------------------------
deparse_svec(c(1,2,4,7,11))

deparse_svec(c(1,2,4,7,11), max_lag = 2)

deparse_svec(c(1,2,4,7,11), max_lag = 3)

## -----------------------------------------------------------------------------
# Most of the CPU Chipset and vendor
get_cpu()

# Total RAM in bytes
get_ram()

# Print-friendly
to_ram_size(get_ram(), 1024)

# WARNING: $free is the total RAM - R usage, is no the actual free RAM
mem_limit2()

## ---- eval=FALSE--------------------------------------------------------------
#  > ask_yesno('Please answer an yes/no question, ok?')
#  ## ♥ Please answer an yes/no question, ok? (Yes/no):
#  > qweee
#  ## ⚠ Please answer Y/yes, N/no, or c to cancel. (Yes/no):
#  > ttt
#  ## ⚠ Please answer Y/yes, N/no, or c to cancel. (Yes/no):
#  > y
#  ## [1] TRUE

## ---- eval=FALSE--------------------------------------------------------------
#  > ask_or_default("What is your password", default = 'I will not tell you!')
#  ## ♥ What is your password
#  ##   [default is ‘I will not tell you!’]
#  >
#  ## [1] "I will not tell you!"

