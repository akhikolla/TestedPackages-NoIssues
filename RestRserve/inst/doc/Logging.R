## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  results = 'markup'
)

## -----------------------------------------------------------------------------
library(RestRserve)
application = Application$new()
application$logger$info("hello from logger")

## -----------------------------------------------------------------------------
application$logger$set_log_level("error")
application$logger$info("you won't see this message")

## -----------------------------------------------------------------------------
application$logger = Logger$new(level = "trace", name = "mylogger")

## -----------------------------------------------------------------------------
application = Application$new()
application$add_get("/sqrt", function(.req, .res) {
  .res$set_body(sqrt(x))
})

## -----------------------------------------------------------------------------
# let's emulate query string "/sqrt?x=10"
request = Request$new(path = "/sqrt", method = "GET", parameters_query = list(x = "10")) 
response = application$process_request(request)
response$body

## -----------------------------------------------------------------------------
fun2 = function(x) {
  sqrt(x)
}
fun1 = function(x) {
  fun2(x)
}
try(fun1('a'))

## -----------------------------------------------------------------------------
application$add_get("/sqrt2", function(.req, .res) {
  .res$set_body(fun1(x))
})
request = Request$new(path = "/sqrt2", method = "GET", parameters_query = list(x = "10")) 
response = application$process_request(request)

