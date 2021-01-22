## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include = FALSE---------------------------------------------------------
run_bg = function(expr) {
  args = c("--vanilla", "-q")
  expr_c = deparse(substitute(expr))
  expr_c = paste(expr_c, collapse = "\n")
  pid = sys::r_background(c(args, "-e", expr_c), std_out = TRUE, std_err = TRUE)
  return(pid)
}

## -----------------------------------------------------------------------------
library(RestRserve)
app = Application$new()

## -----------------------------------------------------------------------------
calc_fib = function(n) {
  if (n < 0L) stop("n should be >= 0")
  if (n == 0L) return(0L)
  if (n == 1L || n == 2L) return(1L)
  x = rep(1L, n)
  
  for (i in 3L:n) {
   x[[i]] = x[[i - 1]] + x[[i - 2]] 
  }
  
  return(x[[n]])
}

## ---- req_res-----------------------------------------------------------------
fib_handler = function(.req, .res) {
  n = as.integer(.req$parameters_query[["n"]])
  if (length(n) == 0L || is.na(n)) {
    raise(HTTPError$bad_request())
  }
  .res$set_body(as.character(calc_fib(n)))
  .res$set_content_type("text/plain")
}

## -----------------------------------------------------------------------------
app$add_get(path = "/fib", FUN = fib_handler)

## -----------------------------------------------------------------------------
request = Request$new(path = "/fib", parameters_query = list(n = 10))
response = app$process_request(request)

cat("Response status:", response$status)
cat("Response body:", response$body)

## -----------------------------------------------------------------------------
yaml_file = system.file("examples", "openapi", "openapi.yaml", package = "RestRserve")
app$add_openapi(path = "/openapi.yaml", file_path = yaml_file)
app$add_swagger_ui(path = "/doc", path_openapi = "/openapi.yaml", use_cdn = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  backend = BackendRserve$new()
#  backend$start(app, http_port = 8080)

## ----echo = FALSE, eval = FALSE-----------------------------------------------
#  pid = run_bg({
#    library(RestRserve)
#    calc_fib = function(n) {
#      if (n < 0L) stop("n should be >= 0")
#      if (n == 0L) return(0L)
#      if (n == 1L || n == 2L) return(1L)
#      x = rep(1L, n)
#      for (i in 3L:n) x[[i]] = x[[i - 1]] + x[[i - 2]]
#      x[[n]]
#    }
#  
#    fib_handler = function(.req, .res) {
#      n = as.integer(.req$parameters_query[["n"]])
#      if (length(n) == 0L || is.na(n)) {
#        raise(HTTPError$bad_request())
#      }
#      .res$set_body(as.character(calc_fib(n)))
#      .res$set_content_type("text/plain")
#    }
#  
#    app = RestRserve::Application$new()
#    app$add_get(path = "/fib", FUN = fib_handler)
#    yaml_file = system.file("examples", "openapi", "openapi.yaml", package = "RestRserve")
#    app$add_openapi(path = "/openapi.yaml", file_path = yaml_file)
#    app$add_swagger_ui(path = "/doc", path_openapi = "/openapi.yaml", use_cdn = TRUE)
#    backend = BackendRserve$new()
#    backend$start(app, http_port = 8080)
#  })
#  tools::pskill(pid)

