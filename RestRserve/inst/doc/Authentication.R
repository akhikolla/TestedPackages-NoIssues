## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  results= 'markup'
)

## -----------------------------------------------------------------------------
library(RestRserve)

allowed_access = list(
  "user-1" = "password-1",
  "user-2" = "password-2"
)

auth_fun = function(user, password) {
  res = FALSE
  try({
    res = identical(allowed_access[[user]], password)
  }, silent = TRUE)
  return(res)
}

## -----------------------------------------------------------------------------
basic_auth_backend = AuthBackendBasic$new(FUN = auth_fun)

## -----------------------------------------------------------------------------

auth_mw = AuthMiddleware$new(
  auth_backend = basic_auth_backend, 
  routes = "/secure/factorial", 
  id = "auth_middleware"
)

app = Application$new(middleware = list(auth_mw))

## -----------------------------------------------------------------------------
factorial_handler = function(.req, .res) {
  x = .req$get_param_query("x")
  x = as.integer(x)
  .res$set_body(factorial(x))
}
app$add_get("/factorial", factorial_handler)
app$add_get("/secure/factorial", factorial_handler)

## -----------------------------------------------------------------------------
req = Request$new(path = "/factorial", parameters_query = list(x = "5"))
res = app$process_request(req)
res$body

## -----------------------------------------------------------------------------
req = Request$new(path = "/secure/factorial", parameters_query = list(x = "5"))
res = app$process_request(req)
res$body

## -----------------------------------------------------------------------------
credentials = jsonlite::base64_enc("user-1:password-1")
headers = list("Authorization" = sprintf("Basic %s", credentials))

req = Request$new(
  path = "/secure/factorial", 
  parameters_query = list(x = "5"), 
  headers = headers
)

res = app$process_request(req)
res$body

## -----------------------------------------------------------------------------
credentials = jsonlite::base64_enc("user-1:password-2")
headers = list("Authorization" = sprintf("Basic %s", credentials))

req = Request$new(
  path = "/secure/factorial", 
  parameters_query = list(x = "5"), 
  headers = headers
)

res = app$process_request(req)
res$body

## -----------------------------------------------------------------------------


allowed_tokens = c(
  "super_secure_token_1",
  "super_secure_token_2"
)

auth_fun = function(token) {
  res = FALSE
  try({
    res = token %in% allowed_tokens
  }, silent = TRUE)
  return(res)
}
basic_auth_backend = AuthBackendBearer$new(FUN = auth_fun)


## -----------------------------------------------------------------------------
auth_mw = AuthMiddleware$new(
  auth_backend = basic_auth_backend, 
  routes = "/secure/", 
  match = "partial",
  id = "auth_middleware"
)
app = Application$new(middleware = list(auth_mw))

## -----------------------------------------------------------------------------
app$add_get("/hello0", function(req, res) {res$body = "OK"})
app$add_get("/secure/hello1", function(req, res) {res$body = "OK"})
app$add_get("/secure/hello2", function(req, res) {res$body = "OK"})

## -----------------------------------------------------------------------------
headers = list("Authorization" = "Bearer super_secure_token_1")
req = Request$new(
  path = "/secure/hello1", 
  headers = headers
)

res = app$process_request(req)
res$body

## -----------------------------------------------------------------------------
headers = list("Authorization" = "Bearer abcd")
req = Request$new(
  path = "/secure/hello2", 
  headers = headers
)

res = app$process_request(req)
res$body

## -----------------------------------------------------------------------------
req = Request$new(path = "/hello0")
res = app$process_request(req)
res$body

