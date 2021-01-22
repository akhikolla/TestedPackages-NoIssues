## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  results = 'markup'
)

## -----------------------------------------------------------------------------
library(RestRserve)
backend = BackendRserve$new()
application = Application$new()
application$add_get(path = "/factorial", function(.req, .res) {
  x = .req$get_param_query("x")
  x = as.integer(x)
  .res$set_body(factorial(x))
})

## -----------------------------------------------------------------------------
request = Request$new(
  path = "/factorial", 
  method = "GET", 
  parameters_query = list(x = 10)
)
response = application$process_request(request)
response

## -----------------------------------------------------------------------------
str(response$body)

## -----------------------------------------------------------------------------
response$content_type

## -----------------------------------------------------------------------------
response$encode

## -----------------------------------------------------------------------------
encode_decode_middleware = EncodeDecodeMiddleware$new()
app1  = Application$new(middleware = list())
app1$append_middleware(encode_decode_middleware)

app2 = Application$new()

## -----------------------------------------------------------------------------

FUN = encode_decode_middleware$ContentHandlers$get_encode('application/json')
FUN

## -----------------------------------------------------------------------------
application$add_get(path = "/factorial-json", function(.req, .res) {
  x = as.integer(.req$get_param_query("x"))
  result = factorial(x)
  .res$set_body(list(result = result))
  .res$set_content_type("application/json")
})

## -----------------------------------------------------------------------------
request = Request$new(
  path = "/factorial-json", 
  method = "GET", 
  parameters_query = list(x = 10)
)
response = application$process_request(request)

## -----------------------------------------------------------------------------
response$body

## -----------------------------------------------------------------------------
application$add_get(path = "/factorial-rds", function(.req, .res) {
  x = as.integer(.req$get_param_query("x"))
  result = factorial(x)
  body_rds = serialize(list(result = result), connection = NULL)
  .res$set_body(body_rds)
  .res$set_content_type("application/x-rds")
})


## -----------------------------------------------------------------------------
request = Request$new(
  path = "/factorial-rds", 
  method = "GET", 
  parameters_query = list(x = 10)
)
response = application$process_request(request)
response$body

## -----------------------------------------------------------------------------
application$add_get(path = "/factorial-rds2", function(.req, .res) {
  x = as.integer(.req$get_param_query("x"))
  result = factorial(x)
  body_rds = serialize(list(result = result), connection = NULL)
  .res$set_body(body_rds)
  .res$set_content_type("application/x-rds")
  .res$encode = identity
})

## -----------------------------------------------------------------------------
request = Request$new(
  path = "/factorial-rds2", 
  method = "GET", 
  parameters_query = list(x = 10)
)
response = application$process_request(request)
unserialize(response$body)

## -----------------------------------------------------------------------------
application = Application$new(content_type = "application/json")
application$add_post("/echo", function(.req, .res) {
  .res$set_body(.req$body)
})

request = Request$new(path = "/echo", method = "POST", body = '{"hello":"world"}', content_type = "application/json")
response = application$process_request(request)
response$body

## -----------------------------------------------------------------------------
encode_decode_middleware = EncodeDecodeMiddleware$new()

encode_decode_middleware$ContentHandlers$set_encode(
  "text/csv", 
  function(x) {
    con = rawConnection(raw(0), "w")
    on.exit(close(con))
    write.csv(x, con, row.names = FALSE)
    rawConnectionValue(con)
  }
)

encode_decode_middleware$ContentHandlers$set_decode(
  "text/csv", 
  function(x) {
    res = try({
      con = textConnection(rawToChar(x), open = "r")
      on.exit(close(con))
      read.csv(con)
    }, silent = TRUE)
    
    if (inherits(res, "try-error")) {
      raise(HTTPError$bad_request(body = attributes(res)$condition$message))
    }
    return(res)
  }
)

## -----------------------------------------------------------------------------
data(iris)
app = Application$new(middleware = list(encode_decode_middleware))

## -----------------------------------------------------------------------------
app$add_get("/iris", FUN = function(.req, .res) {
  .res$set_content_type("text/csv")
  .res$set_body(iris)
})

req = Request$new(path = "/iris", method = "GET")
res = app$process_request(req)

iris_out = read.csv(textConnection(rawToChar(res$body)))
head(iris_out)

## -----------------------------------------------------------------------------
app$add_post("/in", FUN = function(.req, .res) {
  str(.req$body)
})
req = Request$new(path = "/in", method = "POST", body = res$body, content_type = "text/csv")
app$process_request(req)

