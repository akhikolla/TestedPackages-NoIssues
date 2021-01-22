simulate_timewarp <- function (x, stretch = 0, compress = 0,
                               stretch_method = insert_linear_interp,
                               p_index = "rnorm", p_number = "rlnorm",
                               p_index_list = NULL, p_number_list = NULL, 
                               preserve_length = FALSE, seed = NULL, ...) 
{
   
   # stretch ...  numeric >= 0, 
   #                 e.g. if compress = 0 and ...
   #                    stretch = 0 then length(x_new) = length(x), or
   #                    stretch = 0.1 then length(x_new) = length(x) * 1.1, or
   #                    stretch = 1 then length(x_new) = length(x) * 2
   # compress ... numeric >= 0 and < 1, 
   #                 e.g. if stretch = 0 and ...
   #                    compress = 0 then length(x_new) = length(x), or
   #                    compress = 0.2 then length(x_new) = length(x) * 0.8
   
   
   if(stretch < 0) stop("stretch needs to be >= 0")
   if(compress < 0 | compress >= 1) stop("compress needs to be >= 0 and < 1")
   if(compress == 0 & stretch == 0) return(x)
   if(!is.null(seed)) set.seed(seed)
   
   
   if(is.vector(x)){
      return(simulate_timewarp_vec(x, stretch, compress, stretch_method, 
                                   p_index, p_number, p_index_list, p_number_list, 
                                   preserve_length = preserve_length, ...))
   }else if (is.matrix(x)){
      return(simulate_timewarp_mat(x, stretch, compress, stretch_method, 
                                   p_index, p_number, p_index_list, p_number_list,
                                   preserve_length = preserve_length, ...))
   }else{
      stop("x needs to be a vector or matrix")
   }
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


simulate_timewarp_vec <- function (x, stretch = 0, compress = 0,
                                   stretch_method = insert_linear_interp,  
                                   p_index, p_number, p_index_list = NULL, 
                                   p_number_list = NULL, preserve_length = FALSE,...)
{
   
   x_new <- x
   if(preserve_length) {
      if(stretch <= 0) stretch <- compress
      compress <- stretch
   }
   
   
   if (stretch > 0) {
      nx <- length(x_new)
      nx_new <- ceiling(nx * (1 + stretch))
      
      rno2i <- nx_new - nx #remaining number of observations to insert
      while(rno2i > 0){
         
         #draw number of obs to insert
         no2i <- trunc_sample(1, a = 1, b = rno2i, prob = p_number, plist = p_number_list)
         #draw index where to insert obs
         ix   <- trunc_sample(1, a = 2, b = (nx - 1), prob = p_index , plist = p_index_list)
         
         nx_new <- length(x_new)
         x_new <- stretch_method(x = x_new,ix = ix, N = no2i, ...)
         rno2i <- rno2i - no2i
      }
   }
   
   if (compress > 0) {
      if(preserve_length) {
         rno2o <- length(x_new) - length(x)
      }else{
         nx <- length(x_new)
         nx_new <- ceiling(nx * (1 - compress))
         rno2o <- nx - nx_new # remaining number of observations to omit
      }
      
      while(rno2o > 0){
         nx_new <- length(x_new)
         
         # draw number of obs to omit
         no2o <- trunc_sample(1, a = 1, b = rno2o, prob = p_number, plist = p_number_list)
         # draw index where to omit obs
         ix   <- trunc_sample(1, a = 1, b = (nx_new - no2o + 1), prob = p_index , plist = p_index_list)
         
         x_new <- x_new[-c(ix:(ix + no2o - 1))] # omit obs
         rno2o <- rno2o - no2o
      }
      
   }
   return(x_new)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


simulate_timewarp_mat <- function (x, stretch = 0, compress = 0,
                                   stretch_method = insert_linear_interp,  
                                   p_index, p_number, p_index_list = NULL, 
                                   p_number_list = NULL, preserve_length = FALSE,...)
{
   
   x_new <- x
   if(preserve_length) {
      if(stretch <= 0) stretch <- compress
      compress <- stretch
   }
   
   if (stretch > 0) {
      nx <- nrow(x_new)
      nx_new <- ceiling(nx * (1 + stretch))
      
      rno2i <- nx_new - nx #remaining number of observations to insert
      while(rno2i > 0){
         
         #draw number of obs to insert
         no2i <- trunc_sample(1, a = 1, b = rno2i, prob = p_number, plist = p_number_list)
         #draw index where to insert obs
         ix   <- trunc_sample(1, a = 2, b = (nx - 1), prob = p_index , plist = p_index_list)
         
         nx_new <- nrow(x_new)
         x_new <- do.call(cbind, lapply(1:ncol(x), function(i){
            stretch_method(x = x_new[,i], ix = ix, N = no2i, ...)
         }))
         
         rno2i <- rno2i - no2i
      }
   }
   
   if (compress > 0) {
      if(preserve_length) {
         rno2o <- nrow(x_new) - nrow(x)
      }else{
         nx <- nrow(x_new)
         nx_new <- ceiling(nx * (1 - compress))
         rno2o <- nx - nx_new # remaining number of observations to omit
      }
      
      while(rno2o > 0){
         nx_new <- nrow(x_new)
         
         # draw number of obs to omit
         no2o <- trunc_sample(1, a = 1, b = rno2o, prob = p_number, plist = p_number_list)
         # draw index where to omit obs
         ix   <- trunc_sample(1, a = 1, b = (nx_new - no2o + 1), prob = p_index , plist = p_index_list)
         
         x_new <- x_new[-c(ix:(ix + no2o - 1)), , drop = FALSE] # omit obs
         rno2o <- rno2o - no2o
      }
      
   }
   return(x_new)
}



#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

insert_const <- function(x, ix, N, const = NULL){
   if(is.null(const)){
      const <- x[ix]
   }
   c(x[1:ix],  rep(const, N),  x[(ix + 1):length(x)])
}


insert_linear_interp <- function(x, ix, N){
   c(x[1:ix],  
     approx(c(x[ix], x[ix + 1]), n = N)$y,
     x[(ix + 1):length(x)])
}



insert_norm <- function(x, ix, N, mean = 0, sd = 1){
   c(x[1:ix], 
     rep(x[ix], N) + rnorm(N, mean = mean, sd = sd),
     x[(ix + 1):length(x)])
}


insert_linear_norm <- function(x, ix, N, mean = 0, sd = 1){
   c(x[1:ix],  
     approx(c(x[ix], x[ix + 1]), n = N)$y + rnorm(N, mean = mean, sd = sd),
     x[(ix + 1):length(x)])
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


trunc_sample <- function (N, a, b, prob, plist = NULL) {
   
   if(is.null(plist)){
      plist_str <- ")"
   }else{
      plist_str <- paste0( ",", 
                           paste(paste(names(plist), plist, sep = " = "), collapse = ", "),
                           ")")
   }
   tmp_txt <- paste0("y <- ", prob, "(n = ", N, plist_str)
   eval(parse(text = tmp_txt))
   
   # quantile function
   qmin <- 0
   qmax <- 0
   qdist <- paste0("q", substr(prob, start = 2, stop = nchar(prob)))
   tmp_txt <- paste0("qmax <- ", qdist, "(p = 0.999", plist_str)
   eval(parse(text = tmp_txt))
   tmp_txt <- paste0("qmin <- ", qdist, "(p = 0.001", plist_str)
   eval(parse(text = tmp_txt))
   
   if(qmax == qmin){
      return(rep(qmax, N))
   }
   
   y[y > qmax] <- qmax
   y[y < qmin] <- qmin
   
   # min-max scaling and
   # shift to appreciated interval
   return(round((y - qmin )/( qmax - qmin) * (b-a) + a))
}


# trunc_sample(N = 10^3, a=1, b = 30, prob = "runif") %>% qplot
# trunc_sample(N = 10^3, a=1, b = 30, prob = "rnorm") %>% qplot
# trunc_sample(N = 10^3, a=1, b = 30, prob = "rnorm", plist = list(sd = 10)) %>% qplot
