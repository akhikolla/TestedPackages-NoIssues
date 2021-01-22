summary.rundtw <- function(object, ...){
   
   if(is.list(object$dist)){
      return( summary(do.call(c, object$dist), ...) )
   }else{
      return( summary(object$dist, ...) )
   }
   
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


summary.dba <- function(object, ...){
   do.call(rbind, lapply(object$iterDist_m2lot_norm, function(y) {summary(y, ...)}))
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


summary.idtw <- function(object, ...){
   ret <- list()
   
   if(!is.null(object$gcm)){
      ret <- c(ret, list(gcm = summary(as.vector(object$gcm), ...)))
   }
   
   if(!is.null(object$dm)){
      ret <- c(ret, list(dm = table(as.vector(object$dm), ...)))
   }
   
   if(!is.null(object$wp)){
      ret <- c(ret, list(wp = table(object$wp, ...)))
   }
   
   if(!is.null(object$cm)){
      ret <- c(ret, list(cm = summary(as.vector(object$cm), ...)))
   }
   
   if(!is.null(object$diffM)){
      ret <- c(ret, list(diffM = summary(as.vector(object$diffM), ...)))
   }
   
   if(!is.null(object$diffp)){
      ret <- c(ret, list(diffp = summary(as.vector(object$diffp), ...)))
   }
   
   return(ret)
}

# summary_rundtw <- summary.rundtw
