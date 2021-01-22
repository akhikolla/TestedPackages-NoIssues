summary.eBsc <-
function(object,...)
{
    q.hat = object$q.hat
    lambda.hat = object$lambda.hat
    res = object$res
    etq.hat = object$etq.hat    

    if(length(etq.hat)==1){
        cat("Computations performed according to the provided q = ",q.hat,sep="","\n")
        cat("lambda.hat = ",toString(lambda.hat),sep="","\n")
    }else{
        cat("etq.hat = ",etq.hat,sep="","\n")
        cat("q.hat = ",q.hat,sep="","\n")
        cat("lambda.hat = ",toString(lambda.hat),sep="","\n")
    }
    
}
