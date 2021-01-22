plot.eBsc <-
function(x,full = FALSE,...)
{
    etq   = x$etq.hat
    f.hat = x$f.hat
    R.hat = x$R.hat
    q.hat = x$q.hat
    data  = x$data
    out   = x$out

    if(q.hat != 0){                
        if(full){
            oldpar <- par(no.readonly = TRUE)
            on.exit(par(oldpar))

            par(mfrow = c(2,2), mar = c(2,2,3,2))
            plot(data, type='l', col=8, main = "data (grey), estimation (red)", ylab = "y", xlab = "x");
            lines(f.hat, col = 2, lwd = 1);
            plot(1:6, etq, col = 8, xlab = "", ylab = "", type = 'b', pch = 1, main = "Tq");
            abline(h = 0, lwd = 2);
            abline(v = q.hat, lwd = 1, col = 2, lty = 2)
            plot(data - f.hat, xlab = "", ylab = "", type = 'l', col = 8, main = "residuals");
            abline(h = 0, lwd = 2)
            acf(data - f.hat, lag.max = 20, xlab = "", ylab = "", main = "error's acf (black), acf.hat (red)")
            lines(seq(0, 20), R.hat[1:21], col=2, lwd=2)
            par(mfrow=c(1,1))
        }else{
            plot(data, type = 'l', col = 8, main = "data (grey), estimation (red)", ylab = "y", xlab = "x");
            lines(f.hat, col = 2, lwd = 1);
        }
    }
    
    if(q.hat==0){            
        plot(data, type = 'l', col = 8, main = "data (grey), estimations (colors)", ylab = "y", xlab = "x");
        for(i in 1:6){lines(out[[i]]$f.hat,type="l",col=i)}                      
        legend("bottomright",legend=c("q=1","q=2","q=3","q=4","q=5","q=6"),col=1:6,lty=1,box.col=0)
    }        
}
