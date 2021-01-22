cr.rand.max.inner.prod <-
function(XX,Ts,C_i,epp,M = 0,Plot = FALSE,cstar=0.95) {
  #mult.fip <- function (x.in,n,d,epp,C_i,Ts,min.draw,M,cstar) {
  #    out=multi_across_fip(X=x.in, M=M, min_draw=min.draw, tau=C_i,p=cstar, epp=epp,Ts = Ts)
  #    return(out)
  #}
  
 # med <- function(x){
 #     y<-stats:: quantile(x, 0.5, type = 3)[[1]]
 #     return(y[[1]])
  #  }
  
  
  Output = list(NULL)
  if (M == 0) M = floor(dim(XX)[1]*1)
  l=dim(XX)[2]
  n=dim(XX)[1]
  i=0
  min.draw=floor(log(Ts)^2/3)
  cstar=c(0.95,cstar)
  loc=multi_across_fip(X=XX, M=M, min_draw=min.draw, tau=tau.fun(XX[,1]), p=cstar,epp=epp,Ts= log(Ts))
   output=matrix(c(loc[[1]],loc[[2]],loc[[3]],loc[[4]]),M,4)
   max.b = output[which(abs(output[,2]) == max(abs(output[,2])))]
  max.inner = quantile(max.b, 0.5, type = 3)[[1]]#med(max.b)
  value.max.inner = sign(output[,2][abs(output[,2])==max(abs(output[,2]))])* max(abs(output[,2]))
  which.max = which(output[,2]==value.max.inner[1])
  if (Plot == TRUE) plot(x=output[,1],y=abs(output[,2]),xlab="index",ylab="inner products")    
  aux.mat=aggregate(abs(output[,2]),by=list(output[,1]),max)
  max.in.prod.series=matrix(0,dim(XX)[1],2)
  max.in.prod.series[,1] = 1:dim(XX)[1]
  
  for (i in 1:dim(XX)[1]) {
    if (sum(max.in.prod.series[i,1]==aux.mat[,1])==0) {
      max.in.prod.series[i,2]=0
    }  else {
      max.in.prod.series[i,2]=aux.mat[which(aux.mat[,1]==i),2]
    }
  }
  
  Output[[1]] = max.inner
  Output[[2]] = unique(value.max.inner)
  Output[[3]] = quantile(output[which.max,3], 0.5, type = 3)[[1]] #med(output[which.max,3])
  Output[[4]] = quantile(output[which.max,4], 0.5, type = 3)[[1]] #med(output[which.max,4])
  Output[[5]] = max.in.prod.series[,2]
  return(Output)
}
