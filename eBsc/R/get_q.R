 get_q <- function(etq){
    
  if(mean(etq>0)==1) q.hat=1 else if(mean(etq>0)==0) q.hat=6 else if((sum(diff(etq>0))==1)&(sum((diff(etq>0)==1))==1))  
  {
    etq.ind=which(etq>0)
    ind=etq.ind[1]; etq.min=min(abs(etq[ind]),abs(etq[ind-1])); q.hat=which(abs(etq)==etq.min)
  } else q.hat=0
    q.hat
  }
  
