
GenerateResidual=function(x,y,S,suffStat){
  mat=NULL
  if (length(S)!=0) mat=suffStat$dat_in[,S]
  x.= as.matrix(suffStat$dat_in[, x])
  y.= as.matrix(suffStat$dat_in[, y])
  
  if (suffStat$type[x] == "c"){
    
    if (suffStat$level[x]==2){
      if (is.null(mat)){
        fit=glm(x.~1,family = "binomial")
        pi_hat=predict(fit, as.data.frame(rep(1,nrow(suffStat$dat_in))), type="response")
      }else{
        fit=logistf::logistf(x.~mat)
        pi_hat=fit$predict
      }
      residual_x=(x.-pi_hat)/((pi_hat*(1-pi_hat))^0.5)
      
    }else{
        if (is.null(mat)){
          fit=multinom(x.~1,trace=F)
          pi_hat=fit$fitted.values
        }else{
          fit=multinom(x.~mat,trace=F)
          pi_hat=fit$fitted.values
        }
        dummy=sapply(colnames(pi_hat), function(x) as.numeric(as.character(x.)==x) )
        l_1=ncol(pi_hat)-1
        r=(dummy[,1:l_1,drop=F]-pi_hat[,1:l_1,drop=F])/( (pi_hat[,1:l_1,drop=F]*(1-pi_hat[,1:l_1,drop=F]))^0.5 )
        residual_x=rowSums(r)
    }
    
  } else{
    if (is.null(mat)){
      residual_x=resid(lm( x.~1 ))
    }else{
      residual_x=resid(lm.fit(y = x.,x = cbind(1,as.matrix(mat)) ))
    }
  }
  
  if (suffStat$type[y] == "c"){
    
    if (suffStat$level[y]==2){
      
        if (is.null(mat)){
          fit=glm(y.~1,family = "binomial")
          pi_hat=predict(fit, as.data.frame(rep(1,nrow(suffStat$dat_in))), type="response")
        }else{
          fit=logistf::logistf(y.~mat)
          pi_hat=fit$predict
          #pi_hat=predict(fit,as.data.frame(mat),type="response")
        }
        residual_y=(y.-pi_hat)/((pi_hat*(1-pi_hat))^0.5)
      
    }else{
        
        if (is.null(mat)){
          fit=multinom(y.~1,trace=F)
          pi_hat=fit$fitted.values
        }else{
          fit=multinom(y.~mat,trace=F)
          pi_hat=fit$fitted.values
        }
        dummy=sapply(colnames(pi_hat), function(x) as.numeric(as.character(y.)==x) )
        l_1=ncol(pi_hat)-1
        r=(dummy[,1:l_1,drop=F]-pi_hat[,1:l_1,drop=F])/( (pi_hat[,1:l_1,drop=F]*(1-pi_hat[,1:l_1,drop=F]))^0.5 )
        residual_y=rowSums(r)
    }
    
  } else{
    if (is.null(mat)){
      residual_y=resid(lm( y.~1 ))
    }else{
      residual_y=resid(lm.fit(y = y.,x = cbind(1,as.matrix(mat)) ))
    }
    
  }
  
  rxy = cbind(residual_x,residual_y)
  partial_observed = cor.test(rxy[, 1], rxy[, 2])$estimate
  return(list(rxy = rxy, partial_observed = partial_observed))
}


CalcultePvalue=function (rxy, partial_observed, nperm) 
{
  sam = replicate(nperm, sample(1:nrow(rxy), nrow(rxy)))
  permuate_rxy1 = apply(sam, 2, function(x) rxy[x, 1])
  partial_permute = abs(apply(permuate_rxy1, 2, function(x) cor.test(x, 
                                                                     rxy[, 2])$estimate))
  pvalue = length(which(abs(partial_observed) < partial_permute))/nperm
  return(pvalue)
}

ConditionalTestPermute=function (x, y, S, suffStat) 
{
  if (suffStat$type[x] == "c" | suffStat$type[y] == "c") {
    residual = GenerateResidual(x, y, S, suffStat)
    pvalue = CalcultePvalue(residual$rxy, residual$partial_observed, 
                              suffStat$nperm)
  }else pvalue = pcalg::gaussCItest(x, y, S, suffStat)
  return(pvalue)
}



