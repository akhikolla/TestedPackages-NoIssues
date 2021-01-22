
# output of multinom_BIC is -0.5*BIC
multinom_BIC=function(x,y,weights){
  formula=y~x
  fit=nnet::multinom(formula,weights=weights,trace=F)
  p=fit$rank
  AIC=fit$AIC
  out=-0.5*AIC+p-0.5*p*log(length(weights)) # -0.5*BIC
  return( list("out"=out) )
}

