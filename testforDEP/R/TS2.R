
setClass("TS2",

  contains = "GeneralTest"

)

setValidity("TS2", function(object){
  if(object@p.opt == "dist")
    stop('No "dist" option for TS2, please use "MC".')
  if(object@p.opt == "table")
    stop('No "table" option for TS2, please use "MC".')

})

setMethod("test", signature(object = "TS2"), function(object){
  b = function(n ,x){
    if(n==1)
      return(sqrt(3)*(2*x-1))
    else if(n==2)
      return(sqrt(5)*(6*x^2-6*x+1))
    else if(n==3)
      return(sqrt(7)*(20*x^3-30*x^2+12*x-1))
    else if(n==4)
      return(sqrt(9)*(70*x^4-140*x^3+90*x^2-20*x+1))
  }

  p = object@pdata
  data = p[[ls(p)]]
  n = nrow(data)

  Kallenberg_TS2 = function(data){
    rankX = (rank(data[,1])-0.5)/n
    rankY = (rank(data[,2])-0.5)/n

    T_sub = c()
    for(k in 1:4){
      T_sub[k] = ((t(b(k ,rankX)) %*% b(k, rankY))^2)/n
    }

    Tk = c()
    for(i in 1:4){
      Tk[i] = sum(T_sub[1:i])
    }

    S2 = which.max(Tk - c(1:4)*log(n))
    return(Tk[S2])

  }
  #Ts
  TS2 = Kallenberg_TS2(data)

  #MC

    sn = object@num.MC

    TS2s = c()
    for(i in 1:sn){
      if(object@set.seed){set.seed(i)}
      data = cbind(rnorm(n, 0, 1), rnorm(n, 0, 1))
      TS2s[i] = Kallenberg_TS2(data)
    }

    NGE = length(which(TS2s>TS2))
    pv = NGE/(sn+1)


  #Bootstrap

  if(object@BS.CI != 0){
    times = 1000
    TS2s = c()
    for(i in 1:times){
      if(object@set.seed){set.seed(i)}
      BSdata = data[sample(1:n, n, replace = TRUE),]
      TS2s[i] = Kallenberg_TS2(BSdata)
    }
    CI = getCI(TS2, TS2s, object@BS.CI)
    return(new("testforDEP_result", TS = TS2, p_value = pv, CI = CI))
  }
  else{
    return(new("testforDEP_result", TS = TS2, p_value = pv))
  }


})





