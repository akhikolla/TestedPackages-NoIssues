
setClass("V",

  contains = "GeneralTest"

)

setValidity("V", function(object){
  if(object@p.opt == "dist")
    stop('No "dist" option for V, please use "MC".')
  if(object@p.opt == "table")
    stop('No "table" option for V, please use "MC".')

})

setMethod("test", signature(object = "V"), function(object){
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

  Kallenberg_V = function(data){

    getV = function(data, i, j){
      rankX = (rank(data[,1])-0.5)/n
      rankY = (rank(data[,2])-0.5)/n

      return((t(b(i, rankX)) %*% b(j, rankY))^2/n)
    }

    Vij = c()
    c = 1
    for(i in 1:2)
      for(j in 1:2){
        Vij[c] = getV(data, i, j)
        c = c + 1
      }

    if(max(Vij[2:4]) < log(n))
      return(Vij[1])
    else
      return(Vij[1] + max(Vij[2:4]))
  }

    V = Kallenberg_V(data)

    #MC

    sn = object@num.MC

    Vs = c()
    for(i in 1:sn){
      if(object@set.seed){set.seed(i)}
      data = cbind(rnorm(n, 0, 1), rnorm(n, 0, 1))
      Vs[i] = Kallenberg_V(data)
    }

    NGE = length(which(Vs>V))
    pv = NGE/(sn+1)


  #Bootstrap

  if(object@BS.CI != 0){
    times = 1000
    Vs = c()
    for(i in 1:times){
      if(object@set.seed){set.seed(i)}
      BSdata = data[sample(1:n, n, replace = TRUE),]
      Vs[i] = Kallenberg_V(BSdata)
    }
    CI = getCI(V, Vs, object@BS.CI)
    return(new("testforDEP_result", TS = V, p_value = pv, CI = CI))
  }
  else{
    return(new("testforDEP_result", TS = V, p_value = pv))
  }


})





