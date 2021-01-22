methods::setClass (
  "xhat",
  representation=representation(
    xhat="matrix"
  )
)
methods::setClass (
  "W",
  representation=representation(
    W="matrix"
  )
)
methods::setClass (
  "zc",
  representation=representation(
    zc="vector"
  )
)
methods::setClass (
  "zcchain",
  representation=representation(
    zcchain="list"
  )
)
methods::setClass (
  "zrchain",
  representation=representation(
    zrchain="vector"
  )
)
methods::setClass (
  "rho",
  representation=representation(
    rho="vector"
  )
)
methods::setClass (
  "rhochain",
  representation=representation(
    rhochain="list"
  )
)
methods::setClass (
  "pichain",
  representation=representation(
    pichain="vector"
  )
)
methods::setClass (
  "params",
  representation=representation(
    params="list"
  )
)
methods::setClass (
  "paramschain",
  representation=representation(
    paramschain="list"
  )
)
methods::setClass (
  "dlist",
  representation=representation(
    dlist="vector"
  )
)



# Result for co-clustering
methods::setClass (
  "ResultCoclustOrdinal",
  
  # Defining slot type
  representation = representation (
    V = "matrix",
    zr = "vector",
    pi = "vector",
    m = "vector",
    icl = "numeric",
    name="character"
  ),
  contains=c("W","zc","rho","params","paramschain","xhat","pichain","rhochain","zrchain","zcchain")
)


# Result for clustering
methods::setClass (
  "ResultClustOrdinal",
  
  # Defining slot type
  representation = representation (
    V = "matrix",
    zr = "vector",
    pi = "vector",
    m = "vector",
    icl = "numeric",
    name="character"
  ),
  contains=c("params","paramschain","xhat","zrchain","pichain")
)



# Result for classification
methods::setClass (
  "ResultClassifOrdinal",
  
  # Defining slot type
  representation = representation (
    V = "matrix",
    zr = "vector",
    pi = "vector",
    icl = "numeric",
    kr = "integer",
    kc = "vector",
    J = "vector",
    number_distrib = "integer",
    m = "vector",
    nbSEM = "integer",
    name="character"
  ),
  contains=c("W","zc","rho","params","dlist","xhat")
)

methods::setClass (
  "ResultPredictionOrdinal",
  
  # Defining slot type
  representation = representation (
    zr_topredict = "vector",
    V_topredict = "matrix"
  )
)