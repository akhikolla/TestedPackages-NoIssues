#' Distance correlation T-test
#' It uses the result of U-statistic distance correlation to test independence for
#' high dimensional data
#' @param x data of x
#' @param y data of y
#' @importFrom stats pt
#' @examples
#' n = 200
#' x = rnorm(n)
#' y = rnorm(n)
#' res = dcor.ttest(x,y)
#' @export

dcor.ttest<-function(x,y){

  if(!is.matrix(x)) x = as.matrix(x)
  if(!is.matrix(y)) y = as.matrix(y)
  if(nrow(x)!=nrow(y)) stop("x and y should have same numeber of rows.")
  stat = dcor(x,y,'U')
  n = nrow(x)
  v = n*(n-3)/2
  tstat = sqrt(v-1)*stat/sqrt(1-stat^2)
  df = v-1
  pval = 1-pt(tstat,df)
  res = list(statistic=tstat, dcor=stat, df=df, p.value = pval)
  return(res)

}


