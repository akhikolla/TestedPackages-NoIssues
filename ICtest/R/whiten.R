whiten <-
function(x){
  sweep(x, 2, colMeans(x), '-')%*%symmetricPower_C(cov(x), -0.5)
}
