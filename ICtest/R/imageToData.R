imageToData <- function(image.matrix, n){
  side.length <- dim(image.matrix)[1]
  image.data = matrix(0, n, 2)
  for(i in 1:n){
    hit <- FALSE
    # Iterate until a "hit"
    while(!hit){
      single.try <- runif(2)
      # Find the corresponding pixel
      xcoord <- ceiling(side.length*(single.try[1]))
      ycoord <- ceiling(side.length*(single.try[2]))
      # The darker the region, the more probable to "hit"
      if(rbinom(1, 1, image.matrix[xcoord, ycoord]) == 0) hit <- TRUE
    }
    image.data[i, ] <- single.try
  }
  return(image.data)
}