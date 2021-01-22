rMU <- function(n)
{
  img <- readPNG(system.file("images/mu.png", package = "ICtest")[1])[, , 1]
  img <- t(img)
  img <- img[, ncol(img):1]
  x <- imageToData(img, n)
  x
}