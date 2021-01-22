ggscreeplot <- function (x, type = "barplot", main = deparse(substitute(x)),
                              ylab = "criterion", xlab = "component")
{
  type <- match.arg(type, c("barplot", "lines"))
  
  index = 1:length(x$D)
  D = x$D
  dummy <- data.frame(index = index, D = D)
  
  if (type == "barplot") {
    ggplot(dummy, aes(index, D)) +
      geom_bar(stat = "identity") +
      xlab(xlab) +
      ylab(ylab) +
      ggtitle(main)
    
  } else {
    ggplot(dummy, aes(index, D)) +
      geom_line() +
      geom_point() +
      xlab(xlab) +
      ylab(ylab) +
      ggtitle(main)
  }
}

