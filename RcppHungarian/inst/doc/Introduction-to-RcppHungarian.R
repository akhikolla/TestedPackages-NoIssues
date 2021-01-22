## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center"
)

## ------------------------------------------------------------------------
set.seed(593903)

library(RcppHungarian)
library(ggplot2)

# humans are rows (5), dogs are columns (4)
cost <- rbind(c(1,5,2,19), 
              c(4,0,3,10), 
              c(6,9,6,1), 
              c(9,8,1,3), 
              c(6,1,1,1))
HungarianSolver(cost)

## ----include=FALSE-------------------------------------------------------
plot_arrows <- function(dat, col){
  dat <- data.frame(dat, group=col)
  colnames(dat) <- c("Coordinate 1", "Coordinate 2", "Eigenvector")
  ggplot(dat) +
    geom_segment(aes(xend=`Coordinate 1`, yend=`Coordinate 2`, x=0, y=0, color=Eigenvector), 
                 arrow=arrow(length=unit(0.1, "inches")), alpha=0.7) +
    theme_minimal()+
    xlab("Coordinate 1") +
    ylab("Coordinate 2")+
    scale_color_brewer(palette="Set1") +
    ylim(c(-1,1))+
    xlim(c(-1,1))
}

## ----fig.height=5, fig.width=5-------------------------------------------
# Create a distribution of matricies that have identical eigenvectors but there
# is randomness in the eigenvalues and a tiny bit of added noise from a wishart
const <- 1/sqrt(2)
V <- cbind(c(const,const), c(const,-const))
Sigma <- array(0, dim=c(2, 2, 100))
for (i in 1:dim(Sigma)[3]){
  Sigma[,,i] <- V %*% diag(sample(c(1, 1.3))) %*% t(V) + rWishart(1, 50, 0.002*diag(2))[,,1]
}

## ----fig.height=5, fig.width=5-------------------------------------------
plot_arrows(V, c("1", "2"))

## ----fig.height=5, fig.width=5-------------------------------------------
# Decompose each sample
Sigma.vectors <- array(0, dim=dim(Sigma))
for (i in 1:dim(Sigma)[3]){
  decomp <- eigen(Sigma[,,i])
  Sigma.vectors[,,i] <- decomp$vectors
}


# Plot distribution of the first two Eigenvectors
dat <- rbind(t(Sigma.vectors[,1,]), 
             t(Sigma.vectors[,2,]))
col <- c(rep("1", dim(Sigma)[3]), rep("2", dim(Sigma)[3]))
plot_arrows(dat, col)

## ----fig.height=5, fig.width=5-------------------------------------------

# Here is our distance function - calculates the pairwise distance
# between two sets of vectors. 
cost_distance <- function(V1, V2){
  n1 <- ncol(V1)
  n2 <- ncol(V2)
  D <-  as.matrix(dist(t(cbind(V1, V2))))
  return(D[1:n1, (n1+1):(n1+n2)])
}

# Make some empty containers to store results
Sigma.vectors <- array(0, dim=dim(Sigma))

# Initialize n=1 case - this is our reference against which 
# we will match the rest of the samples
decomp <- eigen(Sigma[,,1])
Sigma.vectors[,,1] <- decomp$vectors

# Now match up each sample with the reference
for (i in 2:dim(Sigma)[3]){
  decomp <- eigen(Sigma[,,i])
  Sigma.vectors[,,i] <- decomp$vectors
  
  # Here we add a second (negated) copy of the eigenvectors in the sample
  Sigma.vectors.expanded <-  cbind(Sigma.vectors[,,i], -Sigma.vectors[,,i])
  
  # Here we compute the "cost matrix" (note this is rectangular)
  # because we have the duplicated (negated) vectors from the sample
  cost <- cost_distance(Sigma.vectors[,,1],Sigma.vectors.expanded)

  # This is the key line where we solve the matching problem
  matching <- RcppHungarian::HungarianSolver(cost)
  
  # This is the key line were we reorder the eigenvectors / select if we wanted
  # the negative or positive versions
  Sigma.vectors[,,i] <- Sigma.vectors.expanded[,matching$pairs[,2]]
}

# Now we check our work visually
dat <- rbind(t(Sigma.vectors[,1,]), 
             t(Sigma.vectors[,2,]))
col <- c(rep("1", dim(Sigma)[3]), rep("2", dim(Sigma)[3]))
plot_arrows(dat, col)

## ------------------------------------------------------------------------
plot_arrows

