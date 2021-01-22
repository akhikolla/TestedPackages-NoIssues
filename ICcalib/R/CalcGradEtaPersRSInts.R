# Based on ICsurv:::fast.PH.ICsurv.EM
# @importFrom ICsurv Ispline
CalcGradEtaPersRSInts <- function(d1, d2, d3, Li, Ri, Q, fit.cox.rs.ints, pts.for.ints, tm, n.etas.per.fit)
{
  n <- length(Ri)
  n.fits <- length(pts.for.ints)
  deriv.ell.etas <- matrix(nrow = n, ncol = sum(n.etas.per.fit), 0)
  for (j in 1:n.fits)
  {
  point <- pts.for.ints[j]
  fit.cox.int <- fit.cox.rs.ints[[j]]
  eta.b <- fit.cox.int$b
  eta.g <- fit.cox.int$g
  n.g <- length(eta.g)
  n.b <- length(eta.b)
  knots <- fit.cox.int$knots
  order <- fit.cox.int$order
  in.risk.set <- tm >= point
  n.set <- sum(in.risk.set)
  Li.int <- Li[in.risk.set]
  Ri.int <- Ri[in.risk.set]
  d1.int <- d1[in.risk.set]
  d2.int <- d2[in.risk.set]
  d3.int <- d3[in.risk.set]
  Q.int <- Q[in.risk.set,]
  expQb <- as.vector(exp(Q.int%*%eta.b))
  # Portions of the code are taken from the ICsurv package
  ti <- c(Li.int[d1.int == 0], Ri.int[d3.int == 0])
  ti.max <- max(ti) + 1e-05
  ti.min <- min(ti) - 1e-05
  bRi <- t(ICsurv::Ispline(x = Ri.int, order = order, knots = knots))
  bLi <- t(ICsurv::Ispline(x = Li.int, order = order, knots = knots))
  GRi <- as.vector(bRi %*% eta.g)
  GLi <- as.vector(bLi %*% eta.g)
  HRi <-  as.vector(GRi*expQb)
  HLi <-  as.vector(GLi*expQb) 
  SRi <- exp(-HRi)
  SLi <- exp(-HLi)
  FRi <- 1-SRi
  FLi <- 1-SLi
  
  term.deriv.etab.d1 <- Q.int*(SRi*HRi/FRi)
  term.deriv.etab.d2 <- Q.int*(SRi*HRi -SLi*HLi)/(SLi-SRi)
  term.deriv.etab.d3 <- -Q.int*HLi
  
  term.deriv.etag.d1 <- bRi*(SRi*expQb/FRi)
  term.deriv.etag.d2 <- expQb*(SRi*bRi -SLi*bLi)/(SLi-SRi)
  term.deriv.etag.d3 <- -bLi*expQb
  
  deriv.ell.etag.int <- matrix(nrow = n.set, ncol =  n.g)
  deriv.ell.etab.int <- matrix(nrow = n.set, ncol =  n.b)
  
  deriv.ell.etab.int[d1.int==1,] <- term.deriv.etab.d1[d1.int==1]
  deriv.ell.etab.int[d2.int==1,] <- term.deriv.etab.d2[d2.int==1]
  deriv.ell.etab.int[d3.int==1,] <- term.deriv.etab.d3[d3.int==1]
  
  deriv.ell.etag.int[d1.int==1,] <- term.deriv.etag.d1[d1.int==1]
  deriv.ell.etag.int[d2.int==1,] <- term.deriv.etag.d2[d2.int==1]
  deriv.ell.etag.int[d3.int==1,] <- term.deriv.etag.d3[d3.int==1]
  
  deriv.ell.etas.int <- cbind(deriv.ell.etab.int,deriv.ell.etag.int)
  if (j > 1) {
  deriv.ell.etas[in.risk.set, (sum(n.etas.per.fit[1:(j-1)]) + 1):sum(n.etas.per.fit[1:j])] <- deriv.ell.etas.int 
  } else {
  deriv.ell.etas[in.risk.set, 1:n.etas.per.fit[1]] <- deriv.ell.etas.int 
  }
  }
return(deriv.ell.etas)
}
