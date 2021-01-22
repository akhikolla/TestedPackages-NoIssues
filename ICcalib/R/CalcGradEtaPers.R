# Based on ICsurv:::fast.PH.ICsurv.EM
# @importFrom ICsurv Ispline
CalcGradEtaPers <- function(d1, d2, d3, Li, Ri,  knots, order, eta.g, eta.b, Q)
{
  n <- length(Ri)
  n.g <- length(eta.g)
  n.b <- length(eta.b)
  expQb <- as.vector(exp(Q%*%eta.b))
  # Portions of the code are taken from the ICsurv package
  ti <- c(Li[d1 == 0], Ri[d3 == 0])
  ti.max <- max(ti) + 1e-05
  ti.min <- min(ti) - 1e-05
  bRi <- t(ICsurv::Ispline(x = Ri, order = order, knots = knots))
  bLi <- t(ICsurv::Ispline(x = Li, order = order, knots = knots))
  GRi <- as.vector(bRi %*% eta.g)
  GLi <- as.vector(bLi %*% eta.g)
  HRi <-  as.vector(GRi*expQb )
  HLi <-  as.vector(GLi*expQb) 
  SRi <- exp(-HRi)
  SLi <- exp(-HLi)
  FRi <- 1-SRi
  FLi <- 1-SLi
  
  term.deriv.etab.d1 <- Q*(SRi*HRi/FRi)
  term.deriv.etab.d2 <- Q*(SRi*HRi -SLi*HLi)/(SLi-SRi)
  term.deriv.etab.d3 <- -Q*HLi
  
  term.deriv.etag.d1 <- bRi*(SRi*expQb/FRi)
  term.deriv.etag.d2 <- expQb*(SRi*bRi -SLi*bLi)/(SLi-SRi)
  term.deriv.etag.d3 <- -bLi*expQb
  
  deriv.ell.etag <- matrix(nrow = n, ncol =  n.g)
  deriv.ell.etab <- matrix(nrow = n, ncol =  n.b)
  
  deriv.ell.etab[d1==1,] <- term.deriv.etab.d1[d1==1]
  deriv.ell.etab[d2==1,] <- term.deriv.etab.d2[d2==1]
  deriv.ell.etab[d3==1,] <- term.deriv.etab.d3[d3==1]
  
  deriv.ell.etag[d1==1,] <- term.deriv.etag.d1[d1==1]
  deriv.ell.etag[d2==1,] <- term.deriv.etag.d2[d2==1]
  deriv.ell.etag[d3==1,] <- term.deriv.etag.d3[d3==1]
  
  deriv.ell.etas <- cbind(deriv.ell.etab,deriv.ell.etag)
  return(deriv.ell.etas)
}
