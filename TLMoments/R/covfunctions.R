# This file contains internal functions used to calculate covariances
# of estimations of PWMs, TLMoments, Parameters, and Quantiles.

# Generate matrix for conversion of lambda-cov to tau-cov
CovLambdaToTau  <- function(l) {

  K <- length(l)

  if (K > 3) {
    A <- rbind(
      c(-l[2]/l[1]^2, rep(0, K-2)),
      c(1/l[1], -l[-c(1,2)]/l[2]^2),
      cbind(0, diag(rep(1/l[2], K-2)))
    )
  } else if (K == 3) {
    A <- rbind(
      c(-l[2]/l[1]^2, rep(0, K-2)),
      c(1/l[1], -l[-c(1,2)]/l[2]^2),
      cbind(0, 1/l[2])
    )
  } else {
    A <- rbind(-l[2]/l[1]^2, 1/l[1])
  }

  A
}


v <- function(r, delta, scale, shape) {
  G <- function(x, a, b = a) hypergeo::hypergeo(-b, -2*a, 1-a, -x)

  if (delta == 0) {
    out <- scale^2 * (-shape)^(-2) * (r+1)^(2*shape) * ( gamma(1-2*shape) * G(r/(r+1), shape) - gamma(1-shape)^2 )
  } else if (delta == 1) {
    out <- .5 * scale^2 * (-shape)^(-2) * ( (r+2)^(2*shape) * gamma(1-2*shape) * G(r/(r+2), shape) + (r+1)^shape * ( (r+1)^shape - 2*(r+2)^shape ) * gamma(1-shape)^2 )
  } else {
    s <- delta
    out <- .5 * scale^2 * (-shape)^(-2) * gamma(1-2*shape) * (
      (2*r+s+1)^(2*shape) * G(-(r+s)/(2*r+s+1), shape, -1) -
        (2*r+s+2)^(2*shape) * G(-(r+s+1)/(2*r+s+2), shape, -1) +
        (2*r+s+1)^(2*shape) * G(-(r)/(2*r+s+1), shape, -1) -
        (2*r+s+2)^(2*shape) * G(-(r+1)/(2*r+s+2), shape, -1)
    )
  }
  Re(out)
}

# @return Symmetric Matrix of dimensions determined by \code{betas}.
PWMCov_GEV <- function(betas = 0:3, scale, shape) {
  if (abs(shape) <= 1e-06) shape <- 1e-6
  mirs <- min(betas)
  mars <- max(betas)
  out <- matrix(NA, nrow = mars+1, ncol = mars+1)
  for (r in mirs:mars) for (s in r:mars) {
    if (r == s) {
      out[r+1, s+1] <- v(r, s-r, scale, shape)
    } else {
      out[r+1, s+1] <- out[s+1, r+1] <- v(r, s-r, scale, shape)
    }
  }
  out <- out[betas+1, betas+1, drop = FALSE]

  colnames(out) <- rownames(out) <- paste0("beta", betas)
  out
}

parametricPWMCov <- function(distribution, order = 0:3, ...) {
  p <- list(...)
  if (distribution == "gev") {
    if (!all(c("scale", "shape") %in% names(p))) stop("scale and shape must be given. ")

    return(PWMCov_GEV(order, p[["scale"]], p[["shape"]]))
  }

  stop("No functions available for given distribution. ")
}


# Matrices to calculate Cov of parameter estimations out of TL(0,0)- or TL(0,1)-Covs.
# Need L2 and L3 with corresponding trimmings.
CovTLtoPara <- function(distr, l2, l3, leftrim, rightrim) {
  if (distr != "gev") stop("Only GEV for now. ")

  if (distr == "gev") {
    if (leftrim == 0 & rightrim == 0) {
      GEV_00(l2, l3)
    } else if (leftrim == 0 & rightrim == 1) {
      GEV_01(l2, l3)
    } else if (leftrim == 0 & rightrim == 2) {
      GEV_02(l2, l3)
    }
  } else {
    stop("Unsupported distr or trimming. ")
  }
}

GEV_00 <- function(l2, l3) {

  a <- -7.859
  b <- -2.9554
  kappa <- log(2)/log(3)

  zeta <- 1/(l3+3*l2)
  eta <- 2*l2*zeta - kappa
  theta <- 2*zeta - 6*l2*zeta^2

  # shape
  g1 <- 0
  g2 <- -2*l3*(2*b*kappa*l3-a*l3+6*b*kappa*l2-4*b*l2-3*a*l2)*zeta^3
  g3 <- 2*l2*(2*b*kappa*l3-a*l3+6*b*kappa*l2-4*b*l2-3*a*l2)*zeta^3
  r3 <- c(g1, g2, g3)

  # scale
  g2 <- (log(2)*(-b*eta^2-a*eta)*2^(b*eta^2+a*eta)*l2*(2*b*eta*theta+a*theta))/((1-2^(b*eta^2+a*eta))^2*gamma(-b*eta^2-a*eta+1))-((-b*eta^2-a*eta)*digamma(-b*eta^2-a*eta+1)*l2*(-2*b*eta*theta-a*theta))/((1-2^(b*eta^2+a*eta))*gamma(-b*eta^2-a*eta+1))+(l2*(-2*b*eta*theta-a*theta))/((1-2^(b*eta^2+a*eta))*gamma(-b*eta^2-a*eta+1))+(-b*eta^2-a*eta)/((1-2^(b*eta^2+a*eta))*gamma(-b*eta^2-a*eta+1))
  g3 <- -((-b*eta^2-a*eta)*digamma(-b*eta^2-a*eta+1)*l2*(4*b*eta*l2*zeta^2+2*a*l2*zeta^2))/((1-2^(b*eta^2+a*eta))*gamma(-b*eta^2-a*eta+1))+(l2*(4*b*eta*l2*zeta^2+2*a*l2*zeta^2))/((1-2^(b*eta^2+a*eta))*gamma(-b*eta^2-a*eta+1))+(log(2)*(-b*eta^2-a*eta)*2^(b*eta^2+a*eta)*l2*(-4*b*eta*l2*zeta^2-2*a*l2*zeta^2))/((1-2^(b*eta^2+a*eta))^2*gamma(-b*eta^2-a*eta+1))
  r2 <- c(0, g2, g3)

  # loc
  g2 <- -(log(2)*(-b*eta^2-a*eta)*2^(b*eta^2+a*eta)*(gamma(-b*eta^2-a*eta+1)-1)*l2*(2*b*eta*theta+a*theta))/((b*eta^2+a*eta)*(1-2^(b*eta^2+a*eta))^2*gamma(-b*eta^2-a*eta+1))+((-b*eta^2-a*eta)*(gamma(-b*eta^2-a*eta+1)-1)*l2*(2*b*eta*theta+a*theta))/((b*eta^2+a*eta)^2*(1-2^(b*eta^2+a*eta))*gamma(-b*eta^2-a*eta+1))+
    ((-b*eta^2-a*eta)*digamma(-b*eta^2-a*eta+1)*(gamma(-b*eta^2-a*eta+1)-1)*l2*(-2*b*eta*theta-a*theta))/((b*eta^2+a*eta)*(1-2^(b*eta^2+a*eta))*gamma(-b*eta^2-a*eta+1))-((gamma(-b*eta^2-a*eta+1)-1)*l2*(-2*b*eta*theta-a*theta))/((b*eta^2+a*eta)*(1-2^(b*eta^2+a*eta))*gamma(-b*eta^2-a*eta+1))-((-b*eta^2-a*eta)*digamma(-b*eta^2-a*eta+1)*l2*(-2*b*eta*theta-a*theta))/((b*eta^2+a*eta)*(1-2^(b*eta^2+a*eta)))-
    ((-b*eta^2-a*eta)*(gamma(-b*eta^2-a*eta+1)-1))/((b*eta^2+a*eta)*(1-2^(b*eta^2+a*eta))*gamma(-b*eta^2-a*eta+1))

  g3 <- ((-b*eta^2-a*eta)*digamma(-b*eta^2-a*eta+1)*(gamma(-b*eta^2-a*eta+1)-1)*l2*(4*b*eta*l2*zeta^2+2*a*l2*zeta^2))/((b*eta^2+a*eta)*(1-2^(b*eta^2+a*eta))*gamma(-b*eta^2-a*eta+1))-((gamma(-b*eta^2-a*eta+1)-1)*l2*(4*b*eta*l2*zeta^2+2*a*l2*zeta^2))/((b*eta^2+a*eta)*(1-2^(b*eta^2+a*eta))*gamma(-b*eta^2-a*eta+1))-
    ((-b*eta^2-a*eta)*digamma(-b*eta^2-a*eta+1)*l2*(4*b*eta*l2*zeta^2+2*a*l2*zeta^2))/((b*eta^2+a*eta)*(1-2^(b*eta^2+a*eta)))-(log(2)*(-b*eta^2-a*eta)*2^(b*eta^2+a*eta)*(gamma(-b*eta^2-a*eta+1)-1)*l2*(-4*b*eta*l2*zeta^2-2*a*l2*zeta^2))/((b*eta^2+a*eta)*(1-2^(b*eta^2+a*eta))^2*gamma(-b*eta^2-a*eta+1))+
    ((-b*eta^2-a*eta)*(gamma(-b*eta^2-a*eta+1)-1)*l2*(-4*b*eta*l2*zeta^2-2*a*l2*zeta^2))/((b*eta^2+a*eta)^2*(1-2^(b*eta^2+a*eta))*gamma(-b*eta^2-a*eta+1))

  r1 <- c(1, g2, g3)

  rbind(r1, r2, r3)
}

GEV_01 <- function(l2, l3) {
  kappa <- (2*log(2)-log(3))/(3*log(3)-2*log(4))
  a <- -8.567394
  b <- 0.675969

  zeta <- 10/(9*(l3+2*l2))
  eta <- 10*l2/(9*(l3+2*l2))
  theta <- 20*l2/(9*(l3+2*l2)^2)

  r3 <- (a + 2 * b * eta - kappa) * 10/(9*(l3+2*l2)^2) * c(
    0, l3, -l2
  )

  g1 <- 0
  g2 <- -(2*l2*(log(3)*3^(a*(eta-kappa)+b*(eta-kappa)^2)*(2*b*(eta-kappa)*(zeta-theta)+a*(zeta-theta))-log(2)*2^(a*(eta-kappa)+b*(eta-kappa)^2+1)*(2*b*(eta-kappa)*(zeta-theta)+a*(zeta-theta))))/(3*gamma(-a*(eta-kappa)-b*(eta-kappa)^2)*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1)^2)-(2*digamma(-a*(eta-kappa)-b*(eta-kappa)^2)*l2*(-2*b*(eta-kappa)*(zeta-theta)-a*(zeta-theta)))/(3*gamma(-a*(eta-kappa)-b*(eta-kappa)^2)*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1))+    2/(3*gamma(-a*(eta-kappa)-b*(eta-kappa)^2)*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1))
  g3 <- -(2*l2*(log(3)*3^(a*(eta-kappa)+b*(eta-kappa)^2)*(-(20*b*(eta-kappa)*l2)/(9*(l3+2*l2)^2)-(10*a*l2)/(9*(l3+2*l2)^2))-log(2)*2^(a*(eta-kappa)+b*(eta-kappa)^2+1)*(-(20*b*(eta-kappa)*l2)/(9*(l3+2*l2)^2)-(10*a*l2)/(9*(l3+2*l2)^2))))/(3*gamma(-a*(eta-kappa)-b*(eta-kappa)^2)*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1)^2)  -(2*digamma(-a*(eta-kappa)-b*(eta-kappa)^2)*l2*((20*b*(eta-kappa)*l2)/(9*(l3+2*l2)^2)+(10*a*l2)/(9*(l3+2*l2)^2)))/(3*gamma(-a*(eta-kappa)-b*(eta-kappa)^2)*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1))

  r2 <- c(g1, g2, g3)

  g1 <- 1
  g2 <- (2*(2^(a*(eta-kappa)+b*(eta-kappa)^2)-2)*l2*(log(3)*3^(a*(eta-kappa)+b*(eta-kappa)^2)*(2*b*(eta-kappa)*(zeta-theta)+a*(zeta-theta))-log(2)*2^(a*(eta-kappa)+b*(eta-kappa)^2+1)*(2*b*(eta-kappa)*(zeta-theta)+a*(zeta-theta)))  )/(3*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1)^2)-  (2*l2*(log(3)*3^(a*(eta-kappa)+b*(eta-kappa)^2)*(2*b*(eta-kappa)*(zeta-theta)+a*(zeta-theta))-log(2)*2^(a*(eta-kappa)+b*(eta-kappa)^2+1)*(2*b*(eta-kappa)*(zeta-theta)+a*(zeta-theta))))/(3*gamma(-a*(eta-kappa)-b*(eta-kappa)^2)*(a*(eta-kappa)+b*(eta-kappa)^2)*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1)^2)-  (log(2)*2^(a*(eta-kappa)+b*(eta-kappa)^2+1)*l2*(2*b*(eta-kappa)*(zeta-theta)+a*(zeta-theta)))/(3*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1))-  (2*l2*(2*b*(eta-kappa)*(zeta-theta)+a*(zeta-theta)))/(3*gamma(-a*(eta-kappa)-b*(eta-kappa)^2)*(a*(eta-kappa)+b*(eta-kappa)^2)^2*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1))-  (2*digamma(-a*(eta-kappa)-b*(eta-kappa)^2)*l2*(-2*b*(eta-kappa)*(zeta-theta)-a*(zeta-theta)))/(3*gamma(-a*(eta-kappa)-b*(eta-kappa)^2)*(a*(eta-kappa)+b*(eta-kappa)^2)*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1))-  (2*(2^(a*(eta-kappa)+b*(eta-kappa)^2)-2))/(3*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1))+    2/(3*gamma(-a*(eta-kappa)-b*(eta-kappa)^2)*(a*(eta-kappa)+b*(eta-kappa)^2)*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1))
  g3 <- (2*(2^(a*(eta-kappa)+b*(eta-kappa)^2)-2)*l2*(log(3)*3^(a*(eta-kappa)+b*(eta-kappa)^2)*(-(20*b*(eta-kappa)*l2)/(9*(l3+2*l2)^2)-(10*a*l2)/(9*(l3+2*l2)^2))-log(2)*                                                       2^(a*(eta-kappa)+b*(eta-kappa)^2+1)*(-(20*b*(eta-kappa)*l2)/(9*(l3+2*l2)^2)-(10*a*l2)/(9*(l3+2*l2)^2))))/(3*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1)^2)-  (2*l2*(log(3)*3^(a*(eta-kappa)+b*(eta-kappa)^2)*(-(20*b*(eta-kappa)*l2)/(9*(l3+2*l2)^2)-(10*a*l2)/(9*(l3+2*l2)^2))-log(2)*2^(a*(eta-kappa)+b*(eta-kappa)^2+1)*(-(20*b*(eta-kappa)*l2)/(9*(l3+2*l2)^2)-(10*a*l2)/(9*(l3+2*l2)^2))))/(3*gamma(-a*(eta-kappa)-b*(eta-kappa)^2)*(a*(eta-kappa)+b*(eta-kappa)^2)*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1)^2)  -(2*digamma(-a*(eta-kappa)-b*(eta-kappa)^2)*l2*((20*b*(eta-kappa)*l2)/(9*(l3+2*l2)^2)+(10*a*l2)/(9*(l3+2*l2)^2)))/(3*gamma(-a*(eta-kappa)-b*(eta-kappa)^2)*(a*(eta-kappa)+b*(eta-kappa)^2)*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1))-  (log(2)*2^(a*(eta-kappa)+b*(eta-kappa)^2+1)*l2*(-(20*b*(eta-kappa)*l2)/(9*(l3+2*l2)^2)-(10*a*l2)/(9*(l3+2*l2)^2)))/(3*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1))-  (2*l2*(-(20*b*(eta-kappa)*l2)/(9*(l3+2*l2)^2)-(10*a*l2)/(9*(l3+2*l2)^2)))/(3*gamma(-a*(eta-kappa)-b*(eta-kappa)^2)*(a*(eta-kappa)+b*(eta-kappa)^2)^2*(-2^(a*(eta-kappa)+b*(eta-kappa)^2+1)+3^(a*(eta-kappa)+b*(eta-kappa)^2)+1))

  r1 <- c(g1, g2, g3)

  rbind(r1, r2, r3)
}

GEV_02 <- function(l2, l3) {
  kappa <- (3*log(5)-8*log(4)+6*log(3)) / (log(4)-3*log(3)+3*log(2)) - 2
  a <- 2.468959
  b <- -1.130074
  c <- 0.635912

  zeta <- 6/5 * l3/l2
  eta <- zeta - kappa
  mu <- c*eta^3 + b*eta^2 + a*eta
  theta <- 18*c*eta^2/(5*l2) + 12*b*eta/(5*l2) + 6*a/(5*l2)

  r3 <- theta * c(
    0,
    -l3/l2,
    1
  )

  g1 <- 0
  g2 <- (((-l3*3^(mu+1)+l3*4^mu+3*l3*2^mu-l3)*digamma(-mu)+log(3)*l3*3^(mu+1)-log(4)*l3*4^mu-3*log(2)*l3*2^mu)*theta+3^(mu+1)-4^mu-3*2^mu+1)/((2*3^(2*mu+2)+9*2^(2*mu+1)+4^mu*(3*2^(mu+2)-4*3^(mu+1)-4)-3*2^(mu+2)+3^mu*(12-9*2^(mu+2))+2*4^(2*mu)+2)*gamma(-mu))
  g3 <- -(((-l2*3^(mu+1)+l2*4^mu+3*l2*2^mu-l2)*digamma(-mu)+log(3)*l2*3^(mu+1)-log(4)*l2*4^mu-3*log(2)*l2*2^mu)*theta)/((2*3^(2*mu+2)+9*2^(2*mu+1)+4^mu*(3*2^(mu+2)-4*3^(mu+1)-4)-3*2^(mu+2)+3^mu*(12-9*2^(mu+2))+2*4^(2*mu)+2)*gamma(-mu))

  r2 <- c(g1, g2, g3)

  g1 <- 1
  g2 <- -(((3*log(2)*l3*mu^2*2^(mu+1)+((log(4)-log(3))*l3*mu^2*3^mu+(3*log(2)-3*log(4))*l3*mu^2*2^mu+3*log(4)*l3*mu^2)*4^mu+((6*log(3)-6*log(2))*l3*mu^2*2^mu-8*log(3)*l3*mu^2)*3^mu)*gamma(-mu)+(l3*mu*3^(mu+1)-l3*mu*4^mu-3*l3*mu*2^mu+l3*mu)*digamma(-mu)+(log(4)*l3*mu+l3)*4^mu+(-3*log(3)*l3*mu-3*l3)*3^mu+(3*log(2)*l3*mu+3*l3)*2^mu-l3)*theta+(-mu^2*3^(2*mu+1)+3^mu*(3*mu^2*2^(mu+2)-10*mu^2)+3*mu^2*2^(mu+2)-9*mu^2*2^(2*mu)+(mu^2*3^mu-3*mu^2*2^mu+3*mu^2)*4^mu-3*mu^2)*gamma(-mu)-mu*3^(mu+1)+mu*4^mu+3*mu*2^mu-mu)/((2*mu^2*3^(2*mu+2)+9*mu^2*2^(2*mu+1)+4^mu*(3*mu^2*2^(mu+2)-4*mu^2*3^(mu+1)-4*mu^2)+3^mu*(12*mu^2-9*mu^2*2^(mu+2))-3*mu^2*2^(mu+2)+2*mu^2*4^(2*mu)+2*mu^2)*gamma(-mu))
  g3 <-  (((3*log(2)*l2*mu^2*2^(mu+1)+((log(4)-log(3))*l2*mu^2*3^mu+(3*log(2)-3*log(4))*l2*mu^2*2^mu+3*log(4)*l2*mu^2)*4^mu+((6*log(3)-6*log(2))*l2*mu^2*2^mu-8*log(3)*l2*mu^2)*3^mu)*gamma(-mu)+(l2*mu*3^(mu+1)-l2*mu*4^mu-3*l2*mu*2^mu+l2*mu)*digamma(-mu)+(log(4)*l2*mu+l2)*4^mu+(-3*log(3)*l2*mu-3*l2)*3^mu+(3*log(2)*l2*mu+3*l2)*2^mu-l2)*theta)/((2*mu^2*3^(2*mu+2)+9*mu^2*2^(2*mu+1)+4^mu*(3*mu^2*2^(mu+2)-4*mu^2*3^(mu+1)-4*mu^2)+3^mu*(12*mu^2-9*mu^2*2^(mu+2))-3*mu^2*2^(mu+2)+2*mu^2*4^(2*mu)+2*mu^2)*gamma(-mu))

  r1 <- c(g1, g2, g3)

  rbind(r1, r2, r3)
}


# Matrices to convert parameter estimations to quantile estimations
CovParamtoQuan <- function(distr, param, p) {
  if (distr == "gev") {
    g <- param["shape"]
    sigma <- param["scale"]
    xi <- -log(p)
    return(
      cbind(
        1,
        1/g * (xi^(-g) - 1), #-((-log(p))^g - 1)/(g * (-log(p))^g),
        sigma/g * (g^(-1) - xi^(-g) * (log(xi) + g^(-1))) #-((g*log(-log(p))-(-log(p))^g+1)*sigma)/(g^2*(-log(p))^g)
      )
    )
  } else if (distr == "gumbel") {
    return(
      cbind(
        1,
        -log(-log(p))
      )
    )
  } else {
    stop("distr currently unsupported. ")
  }
}
