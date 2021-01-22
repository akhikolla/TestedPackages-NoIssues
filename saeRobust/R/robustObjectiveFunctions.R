scoreRobustBeta <- function(y, x, matV, psi) {
    # Helper functions
    resid <- function(beta) U$sqrtInv() %*% (y - x %*% beta)
    D <- function(beta) Diagonal(x = psi(resid(beta), deriv = TRUE))

    # Precalculations - they only have to be done once
    U <- matU(matV$V())
    memP0 <- crossprod(x, matV$VInv())
    memP1 <- memP0 %*% U$sqrt()

    f <- function(beta) memP1 %*% psi(resid(beta))
    f1 <- function(beta) - memP0 %*% D(beta) %*% x

    list(f = f, f1 = f1)
}

fixedPointRobustBeta <- function(y, x, matV, psi) {
    makeMatA <- matAConst(y, x, matV, psi)
    function(beta) {
        as.numeric(makeMatA(beta) %*% y)
    }
}

robustObjectiveDelta <- function(y, x, beta, matVFun, psi, K, derivSelect) {
  # This is the squared estimation equation for a variance parameter. It can be
  # used when the fixed point for delta
  function(rho) {

    matV <- matVFun(rho)
    U <- matU(matV$V())
    psiResid <- psi(U$sqrtInv() %*% (y - x %*% beta))

    as.numeric(
      crossprod(psiResid, U$sqrt()) %*% matV$VInv() %*%
        matV$deriv[[derivSelect]]() %*%
        matV$VInv() %*% U$sqrt() %*% psiResid -
        matTrace(K * matV$VInv() %*% matV$deriv[[derivSelect]]())
    )
  }
}

fixedPointNumericDelta <- function(y, x, beta, matVFun, psi, K, derivSelect, stepSize, lowerBound) {
  obDelta <- robustObjectiveDelta(y, x, beta, matVFun, psi, K, derivSelect)
  function(rho) {
    rho + obDelta(rho) / ((obDelta(max(lowerBound, rho - stepSize)) - obDelta(rho)) / stepSize)
  }
}

fixedPointRobustDelta <- function(y, x, beta, matVFun, psi, K, derivSelect = 1) {
    # Precalculations - they only have to be done once
    mem1 <- (y - x %*% beta)

    function(param) {
        matV <- matVFun(param)
        U <- matU(matV$V())
        resid <- U$sqrtInv() %*% mem1
        psiResid <- psi(resid)
        c1 <- K / param * matTrace(matV$VInv() %*% matV$deriv[[derivSelect]]())
        c2 <- crossprod(psiResid, U$sqrt()) %*% matV$VInv() %*%
            matV$deriv[[derivSelect]]() %*% matV$VInv() %*% U$sqrt() %*% psiResid

        as.numeric(c2 / c1)
    }
}

fixedPointRobustDelta2 <- function(y, x, beta, matVFun, psi, K, derivSelect) {
  # Precalculations - they only have to be done once
  mem1 <- (y - x %*% beta)

  function(param) {
    matV <- matVFun(param)
    U <- matU(matV$V())
    resid <- U$sqrtInv() %*% mem1
    psiResid <- psi(resid)

    C1tmp <- K * matV$VInv() %*% matV$deriv[[derivSelect[1]]]() %*% matV$ZVuZInv()
    C2tmp <- K * matV$VInv() %*% matV$deriv[[derivSelect[2]]]() %*% matV$ZVuZInv()

    C1 <- matrix(ncol = 2, c(
      matTrace(C1tmp %*% matV$ZVuBarZ[[derivSelect[1]]]()),
      matTrace(C1tmp %*% matV$ZVuBarZ[[derivSelect[2]]]()),
      matTrace(C2tmp %*% matV$ZVuBarZ[[derivSelect[1]]]()),
      matTrace(C2tmp %*% matV$ZVuBarZ[[derivSelect[2]]]())
    ))

    c2 <- lapply(matV$deriv[derivSelect], function(deriv) {
      as.numeric(
        crossprod(psiResid, U$sqrt()) %*% matV$VInv() %*%
          deriv() %*% matV$VInv() %*% U$sqrt() %*% psiResid
      )}) %>% unlist

    as.numeric(solve(C1) %*% c2)
  }
}

fixedPointRobustRandomEffect <- function(y, x, beta, matV, psi) {

    makeMatB <- matBConst(y, x, beta, matV, psi)
    memResid <- y - x %*% beta
    function(u) as.numeric(as.matrix(makeMatB(u)) %*% memResid)
}
