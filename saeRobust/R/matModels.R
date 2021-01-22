matVFH <- function(.sigma2, .samplingVar) {

  .diag <- function(x) Diagonal(x = x)

  Vu <- getter(rep(.sigma2, length(.samplingVar)), .diag)
  VuInv <- getter(1 / rep(.sigma2, length(.samplingVar)), .diag)
  Ve <- getter(.samplingVar, .diag)
  VeInv <- getter(1 / .samplingVar, .diag)
  V <- getter(Vu() + Ve())
  VInv <- getter(solve(V()))
  Z <- getter(Diagonal(length(.samplingVar)))

  deriv <- list(
    getter(Diagonal(length(.samplingVar)))
  )

  retList()

}

matVSFH <- function(.rho, .sigma2, .W, .samplingVar) {

  .diag <- function(x) Diagonal(x = x)

  Ve <- getter(.samplingVar, .diag)
  VeInv <- getter(1 / .samplingVar, .diag)

  Omega1 <- getter(matOmega1(W = .W, rho = .rho), Matrix)
  Vu <- getter(.sigma2 * Omega1())
  VuInv <- getter(solve(Vu()))

  V <- getter(Vu() + Ve(), as.matrix)
  VInv <- getter(solve(V()), as.matrix)

  Z <- getter(Diagonal(length(.samplingVar)))

  deriv <- list(
    rho = getter(matVDerR1(
      .rho, .sigma2, as.matrix(Z()), as.matrix(Omega1()), as.matrix(.W))),
    sigma2 = getter(as.matrix(Omega1()))
  )

  retList()

}

matVTFH <- function(.rho, .sigma2, .nTime, .samplingVar) {

  .diag <- function(x) Diagonal(x = x)
  .emptyMatrix <- function(dim) Matrix(0, dim, dim)
  .nDomains <- length(.samplingVar) / .nTime


  Ve <- getter(.samplingVar, .diag)
  VeInv <- getter(1 / .samplingVar, .diag)

  Omega1 <- getter(Diagonal(.nDomains))
  Omega2 <- getter(matOmega2(.nTime, .rho))

  .Omega1Inv <- getter(solve(Omega1()))
  .Omega2Inv <- getter(solve(Omega2()))

  Z <- getter(matTZ(.nDomains, .nTime))
  Z1 <- getter(matTZ1(.nDomains, .nTime))

  Vu <- getter(
    bdiag(.sigma2[1] * Omega1(),
          .sigma2[2] * matBlockDiagonal(Omega2(), .nDomains))
  )

  VuInv <- getter(
    bdiag(.Omega1Inv() / .sigma2[1],
          matBlockDiagonal(.Omega2Inv(), .nDomains) / .sigma2[2])
  )

  .V <- getter(matVInvT(
    as.matrix(Omega1()),
    .sigma2[1],
    .rho, .sigma2[2],
    as.matrix(Z1()),
    .samplingVar
  ))

  V <- getter(.V()$V)
  VInv <- getter(.V()$VInv)

  deriv <- list(
    sigma21 = getter(matVDerS1(as.matrix(Omega1()), as.matrix(Z1()))),
    rho = getter(matVDerR2(.rho, .sigma2[2], Omega2(), .nDomains)),
    sigma22 = getter(matVDerS2(Omega2(), .nDomains))
  )

  .ZVuZ <- getter(matVInvT(
    as.matrix(Omega1()),
    .sigma2[1],
    .rho, .sigma2[2],
    as.matrix(Z1()),
    rep(0, length(.samplingVar))
  ))

  ZVuZInv <- getter(.ZVuZ()$VInv)

  ZVuBarZ <- list(
    sigma21 = getter(Z() %*% tcrossprod(
      bdiag(Omega1(), .emptyMatrix(.nDomains * .nTime)), Z())),
    sigma22 = getter(Z() %*% tcrossprod(
      bdiag(.emptyMatrix(.nDomains), matBlockDiagonal(Omega2(), .nDomains)), Z()))
  )

  retList()

}

matVSTFH <- function(.rho, .sigma2, .W, .nTime, .samplingVar) {

  .diag <- function(x) Diagonal(x = x)
  .emptyMatrix <- function(dim) Matrix(0, dim, dim)
  .nDomains <- length(.samplingVar) / .nTime


  Ve <- getter(.samplingVar, .diag)
  VeInv <- getter(1 / .samplingVar, .diag)

  Omega1 <- getter(matOmega1(W = .W, rho = .rho[1]), Matrix)
  Omega2 <- getter(matOmega2(.nTime, .rho[2]))

  .Omega1Inv <- getter(solve(Omega1()))
  .Omega2Inv <- getter(solve(Omega2()))

  Z <- getter(matTZ(.nDomains, .nTime))
  Z1 <- getter(matTZ1(.nDomains, .nTime))

  Vu <- getter(
    bdiag(.sigma2[1] * Omega1(),
          .sigma2[2] * matBlockDiagonal(Omega2(), .nDomains))
  )

  VuInv <- getter(
    bdiag(.Omega1Inv() / .sigma2[1],
          matBlockDiagonal(.Omega2Inv(), .nDomains) / .sigma2[2])
  )

  .V <- getter(matVInvT(
    as.matrix(Omega1()),
    .sigma2[1],
    .rho[2], .sigma2[2],
    as.matrix(Z1()),
    .samplingVar
  ))

  V <- getter(.V()$V)
  VInv <- getter(.V()$VInv)

  deriv <- list(
    rho1 = getter(matVDerR1(
      .rho[1], .sigma2[1], as.matrix(Z1()), as.matrix(Omega1()), as.matrix(.W))),
    rho2 = getter(matVDerR2(.rho[2], .sigma2[2], Omega2(), .nDomains)),
    sigma21 = getter(matVDerS1(as.matrix(Omega1()), as.matrix(Z1()))),
    sigma22 = getter(matVDerS2(Omega2(), .nDomains))
  )

  .ZVuZ <- getter(matVInvT(
    as.matrix(Omega1()),
    .sigma2[1],
    .rho[2], .sigma2[2],
    as.matrix(Z1()),
    rep(0, length(.samplingVar))
  ))

  ZVuZInv <- getter(.ZVuZ()$VInv)

  ZVuBarZ <- list(
    sigma21 = getter(Z() %*% tcrossprod(
      bdiag(Omega1(), .emptyMatrix(.nDomains * .nTime)), Z())),
    sigma22 = getter(Z() %*% tcrossprod(
      bdiag(.emptyMatrix(.nDomains), matBlockDiagonal(Omega2(), .nDomains)), Z()))
  )

  retList()

}
