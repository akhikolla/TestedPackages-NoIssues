
.pkgenv <- new.env(parent=emptyenv())

.onLoad <- function(libname, pkgname) {
    LIB_GSL <- Sys.getenv("LIB_GSL")
    .pkgenv[["zero"]] <- vlivC(1, c(0))
    .pkgenv[["one"]] <- vlivC(1, c(1))
    .pkgenv[["two"]] <- vlivC(1, c(2))
    .pkgenv[["three"]] <- vlivC(1, c(3))
    .pkgenv[["four"]] <- vlivC(1, c(4))
    .pkgenv[["five"]] <- vlivC(1, c(5))
    .pkgenv[["six"]] <- vlivC(1, c(6))
    .pkgenv[["seven"]] <- vlivC(1, c(7))
    .pkgenv[["eight"]] <- vlivC(1, c(8))
}
