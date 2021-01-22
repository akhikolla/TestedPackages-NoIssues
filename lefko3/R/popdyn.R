#' Estimate Mean Projection Matrices
#' 
#' \code{lmean()} estimates mean projection matrices as element-wise arithmetic
#' means.
#' 
#' @param mats A \code{lefkoMat} object.
#' @param matsout A string identifying which means to estimate. Option "pop"
#' indicates population-level only, "patch" indicates patch-level only, and
#' "all" indicates that both patch- and population-level means should be
#' estimated.
#' 
#' @return Yields a \code{lefkoMat} object with the following characteristics:
#' 
#' \item{A}{A list of full mean projection matrices in order of sorted populations,
#' patches, and years. These are typically estimated as the sums of the associated
#' mean \code{U} and \code{F} matrices. All matrices output in the \code{matrix} class.}
#' \item{U}{A list of mean survival-transition matrices sorted as in \code{A}. All 
#' matrices output in the \code{matrix} class.}
#' \item{F}{A list of mean fecundity matrices sorted as in \code{A}. All matrices 
#' output in the \code{matrix} class.}
#' \item{hstages}{A data frame showing the pairing of ahistorical stages used to
#' create historical stage pairs. Given if the MPM is historical.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages.}
#' \item{labels}{A data frame detailing the order of population, patch, and year 
#' of each mean matrix. If \code{pop}, \code{patch}, or \code{year2} are NA in the
#' original \code{labels} set, then these will be re-labeled as \code{A}, \code{1}, or \code{1},
#' respectively.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements in
#' \code{U} and \code{F} mean matrices, and the number of annual matrices.}
#' \item{modelqc}{The \code{qc} portion of the modelsuite input, if provided.}
#' \item{dataqc}{A vector showing the numbers of individuals and rows in the
#' vertical dataset used as input.}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' ehrlen3mean$A[[1]]
#' 
#' @export
lmean <- function(mats, matsout = "all") {
  
  if (class(mats) != "lefkoMat") {
    stop("An object of class lefkoMat is required as input.")
  }
  
  matsout.possible <- c("all", "pop", "patch")
  
  if (!is.element(tolower(matsout), matsout.possible)) {
    stop("The matsout option must take a value of all, pop, or patch.")
  }
  
  #First we create an index based on the input lefkoMat object
  listofyears <- mats$labels
  
  if (all(is.na(listofyears$pop))) {
    listofyears$pop <- 1
  }
  
  if (all(is.na(listofyears$patch))) {
    listofyears$patch <- 1
  }
  
  if (all(is.na(listofyears$year2))) {
    listofyears$year2 <- 1
  }
  
  listofyears$poppatch <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {
    paste(listofyears$pop[X], listofyears$patch[X])
  })
  
  listofyears$popc <- apply(as.matrix(listofyears$pop), 1, function(X) {which(unique(listofyears$pop) == X)-1})
  listofyears$poppatchc <- apply(as.matrix(listofyears$poppatch), 1, function(X) {which(unique(listofyears$poppatch) == X)-1})
  listofyears$year2c <- apply(as.matrix(listofyears$year2), 1, function(X) {which(unique(listofyears$year2) == X)-1})
  
  listofyears$patchesinpop <- apply(as.matrix(c(1:length(listofyears$poppatchc))), 1, function(X) {length(unique(listofyears$poppatchc[which(listofyears$popc == listofyears$popc[X])]))})
  listofyears$yearsinpatch <- apply(as.matrix(c(1:length(listofyears$year2c))), 1, function(X) {length(unique(listofyears$year2c[which(listofyears$poppatchc == listofyears$poppatchc[X])]))})
  
  numofpops <- length(unique(listofyears$popc))
  numofpatches <- length(unique(listofyears$poppatchc))
  numofyears <- length(unique(listofyears$year2c))
  
  listofyears$pop <- as.numeric(listofyears$pop)
  listofyears$patch <- listofyears$poppatchc + 1
  listofyears$year2 <- as.numeric(listofyears$year2)
  listofyears$poppatchc <- as.numeric(listofyears$poppatchc)
  
  if (matsout == "all") {
    poponly <- 1;
    patchonly <- 1;
  } else if (matsout == "patch") {
    poponly <- 0;
    patchonly <- 1;
  } else {
    poponly <- 1;
    patchonly <- 0;
  }
  
  if (length(unique(listofyears$poppatchc)) == 1) {
    poponly <- 0
    patchonly <- 1
  }
  
  if (is.element("modelqc", names(mats))) {
    extraqc <- mats$modelqc
  } else {
    extraqc <- mats$dataqc
  }
  
  if (!all(is.na(mats$hstages))) {
    output <- turbogeodiesel(listofyears, mats$U, mats$F, mats$ahstages, mats$hstages, extraqc, patchonly, poponly)
  } else {
    output <- geodiesel(listofyears, mats$U, mats$F, mats$ahstages, extraqc, patchonly, poponly)
    output$hstages <- NA
  }
  
  if (is.element("dataqc", names(mats))) {
    names(output)[which(names(output) == "modelqc")] <- "dataqc"
  }
  
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Dominant Eigenvalue and Deterministic Population Growth Rate Estimation
#' 
#' \code{lambda3()} is a generic function that returns the dominant eigenvalue of
#' a matrix, and set of dominant eigenvalues of a set of matrices. It can handle
#' very large and sparse matrices supplied as \code{lefkoMat} objects or as
#' individual matrices, and can be used with large historical matrices, IPMs, 
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats A lefkoMat object, or a single projection matrix, for which the
#' dominant eigenvalue is desired.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' See related functions for more details.
#' 
#' @seealso \code{\link{lambda3.lefkoMat}()}
#' @seealso \code{\link{lambda3.matrix}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' lambda3(ehrlen3mean)
#' 
#' @export
lambda3 <- function(mats) UseMethod("lambda3")

#' Estimate Deterministic Population Growth Rates of Matrices in a lefkoMat Object
#' 
#' \code{lambda3.lefkoMat()} returns the dominant eigenvalues of projection
#' matrices supplied within \code{lefkoMat} objects. This function can handle large 
#' and sparse matrices, and so can be used with large historical matrices, IPMs,
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' 
#' @return This function returns the dominant eigenvalue of each \code{$A} matrix in
#' the lefkoMat object input. For square matrices with fewer than 400 rows, this is
#' given as the largest real part of all eigenvalues estimated via the \code{eig_gen}()
#' function in the C++ Armadillo library. For larger matrices, the function assumes
#' that matrices are sparse and uses \code{eigs_gen}() instead. The output includes
#' a data frame showing the population, patch, and lambda estimate for each \code{$A}
#' matrix within the object. Row names correspond to the number of the matrix
#' within the \code{$A} element of the \code{lefkoMat} object.
#' 
#' @seealso \code{\link{lambda3}()}
#' @seealso \code{\link{lambda3.matrix}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' lambda3(ehrlen3mean)
#' 
#' @export
lambda3.lefkoMat <- function(mats) {
  baldrick <- if (any(class(mats$A) == "matrix")) {
    
    if (dim(mats$A)[1] > 400) {
      lambda3matrixsp(mats$A)
    } else {
      lambda3matrix(mats$A)
    }
    
  } else if (class(mats$A) == "list") {
    
    if (dim(mats$A[[1]])[1] > 400) {
      unlist(lapply(mats$A, lambda3matrixsp))
    } else {
      unlist(lapply(mats$A, lambda3matrix))
    }
    
  } else {
    
    stop("Input not recognized.")
    
  }
  
  output <- cbind.data.frame(mats$labels, baldrick)
  rownames(output) <- c(1:length(baldrick))
  
  names(output)[length(names(output))] <- "lambda"
  
  return(output)
}

#' Estimate Deterministic Population Growth Rate of a Projection Matrix
#' 
#' \code{lambda3.matrix()} returns the dominant eigenvalue of a single projection
#' matrix. This function can handle large and sparse matrices, so can be used with
#' large historical matrices, IPMs, age x stage matrices, as well as smaller
#' ahistorical matrices.
#' 
#' @param mats A population projection matrix of class \code{matrix}.
#'
#' @return This function returns the dominant eigenvalue of the matrix. For square
#' matrices with fewer than 400 rows, this is given as the largest real part of all
#' eigenvalues estimated via the \code{eig_gen}() function in package the C++ 
#' Armadillo library. For larger matrices, the matrix is assumed to be sparse and
#' \code{eigs_gen}() is used instead.
#' 
#' @seealso \code{\link{lambda3}()}
#' @seealso \code{\link{lambda3.lefkoMat}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' lambda3(ehrlen3mean$A[[1]])
#' 
#' @export
lambda3.matrix <- function(mats)
{
  if (dim(mats)[1] > 400) {
    lambda <- lambda3matrixsp(mats);
  } else {
    lambda <- lambda3matrix(mats);
  }
  
  return(lambda)
}

#' Stable Stage Distribution Estimation
#' 
#' \code{stablestage3()} is a generic function that returns the stable stage 
#' distribution for a population projection matrix or set of matrices. This
#' function is made to handle very large and sparse matrices supplied as 
#' \code{lefkoMat} objects or as individual matrices, and can be used with large
#' historical matrices, IPMs, age x stage matrices, as well as smaller ahistorical
#' matrices.
#' 
#' @param mats A lefkoMat object, or population projection matrix, for which the
#' stable stage distribution is desired.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' See related functions for details.
#' 
#' @seealso \code{\link{stablestage3.lefkoMat}()}
#' @seealso \code{\link{stablestage3.matrix}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' stablestage3(ehrlen3mean)
#' 
#' @export
stablestage3 <- function(mats) UseMethod("stablestage3")

#' Estimate Stable Stage Distribution for a lefkoMat Object
#' 
#' \code{stablestage3.lefkoMat()} returns the stable stage distributions of all
#' \code{$A} matrices in an object of class \code{lefkoMat}. This function can handle
#' large and sparse matrices, and so can be used with large historical matrices,
#' IPMs, age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' 
#' @return This function returns the stable stage distributions corresponding to
#' the matrices in a \code{lefkoMat} object. For square matrices with fewer than
#' 400 rows, the stable stage distribution is given as the right eigenvector
#' associated with largest real part of all eigenvalues estimated via the \code{eig_gen}()
#' function in the C++ Armadillo library divided by the sum of the associated right
#' eigenvector. For larger matrices, the function assumes that the matrix is sparse
#' and conducts a similar calculation but using the \code{eigs_gen}() for sparse matrix
#' eigen analysis.
#' 
#' The output depends on whether the \code{lefkoMat} object used as input is
#' ahistorical or historical. If the former, then a single data frame is output,
#' which includes the number of the matrix within the \code{$A} element of the input
#' \code{lefkoMat} object, followed by the stage id (numeric and assigned through
#' \code{\link{sf_create}()}), the stage name, and the estimated proportion of the
#' stable stage distribution (\code{ss_prop}).
#' 
#' If a historical matrix is used as input, then two data frames are output
#' into a list object. The \code{$hist} element contains a data frame where the 
#' stable stage distribution is given in terms of across-year stage pairs.
#' The structure includes the matrix number, the numeric stage designations for
#' stages in times \emph{t} and \emph{t}-1, respectively, followed by the
#' respective stage names, and ending with the estimated proportion of the stable
#' stage distribution for that stage within its matrix (\code{ss_prop}). The
#' \code{$ahist} element contains the stable stage distribution in stages
#' as given in the original stageframe. It includes a data frame with the matrix 
#' of origin, the numeric stage designation, stage name, and the stable stage
#' distribution estimated as the sum of distribution elements from \code{$hist}
#' corresponding to the equivalent stage in time \emph{t}, irrespective of stage in
#' time \emph{t}-1.
#' 
#' @seealso \code{\link{stablestage3}()}
#' @seealso \code{\link{stablestage3.matrix}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' stablestage3(ehrlen3mean)
#' 
#' @export
stablestage3.lefkoMat <- function(mats) {
  baldrick <- if (any(class(mats$A) == "matrix")) {
    
    if (dim(mats$A)[1] > 400) {
      ss3matrixsp(mats$A)
    } else {
      ss3matrix(mats$A)
    }
    
  } else if (class(mats$A) == "list") {
    
    if (dim(mats$A[[1]])[1] > 400) {
      unlist(lapply(mats$A, ss3matrixsp))
    } else {
      unlist(lapply(mats$A, ss3matrix))
    }
    
  } else {
    
    stop("Input not recognized.")
    
  }
  
  if (class(mats$A) == "list") {
    multiplier <- length(mats$A)
  } else multiplier <- 1
  
  if (all(is.na(mats$hstages))) {
    labels <- mats$ahstages[,1:2]
    
    modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
    
    output <- cbind.data.frame(modlabels, baldrick)
    names(output) <- c("matrix", "stage_id", "stage", "ss_prop")
    rownames(output) <- c(1:dim(output)[1])
  } else {
    labels <- mats$hstages
    
    modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
    
    outputh <- cbind.data.frame(modlabels, baldrick)
    names(outputh) <- c("matrix", "stage_id_2", "stage_id_1", "stage_2", "stage_1", "ss_prop")
    rownames(outputh) <- c(1:dim(outputh)[1])
    
    ahlabels <- mats$ahstages[,c("stage_id", "stage")]
    ss2 <- c(apply(as.matrix(c(1:multiplier)), 1, function(X) {
      rightset <- subset(outputh, matrix == X)
      apply(as.matrix(ahlabels[,1]), 1, function(Y) {
        sum(rightset$ss_prop[which(rightset$stage_id_2 == Y)])
      })
    }))
    outputah <- cbind.data.frame(rep(c(1:multiplier), each = length(ahlabels[,1])), rep(ahlabels[,1], multiplier), 
                                 rep(ahlabels[,2], multiplier), ss2)
    names(outputah) <- c("matrix", "stage_id", "stage", "ss_prop")
    rownames(outputah) <- c(1:dim(outputah)[1])
    
    output <- list(hist = outputh, ahist = outputah)
  }
  
  return(output)
}

#' Estimate Stable Stage Distribution for a Population Projection Matrix
#' 
#' \code{stablestage3.matrix()} returns the stable stage distribution for a 
#' population projection matrix. This function can handle large and sparse
#' matrices, and so can be used with large historical matrices, IPMs, 
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats A population projection matrix of class \code{matrix}.
#' 
#' @return This function returns the stable stage distribution corresponding to
#' the input matrix. For square matrices with fewer than 400 rows, the stable stage
#' distribution is given as the right eigenvector associated with largest real part
#' of the eigenvalues estimated for the matrix via the \code{eig_gen}() function in
#' the C++ Armadillo library, divided by the sum of the associated right
#' eigenvector. For larger matrices, the matrix is assumed to be sparse and the
#' calculation is conducted similarly but using \code{eig_gens}() instead.
#' 
#' @seealso \code{\link{stablestage3}()}
#' @seealso \code{\link{stablestage3.lefkoMat}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' stablestage3(ehrlen3mean$A[[1]])
#' 
#' @export
stablestage3.matrix <- function(mats)
{
  if (dim(mats)[1] > 400) {
    wcorr <- ss3matrixsp(mats)
  } else {
    wcorr <- ss3matrix(mats)
  }
  
  return(wcorr)
}

#' Reproductive Value Estimation
#' 
#' \code{repvalue3()} is a generic function that estimates returns the reproductive
#' values of stages in a population projection matrix or a set of matrices. The
#' specifics of estimation vary with the class of input object. This function is
#' made to handle very large and sparse matrices supplied as \code{lefkoMat} objects
#' or as individual matrices, and can be used with large historical matrices, IPMs,
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats A lefkoMat object, or population projection matrix.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' See related functions for details.
#' 
#' @seealso \code{\link{repvalue3.lefkoMat}()}
#' @seealso \code{\link{repvalue3.matrix}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ", reduce = TRUE)
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' repvalue3(ehrlen3mean)
#' 
#' @export
repvalue3 <- function(mats) UseMethod("repvalue3")

#' Estimate Reproductive Value for a lefkoMat Object
#' 
#' \code{repvalue3.lefkoMat()} returns the reproductive values for stages in a set
#' set of population projection matrices provided as a \code{lefkoMat} object. This
#' function can handle large and sparse matrices, and so can be used with large
#' historical matrices, IPMs, age x stage matrices, as well as smaller ahistorical
#' matrices.
#' 
#' @param mats An object of class \code{lefkoMat} object.
#' 
#' @return This function returns the reproductive values for stages of a matrices
#' within a \code{lefkoMat} object. The nature of the output depends on whether the
#' \code{lefkoMat} object used as input is ahistorical or historical. In both cases,
#' raw reproductive values are estimated as the left eigenvector associated with
#' the largest real part of the dominant eigenvalue, divided by the first non-zero
#' element of the left eigenvector. Eigen analysis is handled by the \code{eig_gen}()
#' function in the C++ Armadillo library for square matrices with fewer than 400
#' rows, and \code{eigs_gen}() for larger matrices.
#' 
#' If an ahistorical matrix set is used as input, then the output is a data frame
#' that includes the number of the matrix within the \code{$A} element of the input
#' \code{lefkoMat} object, followed by the numeric stage designation, the stage name,
#' and the reproductive value estimated (\code{repvalue}). 
#' 
#' If a historical matrix set is used as input, then a list with two elements is
#' output. The first element is a data frame showing the reproductive values given
#' given in terms of across-year stage pairs, as estimated in the procedure
#' described above. The order of variables in the data frame is: the matrix of
#' origin, the stage names for stages in times \emph{t} and \emph{t}-1 respectively,
#' the numeric stage designations corresponding to these stages, and the
#' reproductive value (\code{rep_value}). The second element is a data frame showing the
#' reproductive values of the basic stages in the associated stageframe. The 
#' reproductive values in this second data frame are estimated via the approach
#' developed in Ehrlen (2000), in which each ahistorical stage's reproductive 
#' value is the average of the RVs summed by stage at time \emph{t} weighted by 
#' the proportion of that stage pair within the historical stable stage
#' distribution associated with the matrix.
#' 
#' @seealso \code{\link{repvalue3}()}
#' @seealso \code{\link{repvalue3.matrix}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ", reduce = TRUE)
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' repvalue3(ehrlen3mean)
#' 
#' @export
repvalue3.lefkoMat <- function(mats) {
  
  baldrick <- if (any(class(mats$A) == "matrix")) {
    
    if (dim(mats$A)[1] > 400) {
      rv3matrixsp(mats$A)
    } else {
      rv3matrix(mats$A)
    }
    
  } else if (class(mats$A) == "list") {
    
    if (dim(mats$A[[1]])[1] > 400) {
      unlist(lapply(mats$A, rv3matrixsp))
    } else {
      unlist(lapply(mats$A, rv3matrix))
    }
    
  } else {
    
    stop("Input not recognized.")
    
  }
  
  if (class(mats$A) == "list") {
    multiplier <- length(mats$A)
  } else multiplier <- 1
  
  if (all(is.na(mats$hstages))) {
    labels <- mats$ahstages[,1:2]
    
    modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
    
    output <- cbind.data.frame(modlabels, baldrick)
    names(output) <- c("matrix", "stage_id", "stage", "rep_value")
    rownames(output) <- c(1:dim(output)[1])
    
  } else {
    ss3 <- stablestage3.lefkoMat(mats)
    rahist <- ss3$ahist
    rhist <-ss3$hist
    rhist$ss3sum <- apply(as.matrix(c(1:dim(rhist)[1])), 1, function(X) {
      rahist$ss_prop[intersect(which(rahist$stage_id == rhist$stage_id_2[X]), 
                               which(rahist$matrix == rhist$matrix[X]))]
    })
    rhist$sscorr <- rhist$ss_prop / rhist$ss3sum
    rhist$sscorr[which(is.na(rhist$sscorr))] <- 0
    rhist$rv3raw <- baldrick
    
    rhist$rep_value <- rhist$sscorr * rhist$rv3raw
    
    outputh <- rhist[,c("matrix", "stage_2", "stage_1", "stage_id_2", "stage_id_1", "rep_value")]
    
    ahlabels <- mats$ahstages[,c("stage_id", "stage")]
    rv2 <- Re(c(apply(as.matrix(c(1:multiplier)), 1, function(X) {
      rightset <- subset(outputh, matrix == X)
      apply(as.matrix(ahlabels[,1]), 1, function(Y) {
        sum(rightset$rep_value[which(rightset$stage_id_2 == Y)])
      })
    })))
    outputah <- cbind.data.frame(rep(c(1:multiplier), each = length(ahlabels[,1])), rep(ahlabels[,1], multiplier), 
                                 rep(ahlabels[,2], multiplier), rv2)
    names(outputah) <- c("matrix", "stage_id", "stage", "rep_value_unc")
    
    outputah$rep_value <- apply(as.matrix(c(1:dim(outputah)[1])), 1, function(X) {
      matsub <- subset(outputah, matrix == outputah$matrix[X])
      entrystage <- min(which(abs(matsub$rep_value_unc) > 0))
      return(outputah$rep_value_unc[X] / matsub$rep_value_unc[entrystage])
    })
    outputah <- outputah[,c("matrix", "stage_id", "stage", "rep_value")]
    rownames(outputah) <- c(1:dim(outputah)[1])
    
    output <-list(hist = outputh, ahist = outputah)
    
  }
  
  return(output)
}

#' Estimate Reproductive Value for a Population Projection Matrix
#' 
#' \code{repvalue3.matrix()} returns the reproductive values for stages in a 
#' population projection matrix. The function makes no assumptions about whether
#' the matrix is ahistorical and simply provides standard reproductive values
#' corresponding to each row, meaning that the overall reproductive values of basic
#' life history stages in a historical matrix are not provided (the 
#' \code{\link{repvalue3.lefkoMat}()} function estimates these on the basis of stage
#' description information provided in the \code{lefkoMat} object used as input in
#' that function). This function can handle large and sparse matrices, and so can
#' be used with large historical matrices, IPMs, age x stage matrices, as well as
#' smaller ahistorical matrices.
#' 
#' @param mats A population projection matrix.
#' 
#' @return This function returns a vector data frame characterizing the 
#' reproductive values for stages of a population projection matrix. This is 
#' given as the left eigenvector associated with largest real part of the
#' dominant eigenvalue, divided by the first non-zero element of the left 
#' eigenvector. Eigen analysis is handled by the \code{eig_gen}() function in the
#' C++ Armadillo library for square matrices with fewer than 400 rows, and using
#' the \code{eigs_gen}() function for larger matrices.
#' 
#' @seealso \code{\link{repvalue3}()}
#' @seealso \code{\link{repvalue3.lefkoMat}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ", reduce = TRUE)
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' repvalue3(ehrlen3mean$A[[1]])
#' 
#' @export
repvalue3.matrix <- function(mats)
{
  if (dim(mats)[1] > 400) {
    v <- rv3matrixsp(mats)
  } else {
    v <- rv3matrix(mats)
  }
  
  return(v)
}

#' Calculate Sensitivity of Lambda to Matrix Elements
#' 
#' \code{sensitivity3()} is a generic function that returns the sensitivity of
#' the deterministic population growth rate, lambda, to the elements of the
#' matrix population model. This function is made to handle very large and
#' sparse matrices supplied as \code{lefkoMat} objects or as individual
#' matrices.
#' 
#' @param mats A lefkoMat object, or population projection matrix, for which
#' the stable stage distribution is desired.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' 
#' @seealso \code{\link{sensitivity3.lefkoMat}()}
#' @seealso \code{\link{sensitivity3.matrix}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' sensitivity3(ehrlen3mean)
#' 
#' @export
sensitivity3 <- function(mats) UseMethod("sensitivity3")

#' Calculate Sensitivity of Lambda to Matrix Elements for a lefkoMat Object
#' 
#' \code{sensitivity3.lefkoMat()} returns the sensitivities of lambda to elements
#' of all \code{$A} matrices in an object of class \code{lefkoMat}. This function
#' can handle large and sparse matrices, and so can be used with large historical 
#' matrices, IPMs, age x stage matrices, as well as smaller ahistorical 
#' matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' 
#' @return The output from this function depends on whether the input is a
#' historical or ahistorical \code{lefkoMat} object. If the latter, then this function
#' returns a list with two elements, the first being a list of sensitivity matrices
#' corresponding to the \code{$A} matrices in a \code{lefkoMat} object, and the second being
#' a data frame showing the order of stages. If the former, then the output also
#' includes the equivalent sensitivities to ahistorical transitions estimated using
#' the historical matrix, and a data frame detailing the order of historical paired 
#' stages.
#' 
#' @seealso \code{\link{sensitivity3}()}
#' @seealso \code{\link{sensitivity3.matrix}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' sensitivity3(ehrlen3mean)
#' 
#' @export
sensitivity3.lefkoMat <- function(mats) {
  baldrick <- if (any(class(mats$A) == "matrix")) {
    
    if (dim(mats$A)[1] > 400) {
      sens3matrixsp(mats$A)
    } else {
      sens3matrix(mats$A)
    }
    
  } else if (class(mats$A) == "list") {
    
    if (all(is.na(mats$hstages))) {
      
      if (dim(mats$A[[1]])[1] > 400) {
        lapply(mats$A, sens3matrixsp)
      } else {
        lapply(mats$A, sens3matrix)
      }
      
    } else {
      
      lapply(mats$A, sens3hlefko, mats$ahstages, mats$hstages)
      
    }
    
  } else {
    
    stop("Input not recognized.")
    
  }
  
  if (all(is.na(mats$hstages))) {
    
    ahlabels <- mats$ahstages[,1:2]
    names(ahlabels) <- c("stage_id", "stage")
    
    output <- list(sensmats = baldrick, stages = ahlabels)
    
  } else {
    
    he_list <- lapply(baldrick, function(X) {X$h_smat})
    ahe_list <- lapply(baldrick, function(X) {X$ah_smat})
    
    hlabels <- mats$hstages
    
    ahlabels <- mats$ahstages[,1:2]
    names(ahlabels) <- c("stage_id", "stage")
    
    output <- list(h_sensmats = he_list, ah_sensmats = ahe_list, h_stages = hlabels, 
                   ah_stages = ahlabels)
    
  }
  
  return(output)
}

#' Calculate Sensitivity of Lambda to Matrix Elements for a Matrix
#' 
#' \code{sensitivity3.matrix()} returns the sensitivities of lambda to elements
#' of a single matrix. This function can handle large and sparse matrices, and 
#' so can be used with large historical matrices, IPMs, age x stage matrices,
#' as well as smaller ahistorical matrices.
#' 
#' @param mats An object of class \code{matrix}.
#' 
#' @return This function returns a single sensitivity matrix.
#' 
#' @seealso \code{\link{sensitivity3}()}
#' @seealso \code{\link{sensitivity3.lefkoMat}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' sensitivity3(ehrlen3mean$A[[1]])
#' 
#' @export
sensitivity3.matrix <- function(mats)
{
  if (dim(mats)[1] > 400) {
    wcorr <- sens3matrixsp(mats)
  } else {
    wcorr <- sens3matrix(mats)
  }
  
  return(wcorr)
}

#' Calculate Elasticity of Lambda to Matrix Elements
#' 
#' \code{elasticity3()} is a generic function that returns the elasticity of
#' the deterministic population growth rate, lambda, to the elements of the
#' matrix population model. This function is made to handle very large and
#' sparse matrices supplied as \code{lefkoMat} objects or as individual
#' matrices.
#' 
#' @param mats A lefkoMat object, or population projection matrix, for which
#' the stable stage distribution is desired.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' 
#' @seealso \code{\link{elasticity3.lefkoMat}()}
#' @seealso \code{\link{elasticity3.matrix}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' elasticity3(ehrlen3mean)
#' 
#' @export
elasticity3 <- function(mats) UseMethod("elasticity3")

#' Calculate Elasticity of Lambda to Matrix Elements for a lefkoMat Object
#' 
#' \code{elasticity3.lefkoMat()} returns the elasticities of lambda to elements
#' of all \code{$A} matrices in an object of class \code{lefkoMat}. This function
#' can handle large and sparse matrices, and so can be used with large historical 
#' matrices, IPMs, age x stage matrices, as well as smaller ahistorical 
#' matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' 
#' @return This function's output depends on whether the lefkoMat object is historical
#' or ahistorical. If the former, then it is a list with four elements. The first 
#' (\code{h_elasmats}) is a list of historical elasticity matrices, the second 
#' (\code{ah_elasmats}) is a list of sensitivity matrices in which historical
#' elasticities have been summed by the stage in times \emph{t} and \emph{t}+1, to produce 
#' elasticity matrices equivalent in principle to ahistorical elasticity matrices but 
#' reflecting the effects of stage in time \emph{t}-1. The third element (\code{h_stages}) is a 
#' data frame showing historical stage pairs, and the fourth (\code{ah_stages}) is a data 
#' frame showing the ahistorical stages.
#' 
#' @seealso \code{\link{elasticity3}()}
#' @seealso \code{\link{elasticity3.matrix}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' elasticity3(ehrlen3mean)
#' 
#' @export
elasticity3.lefkoMat <- function(mats) {
  baldrick <- if (any(class(mats$A) == "matrix")) {
    
    if (dim(mats$A)[1] > 400) {
      elas3matrixsp(mats$A)
    } else {
      elas3matrix(mats$A)
    }
    
  } else if (class(mats$A) == "list") {
    
    if (all(is.na(mats$hstages))) {
      
      if (dim(mats$A[[1]])[1] > 400) {
        lapply(mats$A, elas3matrixsp)
      } else {
        lapply(mats$A, elas3matrix)
      }
      
    } else {
      
      lapply(mats$A, elas3hlefko, mats$ahstages, mats$hstages)
      
    }
    
  } else {
    
    stop("Input not recognized.")
    
  }
  
  if (class(mats$A) == "list") {
    multiplier <- length(mats$A)
  } else multiplier <- 1
  
  if (all(is.na(mats$hstages))) {
    
    ahlabels <- mats$ahstages[,1:2]
    names(ahlabels) <- c("stage_id", "stage")
    
    output <- list(elasmats = baldrick, stages = ahlabels)
    
  } else {
    
    he_list <- lapply(baldrick, function(X) {X$h_emat})
    ahe_list <- lapply(baldrick, function(X) {X$ah_emat})
    
    hlabels <- mats$hstages
    
    ahlabels <- mats$ahstages[,1:2]
    names(ahlabels) <- c("stage_id", "stage")
    
    output <- list(h_elasmats = he_list, ah_elasmats = ahe_list, h_stages = hlabels, 
                   ah_stages = ahlabels)
    
  }
  
  return(output)
}

#' Calculate Elasticity of Lambda to Matrix Elements for a Matrix
#' 
#' \code{elasticity3.matrix()} returns the elasticities of lambda to elements
#' of a single matrix. This function can handle large and sparse matrices, and 
#' so can be used with large historical matrices, IPMs, age x stage matrices,
#' as well as smaller ahistorical matrices.
#' 
#' @param mats An object of class \code{matrix}.
#' 
#' @return This function returns a single elasticity matrix.
#' 
#' @seealso \code{\link{elasticity3}()}
#' @seealso \code{\link{elasticity3.lefkoMat}()}
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          sizeacol = "Volume88", repstracol = "FCODE88", 
#'                          fecacol = "Intactseed88", deadacol = "Dead1988", 
#'                          nonobsacol = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' elasticity3(ehrlen3mean$A[[1]])
#' 
#' @export
elasticity3.matrix <- function(mats)
{
  if (dim(mats)[1] > 400) {
    wcorr <- elas3matrixsp(mats)
  } else {
    wcorr <- elas3matrix(mats)
  }
  
  return(wcorr)
}

